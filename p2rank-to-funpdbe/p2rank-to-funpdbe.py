#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import csv
import datetime
import json
import sys
import shutil

from validator.validator import Validator
from validator.residue_index import ResidueIndexes

DATA_RESOURCE = "p2rank"

RESOURCE_VERSION = "1.0"

P2RANK_VERSION = "2.0"

P2RANK_WEB_URL = "http://prank.projekty.ms.mff.cuni.cz/analyze/id/{}"

RELEASE_DATE = datetime.date.today().strftime("%d/%m/%Y")

EVIDENCE_CODE_ONTOLOGY = [
    {  # TODO Decide what to use,
        "eco_term": "computational combinatorial evidence",
        "eco_code": "ECO_0000246"
    }
]

PROBABILITY_TO_CONFIDENCE_THRESHOLDS = [
    (0.5, "low"),
    (0.8, "medium"),
    (1.0, "high"),
]

THIS_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

INPUT_DIRECTORY = os.path.join(
    THIS_DIRECTORY, "..", "data", "p2rank-outputs")

OUTPUT_DIRECTORY = os.path.join(
    THIS_DIRECTORY, "..", "data", "funpdbe")

FUNPDBE_SCHEMA = os.path.join(
    THIS_DIRECTORY, "..", "funpdbe-validator", "data", "funpdbe_schema.json")

ERROR_DIRECTORY = os.path.join(
    THIS_DIRECTORY, "..", "data", "error")


def main():
    init_logging()
    initialize_directories()
    check_directories()

    pdb_ids = collect_pdb_ids()
    failed = []
    for index, pdb_id in enumerate(pdb_ids):
        logging.info("%i/%i : %s", index, len(pdb_ids), pdb_id)
        output_path = os.path.join(OUTPUT_DIRECTORY, pdb_id + ".json")
        try:
            convert_file(
                pdb_id.upper(),
                os.path.join(INPUT_DIRECTORY, pdb_id + ".pdb_predictions.csv"),
                os.path.join(INPUT_DIRECTORY, pdb_id + ".pdb_residues.csv"),
                output_path)
        except Exception as error:
            failed.append(pdb_id)
            logging.exception("Conversion failed for %s : %s", pdb_id, error)
            if os.path.exists(output_path):
                error_path = os.path.join(
                    ERROR_DIRECTORY, os.path.basename(output_path))
                shutil.move(output_path, error_path)
    logging.info("Converted: %s failed: %s",
                 len(pdb_ids) - len(failed), len(failed))


def init_logging(level=logging.DEBUG):
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%H:%M:%S")


def initialize_directories():
    os.makedirs(INPUT_DIRECTORY, exist_ok=True)
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)
    os.makedirs(ERROR_DIRECTORY, exist_ok=True)


def check_directories():
    if not os.path.exists(FUNPDBE_SCHEMA):
        raise Exception(
            "Missing FunPDBe schema. Please read installation guide.")


def collect_pdb_ids():
    result = set()
    for file in os.listdir(INPUT_DIRECTORY):
        if "pdb_predictions.csv" in file:
            result.add(file[:file.index(".")])
    return sorted(list(result))


def convert_file(pdb_id, pocket_path, residues_path, output_path):
    residues = read_residues(residues_path)
    pockets = read_pockets(pocket_path)

    sites = []
    chains_dictionary = {}

    for pocket in pockets:
        site = create_site(pocket)
        sites.append(site)
        for residue in iterate_site_residues(pocket, residues, site["site_id"]):
            add_residue_to_chains(residue, chains_dictionary)
    chains = flat_chains_dictionary(chains_dictionary)

    funpdbe = create_output_file(pdb_id, sites, chains)
    with open(output_path, "w") as out_stream:
        json.dump(funpdbe, out_stream, indent=2)

    validate_file(output_path)


def read_residues(residues_path):
    # chain, residue_label, residue_name, score, zscore, probability, pocket
    return [{
        "chains": row["chain"],
        "type": row["residue_name"],
        "label": row["residue_label"],
        "score": float(row["score"]),
        "probability": float(row["probability"]),
        "pocket": row["pocket"]
    } for row in iterate_csv_file(residues_path)]


def iterate_csv_file(path):
    with open(path) as input_stream:
        csv_reader = csv.reader(
            input_stream,
            delimiter=",",
            skipinitialspace=True)
        header = next(csv_reader)
        for row in csv_reader:
            yield {key: value for key, value in zip(header, row)}


def read_pockets(pocket_path):
    # name, rank, score, sas_points, surf_atoms,
    # center_x, center_y, center_z,
    # residue_ids, surf_atom_ids
    return [{
        "name": row["name"],
        "score": row["score"],
        "rank": row["rank"]
    } for row in iterate_csv_file(pocket_path)]


def create_site(pocket):
    site = {
        "site_id": int(pocket["name"].replace("pocket", "")),
        "label": pocket["name"]
    }
    return site


def iterate_site_residues(pocket, residues, site_id):
    for residue in residues:
        if residue["pocket"] == pocket["rank"]:
            yield residue_to_site_data(residue, site_id)


def residue_to_site_data(residue, site_id):
    confidence_class = probability_to_confidence_class(residue["probability"])
    return {
        "chain": residue["chains"],
        "res": residue["label"],
        "aa": residue["type"],
        "site_data": {
            "site_id_ref": site_id,
            "confidence_score": residue["probability"],
            "confidence_classification": confidence_class,
            "raw_score": residue["score"],
        }
    }


def probability_to_confidence_class(p):
    for threshold, value in PROBABILITY_TO_CONFIDENCE_THRESHOLDS:
        if p <= threshold:
            return value
    raise RuntimeError("Unexpected probability: {}".format(p))


def add_residue_to_chains(residue, chains_dictionary):
    if residue["chain"] not in chains_dictionary:
        chains_dictionary[residue["chain"]] = {}
    chain_residues = chains_dictionary[residue["chain"]]

    residue_key = (residue["res"], residue["aa"])
    if residue_key not in chain_residues:
        chain_residues[residue_key] = {
            "pdb_res_label": residue["res"],
            "aa_type": residue["aa"],
            "site_data": []
        }
        chain_residues[residue_key]["site_data"].append(residue["site_data"])


def flat_chains_dictionary(chains_dictionary):
    return [{
        "chain_label": chain_label,
        "residues": list(residues.values())
    } for chain_label, residues in chains_dictionary.items()]


def create_output_file(pdb_id, sites, chains):
    return {
        "data_resource": DATA_RESOURCE,
        "resource_version": RESOURCE_VERSION,
        "software_version": P2RANK_VERSION,
        "resource_entry_url": P2RANK_WEB_URL.format(pdb_id.upper()),
        "release_date": RELEASE_DATE,
        "pdb_id": pdb_id.lower(),
        "chains": chains,
        "sites": sites,
        "evidence_code_ontology": EVIDENCE_CODE_ONTOLOGY
    }


def validate_file(path):
    logging.debug("Starting validation.")
    validator = Validator(DATA_RESOURCE)
    validator.load_schema(FUNPDBE_SCHEMA)
    validator.load_json(path)
    if not validator.basic_checks():
        raise RuntimeError("Basic checks failed for {}".format(path))
    if not validator.validate_against_schema():
        logging.error(validator.error_log)
        raise RuntimeError("Invalid schema for {}".format(path))

    residue_indexes = ResidueIndexes(validator.json_data)
    if not residue_indexes.check_every_residue():
        logging.error(residue_indexes.mismatches)
        raise RuntimeError("Invalid residues: {}".format(path))

    logging.debug("Output file is valid.")


if __name__ == "__main__":
    main()
