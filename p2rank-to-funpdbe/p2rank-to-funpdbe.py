#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import csv
import datetime
import json
import sys
import shutil
import multiprocessing
import argparse
import enum

from validator.validator import Validator
from validator.residue_index import ResidueIndexes

DATA_RESOURCE = "p2rank"

RESOURCE_VERSION = "1.0"

P2RANK_VERSION = "2.0"

P2RANK_WEB_URL = "http://prankweb.cz/analyze/id_noconser/{}"

RELEASE_DATE = datetime.date.today().strftime("%d/%m/%Y")

EVIDENCE_CODE_ONTOLOGY = [
    {
        "eco_term": "computational combinatorial evidence",
        "eco_code": "ECO_0000246"
    }
]

PROBABILITY_TO_CONFIDENCE_THRESHOLDS = [
    (0.33, "low"),
    (0.6, "medium"),
    (1.0, "high"),
]

THIS_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

FUNPDBE_SCHEMA = os.path.join(
    THIS_DIRECTORY, "..", "funpdbe-validator", "data", "funpdbe_schema.json")


class ConversionStatus(enum.Enum):
    DONE = 0
    SKIPPED = 1
    FAILED = 2


def main(arguments):
    init_logging()
    initialize_directories(arguments)
    check_directories()
    pdb_ids = collect_pdb_ids(arguments["input"])
    if arguments["threads"] < 2:
        convert_single_thread(pdb_ids, arguments)
    else:
        convert_multiple_threads(pdb_ids, arguments)


def init_logging(level=logging.DEBUG):
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%H:%M:%S")


def initialize_directories(arguments):
    os.makedirs(arguments["input"], exist_ok=True)
    os.makedirs(arguments["output"], exist_ok=True)
    os.makedirs(arguments["errorOutput"], exist_ok=True)


def check_directories():
    if not os.path.exists(FUNPDBE_SCHEMA):
        raise Exception(
            "Missing FunPDBe schema. Please read installation guide.")


def collect_pdb_ids(input_dir):
    result = set()
    for file in os.listdir(input_dir):
        if "pdb" in file:
            result.add(file[:file.index(".")])
    return sorted(list(result))


def convert_single_thread(pdb_ids, arguments):
    input_dir = arguments["input"]
    output_dir = arguments["output"]
    error_dir = arguments["errorOutput"]

    failed = []
    skipped = []
    for index, pdb_id in enumerate(pdb_ids):
        logging.info("%i/%i : %s", index, len(pdb_ids), pdb_id)
        output_path = os.path.join(output_dir, pdb_id + ".json")
        try:
            result = convert_file(
                pdb_id.upper(),
                os.path.join(input_dir, pdb_id + ".pdb.gz_predictions.csv"),
                os.path.join(input_dir, pdb_id + ".pdb.gz_residues.csv"),
                output_path)
            if result == ConversionStatus.SKIPPED:
                skipped.append(pdb_id)
        except Exception as error:
            failed.append(pdb_id)
            logging.exception("Conversion failed for %s : %s", pdb_id, error)
            if os.path.exists(output_path):
                error_path = os.path.join(
                    error_dir, os.path.basename(output_path))
                shutil.move(output_path, error_path)
    logging.info("Converted: %s skipped: %s failed: %s",
                 len(pdb_ids) - len(failed) - len(skipped),
                 len(skipped), len(failed))


def convert_multiple_threads(pdb_ids, arguments):
    pool = multiprocessing.Pool(arguments["threads"])
    logging.info("Converting files ...")
    tasks = [{
        "pdb": pdb,
        "input": arguments["input"],
        "output": arguments["output"],
        "errorOutput": arguments["errorOutput"]
    } for pdb in pdb_ids]
    results = pool.map(convert_pdb_file, tasks)
    converted_count = len([x for x in results if x == ConversionStatus.DONE])
    skipped_count = len([x for x in results if x == ConversionStatus.SKIPPED])
    failed_count = len([x for x in results if x == ConversionStatus.FAILED])
    logging.info("Converting files ... done")
    logging.info("Converted: %s skipped: %s failed: %s",
                 converted_count, skipped_count, failed_count)


def convert_pdb_file(task):
    pdb_id = task["pdb"]
    input_dir = task["input"]
    output_dir = task["output"]
    error_dir = task["errorOutput"]
    #
    output_path = get_output_path(output_dir, pdb_id)
    try:
        return convert_file(
            pdb_id.upper(),
            os.path.join(input_dir, pdb_id + ".pdb.gz_predictions.csv"),
            os.path.join(input_dir, pdb_id + ".pdb.gz_residues.csv"),
            output_path)
    except Exception as error:
        logging.exception("Conversion failed for %s : %s", pdb_id, error)
        if os.path.exists(output_path):
            error_path = os.path.join(
                error_dir, os.path.basename(output_path))
            shutil.move(output_path, error_path)
        return ConversionStatus.FAILED


def get_output_path(output_dir, pdb_id):
    sub_dir = pdb_id[1:-1]
    return os.path.join(output_dir, sub_dir, pdb_id + ".json")


def convert_file(pdb_id, pocket_path, residues_path, output_path) \
        -> ConversionStatus:
    residues = read_residues(residues_path)
    pockets = read_pockets(pocket_path)

    if len(pockets) == 0:
        logging.debug("Skipping file with no predictions: %s", pdb_id)
        return ConversionStatus.SKIPPED

    sites = []
    chains_dictionary = {}

    for pocket in pockets:
        site = create_site(pocket)
        sites.append(site)
        for residue in iterate_site_residues(pocket, residues, site["site_id"]):
            add_residue_to_chains(residue, chains_dictionary)
    chains = flat_chains_dictionary(chains_dictionary)

    funpdbe = create_output_file(pdb_id, sites, chains)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as out_stream:
        json.dump(funpdbe, out_stream, indent=2)

    validate_file(output_path)
    return ConversionStatus.DONE


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
        "rank": row["rank"],
        "center_x": row["center_x"],
        "center_y": row["center_y"],
        "center_z": row["center_z"]
    } for row in iterate_csv_file(pocket_path)]


def create_site(pocket):
    site = {
        "site_id": int(pocket["name"].replace("pocket", "")),
        "label": pocket["name"],
        "additional_site_annotations": {
            "score": pocket["score"],
            "center": {
                "x": pocket["center_x"],
                "y": pocket["center_y"],
                "z": pocket["center_z"]
            }
        }
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--input", type=str,
                        default=os.path.join(
                            THIS_DIRECTORY, "..", "data", "p2rank-outputs"))
    parser.add_argument("--output", type=str,
                        default=os.path.join(
                            THIS_DIRECTORY, "..", "data", "funpdbe"))
    parser.add_argument("--errorOutput", type=str,
                        default=os.path.join(
                            THIS_DIRECTORY, "..", "data", "error"))
    main(vars(parser.parse_args()))
