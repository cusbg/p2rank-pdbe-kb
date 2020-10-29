# P2Rank to PDBe-KB export

This repository holds scripts which are used to export P2Rank predictions to PDBe-KB.

## PDBe-KB
More information about PDBe-KB can be found at:
* https://gitlab.ebi.ac.uk/pdbe-kb

You can register into the [PDBe-KB mailing list](https://listserver.ebi.ac.uk/mailman/listinfo/pdbe-kb-depositor) 
for email notifications.

## Quick start

### Prerequisites
* Python 3
* Git

### Installing
```
git clone https://github.com/skodapetr/p2rank-funpdbe.git
cd p2rank-funpdbe
git clone https://gitlabci.ebi.ac.uk/pdbe-kb/funpdbe/funpdbe-client.git
git clone https://gitlabci.ebi.ac.uk/pdbe-kb/funpdbe/funpdbe-schema.git
git clone https://gitlabci.ebi.ac.uk/pdbe-kb/funpdbe/funpdbe-validator.git
```
```
pip install jsonschema
pip install requests
```

### Execution
* Copy input files into ```/data/p2rank-outputs```.
* Setup python paths:
  * Windows: ```set PYTHONPATH=%cd%\funpdbe-validator\;%PYTHONPATH%```
  * Others: ```export PYTHONPATH="$PWD/funpdbe-validator/:$PYTHONPATH"```
* Run the conversion script:
```
python p2rank-to-funpdbe/p2rank-to-funpdbe.py
python p2rank-to-funpdbe/p2rank-to-funpdbe.py --threads 4
```
