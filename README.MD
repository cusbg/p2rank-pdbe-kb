# P2Rank to FunPDBe conversion

## FunPDBe
More information about FunPDBe can be found at:
* https://www.ebi.ac.uk/pdbe/funpdbe/deposition
* https://gitlab.ebi.ac.uk/pdbe-kb

You can register into [mailing list](https://listserver.ebi.ac.uk/mailman/listinfo/pdbe-kb-depositor) 
for email notifications at:

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
  * Others: ```export PYTHONPATH=$PWD/funpdbe-validator:$PYTHONPATH"```
* Run the conversion script:
```
python p2rank-to-funpdbe/p2rank-to-funpdbe.py
python p2rank-to-funpdbe/p2rank-to-funpdbe.py --threads 4
```
