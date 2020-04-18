# metabolomics-parser
**Created by:** Scott Campit

The `metabolomics-parser` function maps metabolites from a metabolomics dataset to the metabolite positions in a genome-scale metabolic reconstruction. 

## Project Contents
```bash
.
|__ test
|     |__ testMetaboAnalystAPI.py 
|     |__ testPubChemAPI.py
|__ __init__.py     
|__ install.sh
|__ metabolomics-parser.py
|__ metabolomics-parser.sh
|__ README.md
|__ requirements.txt
```

## Installing `metabolomics-parser`
To install `metabolomics-parser`, run the `install.sh` file in the terminal. It will automatically create a new
 Python virtual environment and install the necessary packages to run the parser. 
 
 This package was tested and built using Python 3.7.

## Usage
### Linux
To use the `metabolomics-parser.sh` function in the Linux terminal, the following syntax can be used to run the
 function:
```bash
./metabolomics-parser.sh [path-to-metabolic-model.xml] [path-to-metabolomics-data] [OPTIONAL: excel file sheet name]
```

### Python
The main functions are written in the `metabolomics-parser.py` script. The required libraries are documented in the
 `requirement.txt` file. 

## Recent updates
**April 17, 2020: Scott Campit**
  * `metabolomics-parser` now takes in Excel files and can read in specific sheet names as the third argument.
  * This pipeline now synonym mapping for common metabolite names using the PubChem API to improve metabolite name
   coverage when mapping to the metabolic model. 
     * Note that this additional mapping step increases the mapping time significantly.

## Additional notes: 
  * The metabolic model must be an extensible markup file (.xml or .sbml).
  * The parser has only been tested using the RECON1 model, which navigates `xml` namespaces. 
    * Other metabolic models
   may have different `xml` namespaces, which you would have to change in the `metabolomics-parser.py` file.
