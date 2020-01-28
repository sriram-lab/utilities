# metabolomics-parser

The `metabolomics-parser` function maps metabolites from a metabolomics dataset to the metabolite positions in a genome-scale metabolic reconstruction. 

## Installing `metabolomics-parser`
To install `metabolomics-parser`, run the `install.sh` file in the terminal. It will automatically create a new Python virtual environment and install the necessary packages to run the parser. This package was tested and built on Python 3.7.

## Usage
To use the `metabolomics-parser` function, type the following command into the terminal:
```bash
./metabolomics-parser.sh [path-to-metabolomics-data.csv] [path-to-metabolic-model.xml]
```

Note that the metabolomics data currently must be as comma-delimited file (.csv), and the metabolic model must be an extensible markup file (.xml or .sbml).

Currently, the parser has been tested using RECON1. Other metabolic models may have different xml namespaces, which you would have to change in the `metabolomics-parser.py` file.
