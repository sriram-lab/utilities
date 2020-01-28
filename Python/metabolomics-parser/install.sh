#!/bin/bash

# Make virtualenv
virtualenv metabolomics-parser --Python=python3 --user
source metabolomics-parser/bin/activate

# Install necessary packages if they are not installed already
pip3 install -U -r requirements.txt --user

echo "Finished making virtual environment and installing necessary depedencies."
echo "To run the metabolomics parser, run './metabolomicsParser.sh [pathToMetabolomicsFile] [pathToMetabolicModel]' in the terminal."
