"""
testPubChemAPI.py runs tests on different components of the metabolomics-parse.py
@author: Scott Campit
"""

import pubchempy as pcp
import pandas as pd

def test_parse_pubchem_compounds(filename):
    """
    :param:  filename:      A string denoting the path to the metabolomics file
    :return: all_compounds: A Pandas dataframe containing the PubMed mappings
    """
    fileData = pd.read_excel(filename)
    fileData = fileData.drop('Compound Method', axis=1).join(fileData['Compound Method']
                                                  .str.split('/', expand=True)
                                                  .stack().reset_index(level=1, drop=True)
                                                  .rename('Metabolite'))

    fileData['Metabolite'] = fileData['Metabolite'].str.lower()
    pubChemQuery = fileData['Metabolite'].tolist()

    all_compounds = []
    for metabolite in pubChemQuery:
        try:
            df = pcp.get_substances(identifier=metabolite,
                                    namespace='name',
                                    as_dataframe=True)
            df['Name'] = metabolite
            print(df)
            all_compounds.append(df)
        except (KeyError, TimeoutError, pcp.TimeoutError):
            continue
    all_compounds = pd.concat(all_compounds)
    return all_compounds
    #all_compounds.to_csv('testPubChemAPIResults.csv')

def test_explode_pubchem_query(pubchemQuery):
    """
    text_explode_pubchem_query takes the synoynms from the PubChem query and explodes the list as several new rows in the pandas dataframe
    :param   pubchemQuery: A Pandas dataframe containing the PubChem mapping
    :return: queryList:
    """
    # Explode rows from PubChem query
    pubchemQuery = pubchemQuery.drop('synonyms', axis=1).join(pubchemQuery['synonyms']
                                                              .str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('synonyms'))
    # Clean up regexes
    patterns = ["[", "]", "'", '"']
    for p in patterns:
        pubchemQuery['synonyms'] = pubchemQuery['synonyms'].str.replace(p, '')
    pubchemQuery['synonyms'] = pubchemQuery['synonyms'].str.lower()
    queryList = ';'.join(pubchemQuery['synonyms'].unique())
    print(queryList)
    return queryList

# Tests to perform

# 1. test_parse_pubchem_compounds
#filename = '~/Data/Expression/Metabolomics/ME1/raw/raw.xlsx'
#test_parse_pubchem_compounds(filename)

# 2. test_explode_pubchem_query
filename = 'testPubChemAPIResults.csv'
pubchemQuery = pd.read_csv(filename)
test_explode_pubchem_query(pubchemQuery)