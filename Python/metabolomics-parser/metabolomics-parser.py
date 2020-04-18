"""
metabolomics-parser.py uses the common metabolite names in a metabolomics dataset and maps them to standard
metabolite identifiers in established databases. This is used for mapping metabolites to constraint-based metabolic
models (sbml-based namespaces only).
@author: Scott Campit
"""
import sys, os
import re
import pandas as pd
import numpy as np

def queryPubChem(data, filename, sheet=''):
    """
    queryPubChem maps the metabolite name from a pandas Dataframe in the 'Compound Method' column and extracts
    synoynms from several databases using the PubChem API.
    :param data:           A Pandas Dataframe of the metabolomics dataframe with the common metabolite identifiers
                           under the 'Compound Method' column
    :return all_compounds: A Pandas Dataframe from the PubChem API containing the metabolite map. This dataframe is
                           saved as a .csv file.
    :return queryList:     A string with semicolon delimters to be fed into a REST-API
    """

    import pubchempy as pcp

    # Split 'Compound Method' column by the '/' regex and clean up some data
    data = data.drop('Compound Method', axis=1).join(data['Compound Method']
                                                 .str.split('/', expand=True)
                                                 .stack().reset_index(level=1, drop=True)
                                                 .rename('Metabolite'))
    data['Metabolite'] = data['Metabolite'].str.lower()
    pubChemQuery = data['Metabolite'].tolist()

    # Mine the PubChem database for synonyms
    all_compounds = []
    print("Mapping metabolite names to PubChem database for synonym matching and ID retrieval.")
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
    all_compounds.to_csv(filename + '_' + sheet +'_pubchem.csv')

    # Read in assembled database as pandas dataframe. It transforms the nested list in the Pandas Dataframe to a
    # string that can be manipulated using regexs
    all_compounds = pd.read_csv(filename + '_' + sheet +'_pubchem.csv')

    # Explode synonyms column
    all_compounds = all_compounds.drop('synonyms', axis=1).join(all_compounds['synonyms']
                                                              .str.split(',', expand=True)
                                                              .stack().reset_index(level=1, drop=True)
                                                              .rename('synonyms'))

    # Clean up regexes and get unique compounds
    patterns = ["[", "]", "'"]
    for p in patterns:
        all_compounds['synonyms'] = all_compounds['synonyms'].str.replace(p, '')
    all_compounds['synonyms'] = all_compounds['synonyms'].str.lower()

    # Make this JSON serializable
    queryList = ';'.join(all_compounds['synonyms'].unique())
    return all_compounds, queryList

def queryMetaboAnalyst(filename, sheet='Sheet1'):
    """
    queryMetaboAnalyst takes the file (currently only written for csv files),
    and uses the first column as metabolite names in the file. 
    
    The metabolite names are queried using MetaboAnalyst. The output is a metabolite map consisting of
    several different identifiers for the given metabolite.

    :param   filename: A string denoting the path to a text delimited or Excel file containing the metabolomics data
    :return: data:     A Pandas dataframe with a metabolite map containing several different metabolite identifiers
    """

    import requests
    import json

    print('Mapping metabolomics data to additional identifiers')
    if os.path.splitext(filename)[1] is '.xls' or 'xlsx':
        fileData = pd.read_excel(filename, sheet_name=sheet)
    else:
        fileData = pd.read_csv(filename)

    # Get synonyms of each metabolite from PubChem
    all_compounds, queryList = queryPubChem(fileData, filename, sheet)

    # Query metaboAnalyst for additional database identifiers
    queryDict = {"queryList" : queryList,
             "inputType" : "name"}
    metabolites = json.dumps(queryDict)

    # Query MetaboAnalyst using data set metabolites
    url = "http://api.xialab.ca/mapcompounds"
    headers = {
        'Content-Type': "application/json",
        'cache-control': "no-cache",
    }
    response = requests.request("POST", url, data=metabolites, headers=headers).json()
    data = pd.DataFrame(response)

    identifiers = ['hmdb_id', 'kegg_id', 'pubchem_id', 'chebi_id', 'chebi_id', 'metlin_id']
    for id in identifiers:
        data[id] = data[id].replace('-', np.nan)
    data['chebi_id'] = data['chebi_id'].astype(float)

    print('MetaboAnalyst query done!')
    print(data)

    return data

def mapMetabolicModel(model):
    """
    mapMetabolicModel takes in a metabolic model (xml format only), and parses the
    identifiers in the map. 
    
    Currently only tested with RECON1 - xml namespaces will change depending on metabolic models.

    :param  model:    A string describing the path to the metabolic model file (.xml or .sbml file)
    :return modelMap: A Pandas dataframe containing metabolite identifiers from the metabolic model
    """

    from lxml import etree
    print('Parsing metabolic model to get metabolite and associated identifiers')
    context = etree.parse(model)

    # Namespaces for different databases in the sbml models
    chebi_pattern = re.compile(r'.*http://identifiers.org/chebi/CHEBI:(\d+)')
    hmdb_pattern = re.compile(r'.*http://identifiers.org/hmdb/HMDB(\d+)')
    kegg_pattern = re.compile(r'.*http://identifiers.org/kegg.compound/C(\d+)')

    modelMap = pd.DataFrame()
    for metabolite in context.iter(tag='{http://www.sbml.org/sbml/level3/version1/core}species'):
        name = metabolite.get("name")
        bigg = metabolite.get("metaid")
        for element in metabolite.iter():
            if element.tag == '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li':
                content = etree.tostring(element)
                content_string = content.decode("utf-8")

                metabolite_name = []
                chebi = []
                hmdb = []
                kegg = []

                metabolite_name.append(name)
                metabolite_name = str(metabolite_name).replace(
                    '[', '').replace(']', '')

                if re.match(chebi_pattern, content_string):
                    chebi_num = re.findall(chebi_pattern, content_string)
                    chebi.append(chebi_num)
                else:
                    chebi.append(np.nan)
                chebi = str(chebi).replace('[', '').replace(']', '')

                if re.match(hmdb_pattern, content_string):
                    hmdb_num = re.findall(hmdb_pattern, content_string)
                    hmdb_num = ["HMDB" + h for h in hmdb_num]
                    hmdb.append(hmdb_num)
                else:
                    hmdb.append(np.nan)
                hmdb = str(hmdb).replace('[', '').replace(']', '')

                if re.match(kegg_pattern, content_string):
                    kegg_num = re.findall(kegg_pattern, content_string)
                    kegg_num = ["C" + k for k in kegg_num]
                    kegg.append(kegg_num)
                else:
                    kegg.append(np.nan)
                kegg = str(kegg).replace('[', '').replace(']', '')

                # Format stuff correctly before saving
                modelMap = modelMap.append({'Metabolite':metabolite_name,
                                            'BIGG':bigg, 'HMDB':hmdb,
                                            'CHEBI':chebi, 'KEGG':kegg},
                                            ignore_index=True)
        element.clear()

    modelMap['CHEBI'] = modelMap['CHEBI'].str.replace("'", "")
    modelMap['CHEBI'] = modelMap['CHEBI'].astype(float)
    print("Metabolic map complete")
    print(modelMap)
    return modelMap

def matchModelAndData(data, modelMap):
    """
    matchModelAndData uses the identifiers from the metabolomics data and the model
    to find matches between them.

    :param  data:               A Pandas dataframe containing data queried from MetaboAnalyst.
    :param  modelMap:           A Pandas dataframe queried from mapping the metabolite names to the COBRA metabolic
                                model.
    :return mergedModelDataMap: A Pandas dataframe of the merged map between the metabolomics data and the
                                metabolic model.
    """

    print('Match metabolomics identifiers and model identifiers by ChEBI and KEGG IDs')
    chebi = pd.merge(modelMap, data, left_on='CHEBI', right_on='chebi_id')
    chebi = chebi[np.isfinite(chebi['CHEBI'])]
    chebi = chebi[["Metabolite", "query", "BIGG", "CHEBI"]]

    kegg = pd.merge(modelMap, data, left_on='KEGG', right_on='kegg_id')
    kegg = kegg[pd.notnull(kegg['KEGG'])]
    kegg = kegg[['Metabolite', 'query', 'BIGG', 'KEGG']]

    mergedModelDataMap = pd.merge(chebi, kegg,
                                  how='outer', on=['Metabolite', 'query', 'BIGG'])
    mergedModelDataMap = mergedModelDataMap.drop_duplicates(keep='first')

    print("Found matching metabolites based on ChEBI and KEGG identities!")
    print(mergedModelDataMap)
    return mergedModelDataMap

def mapMetabolitePositionsInModel(mergedModelDataMap, model):
    """
    mapMetabolitePositionsInModel gets the metabolite positions using the cobrapy library.

    :param  mergedModelDataMap: A Pandas dataframe of the merged map between the metabolomics data and the metabolic
                                model.
    :param  model:              A string denoting the path to the metabolic model (.xml or .sbml file).
    :return PositionModel:      A Pandas dataframe containing the positions for each metabolite in the metabolic model.
    """

    import cobra
    print("Mapping metabolite positions in metabolic model")
    bigg_ids = mergedModelDataMap['BIGG'].replace('M_', '')
    bigg_ids = list(bigg_ids)
    biggRxn = pd.DataFrame()
    mdl = cobra.io.read_sbml_model(model)

    for index, met in enumerate(mdl.metabolites):
        if any(str(met) in m for m in bigg_ids):
            biggRxn = biggRxn.append(pd.DataFrame({"Position": index,
                                                   "Metabolite": [met]}
                                                  ), ignore_index=True)

    biggRxn['Metabolite'] = biggRxn['Metabolite'].astype(str)
    biggRxn["Compartment"] = biggRxn["Metabolite"].str.rsplit('_').str[-1]
    biggRxn["Metabolite"] = biggRxn["Metabolite"].str.split('_').str[0]

    biggRxn = biggRxn.pivot_table(index='Metabolite',
                            columns='Compartment',
                            values='Position',
                            aggfunc='mean', fill_value=0)

    mergedModelDataMap['BIGG'] = mergedModelDataMap['BIGG'].str.split('_', n=1).str[-1]
    mergedModelDataMap['BIGG'] = mergedModelDataMap['BIGG'].str.rsplit('_', n=1).str[0]
    mergedModelDataMap = mergedModelDataMap.drop(['CHEBI', 'KEGG', 'Metabolite'], axis=1)
    mergedModelDataMap = mergedModelDataMap.drop_duplicates(keep='first')

    PositionModel = pd.merge(biggRxn, mergedModelDataMap,
                              left_index=True, right_on='BIGG')
    PositionModel = PositionModel.drop(['BIGG'], axis=1)
    PositionModel = PositionModel.set_index(['query'])

    print('Mapped metabolite positions in metabolic model to metabolite name')
    print(PositionModel)
    return PositionModel

def constructFinalDataset(PositionModel, filename, sheet='Sheet1'):
    """
    constructFinalDataset merges the array of metabolite positions and the file name together.
    It automatically outputs an Excel file, ready for piping into DFA.
    
    :param  PositionModel: A Pandas dataframe containing the metabolite positions in the metabolic model.
    :param  filename:      A string denoting the path of the metabolomics data.
    :return df:            A Pandas dataframe containing the metabolomics data that intersects with the metabolic
                           model and the metabolite positions in the metabolomic model.
    """

    print("Merging metabolomics data to model map")
    if os.path.splitext(filename)[1] is '.xls' or 'xlsx':
        fileData = pd.read_excel(filename, sheet_name=sheet)
    else:
        fileData = pd.read_csv(filename)
    df = pd.merge(PositionModel, fileData,
                  left_index=True, right_on=fileData.iloc[:, 0],
                  how='inner')
    #df = df.drop(['Compound Method'], axis=1)
    print(df)
    df = df['key_0'].rename(columns={'key_0':'Metabolites'})
    name = filename.split('.')[0] + "_processed"
    df.to_excel(name+'.xlsx', sheet_name=sheet, index=False)
    print("Finished merging metabolomics data to model map. Results are saved as an Excel file!")

if __name__=='__main__':
   data = queryMetaboAnalyst(filename=sys.argv[2], sheet=sys.argv[3])
   modelMap = mapMetabolicModel(model=sys.argv[1])
   mergedModelDataMap = matchModelAndData(data, modelMap)
   PositionModel = mapMetabolitePositionsInModel(mergedModelDataMap, model=sys.argv[1])
   df = constructFinalDataset(PositionModel, filename=sys.argv[2], sheet=sys.argv[3])
