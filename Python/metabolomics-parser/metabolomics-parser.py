"""
metabolomics-parser.py uses the common metabolite names in a metabolomics dataset and maps them to standard
metabolite identifiers in established databases. This is used for mapping metabolites to constraint-based metabolic
models (sbml-based namespaces only).

TODO:
  * Figure out robust way to query large datasets

@author: Scott Campit
"""
import sys, os
import re
import pandas as pd
import numpy as np


def queryPubChem(data):
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

    # The data is too large to keep in memory. So I wrote it into a csv file, and will read in.
    for metabolite in pubChemQuery:
        try:
            df = pcp.get_substances(identifier=metabolite,
                                    namespace='name',
                                    as_dataframe=True)
            df['Name'] = metabolite
            df = df.applymap(str)
            df = df.drop('synonyms', axis=1).join(df['synonyms']
                                            .str.split(',', expand=True)
                                            .stack().reset_index(level=1, drop=True)
                                            .rename('synonyms'))
            df.to_csv('~/Data/Mappings/ME1/pubmed_me1_query.csv',
                      mode='a',
                      header=False, index=False)

        except (KeyError, TimeoutError, pcp.TimeoutError):
            continue

    print("Finished metabolite common name -> identifier synoynm matching!")


def queryMetaboAnalyst(filename='', sheet='Sheet1', synmatch=True):
    """
    queryMetaboAnalyst takes the file (currently only written for csv files),
    and uses the first column as metabolite names in the file.

    The metabolite names are queried using MetaboAnalyst. The output is a metabolite map consisting of
    several different identifiers for the given metabolite.

    INPUTS:
    :param   filename: A string denoting the path to a text delimited or Excel file containing the metabolomics data
    :param   sheet:    A string denoting the tab name to read in
    :param   synmatch: A boolean flag determining whether to perform synonym matching or not

    OUTPUT:
    :return: data:     A Pandas dataframe with a metabolite map containing several different metabolite identifiers
    """

    import requests
    import json

    print('Mapping metabolomics data to additional identifiers')
    if os.path.splitext(filename)[1] is '.xls' or 'xlsx':
        fileData = pd.read_excel(filename, sheet_name=sheet)
    else:
        fileData = pd.read_csv(filename)

    # Load from Google Drive
    # wb = gc.open_by_url(filename)
    # wks = wb.worksheet(sheet)
    # data = wks.get_all_values()
    # fileData = pd.DataFrame(data)
    # header = fileData.iloc[0]
    # fileData = fileData[1:]
    # fileData.columns = header

    patterns = ["[", "]", "'", '"', '.']
    if synmatch is True:
        # Get synonyms of each metabolite from PubChem
        #queryPubChem(fileData)
        all_compounds = pd.read_csv(r'~/Data/Mappings/ME1/pubmed_me1_query.csv', chunksize=1E3)
        for chunk in all_compounds:
            chunk['Name'] = chunk['Name'].astype(str)

            # Clean up regexes and get unique compounds
            for p in patterns:
                chunk['Name'] = chunk['Name'].str.replace(p, '')
                chunk['Name'] = chunk['Name'].str.strip()
            chunk['Name'] = chunk['Name'].str.lower()
            chunk['Compound Method'] = chunk['Name']

            queryList = ';'.join(chunk['Compound Method'].unique())

            # Query metaboAnalyst for additional database identifiers
            queryDict = {"queryList": queryList,
                         "inputType": "name"}
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
                data[id] = data[id].astype(str)

            data.to_csv('~/Data/Mappings/ME1/metaboanalyst_me1_query.csv', mode='a', header=False)
            print('MetaboAnalyst query done!')
        data = list()

    else:
        all_compounds = fileData
        for p in patterns:
            all_compounds['Compound Method'] = all_compounds['Compound Method'].str.replace(p, '')
            all_compounds['Compound Method'] = all_compounds['Compound Method'].str.strip()
        all_compounds['Compound Method'] = all_compounds['Compound Method'].str.lower()

        queryList = ';'.join(all_compounds['Compound Method'].unique())

        # Query metaboAnalyst for additional database identifiers
        queryDict = {"queryList": queryList,
                     "inputType": "name"}
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
            data[id] = data[id].astype(str)
        print('MetaboAnalyst query done!')

    return data


def mapMetabolicModel(model):
    """
    mapMetabolicModel takes in a metabolic model (xml format only), and parses the
    identifiers in the map. 

    Currently only tested with RECON1 - xml namespaces will change depending on metabolic models.

    :param  model:    A string describing the path to the metabolic model file (`.xml` or `.sbml` file types supported only.
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

                # Get metabolite identifiers for CHEBI, HMDB and KEGG
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
                modelMap = modelMap.append({'Metabolite': metabolite_name,
                                            'BIGG': bigg, 'HMDB': hmdb,
                                            'CHEBI': chebi, 'KEGG': kegg},
                                           ignore_index=True)
        element.clear()

    for col in modelMap:
        modelMap[col] = modelMap[col].astype(str)
        modelMap[col] = modelMap[col].str.replace("'", "")
    #modelMap = modelMap.drop_duplicates(keep='first')
    modelMap.to_csv('~/Data/Mappings/ME1/RECON1_ID_Map.csv', index=False)
    print("Metabolic map complete")
    return modelMap


def matchModelAndData(data, modelMap, synmatch=True):
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

    if synmatch is True:
        data = pd.read_csv(r'/home/scampit/Data/Mappings/ME1/metaboanalyst_me1_query.csv', dtype=str, chunksize=1E1)
        for chunk in data:
            chebi = pd.merge(modelMap, chunk, left_on='CHEBI', right_on='chebi_id', how='inner')
            chebi = chebi[["Metabolite", "query", "BIGG", "CHEBI"]]
            chebi = chebi.dropna()

            kegg = pd.merge(modelMap, chunk, left_on='KEGG', right_on='kegg_id', how='inner')
            kegg = kegg[['Metabolite', 'query', 'BIGG', 'KEGG']]
            kegg = kegg.dropna()
            merged_data = pd.merge(chebi, kegg,
                                          how='inner', on=['Metabolite', 'query', 'BIGG'])
            merged_data.to_csv(r'~/Data/Mappings/ME1/metaboanalyst_recon1_map.csv', mode='a', header=False, index=False)
    else:
        chebi = pd.merge(modelMap, data, left_on='CHEBI', right_on='chebi_id')
        chebi = chebi[["Metabolite", "query", "BIGG", "CHEBI"]]

        kegg = pd.merge(modelMap, data, left_on='KEGG', right_on='kegg_id')
        kegg = kegg[['Metabolite', 'query', 'BIGG', 'KEGG']]
        merged_data = pd.merge(chebi, kegg,
                                      how='inner', on=['Metabolite', 'query', 'BIGG'])
        merged_data = merged_data.drop_duplicates('query', keep='first')
    print("Found matching metabolites based on ChEBI and KEGG identities!")
    return merged_data

def mapMetabolitePositionsInModel(mergedModelDataMap, model):
    """
    mapMetabolitePositionsInModel gets the metabolite positions using the cobrapy library.

    :param  mergedModelDataMap: A Pandas dataframe of the merged map between the metabolomics data and the metabolic
                                model.
    :param  model:              A string denoting the path to the metabolic model (`.xml` or `.sbml` file types supported only.
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
    #mergedModelDataMap = mergedModelDataMap.drop(['CHEBI', 'KEGG'], axis=1)
    mergedModelDataMap = mergedModelDataMap.drop_duplicates('Metabolite', keep='first')

    PositionModel = pd.merge(biggRxn, mergedModelDataMap,
                              left_index=True, right_on='BIGG')
    PositionModel = PositionModel.drop(['BIGG'], axis=1)
    PositionModel = PositionModel.set_index(['Query'])
    PositionModel.to_csv(r'~/Data/Mappings/ME1/RECON1_position_map.csv', index=True)

    print('Mapped metabolite positions in metabolic model to metabolite name')
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

    # Load from Google Drive
    # wb = gc.open_by_url(filename)
    # wks = wb.worksheet(sheet)
    # data = wks.get_all_values()
    # fileData = pd.DataFrame(data)
    # header = fileData.iloc[0]
    # fileData = fileData[1:]
    # fileData.columns = header

    # Clean up regexes
    all_compounds = fileData
    patterns = ["[", "]", "'", '"', '.']
    for p in patterns:
        all_compounds['Compound Method'] = all_compounds['Compound Method'].str.replace(p, '')
        all_compounds['Compound Method'] = all_compounds['Compound Method'].str.strip()
        PositionModel.index = PositionModel.index.str.replace(p, '')
        PositionModel.index = PositionModel.index.str.strip()

    all_compounds['Compound Method'] = all_compounds['Compound Method'].str.lower()
    PositionModel.index = PositionModel.index.str.lower()
    print(PositionModel.index)
    print(all_compounds['Compound Method'].values)

    df = pd.merge(PositionModel, all_compounds,
                  left_index=True, right_on='Compound Method',
                  how='inner')
    #print(df)
    print("Finished merging metabolomics data to model map!")
    return df

if __name__=='__main__':
    name = r'/home/scampit/Data/Expression/Metabolomics/ME1/raw/ME1_Metabolomics.xlsx'
    model = r'/home/scampit/Data/CBM/MetabolicModels/RECON1/RECON1.xml'

    # Create a spreadsheet that will save all of the data
    from openpyxl import load_workbook
    writer = pd.ExcelWriter(r'/home/scampit/Data/Expression/Metabolomics/ME1/processed/ME1_mapped_metabolomics.xlsx',
                            engine='openpyxl')

    xl = pd.ExcelFile(name)
    sheetNames = xl.sheet_names
    sheetNames = sheetNames[1:]

    # Parts of the parser that will take a long time, especially if `synmatch` is True
    #modelMap = mapMetabolicModel(model)
    modelMap = pd.read_csv(r'~/Data/Mappings/ME1/RECON1_ID_Map.csv', dtype=str)
    #data = queryMetaboAnalyst(filename=name, sheet=sheetNames[0], synmatch=True)
    data = []
    #data = pd.read_csv(r'/home/scampit/Data/Mappings/ME1/metaboanalyst_me1_query.csv')
    #mergedModelDataMap = matchModelAndData(data, modelMap, synmatch=True)
    mergedModelDataMap = pd.read_csv(r'~/Data/Mappings/ME1/metaboanalyst_recon1_map.csv')
    PositionModel = mapMetabolitePositionsInModel(mergedModelDataMap, model)
    #PositionModel = pd.read_csv(r'~/Data/Mappings/ME1/RECON1_position_map.csv', index_col='Query')
    #df = constructFinalDataset(PositionModel, name, sheetNames[0])
    #print(df)

    # Save multiple sheets
    #for sht in sheetNames:
    #    df = constructFinalDataset(PositionModel, name, sht)
    #    df.to_excel(writer, sheet_name=sht, index=False)
    #writer.save()
