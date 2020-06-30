"""
identifier_mapping gets several different ids using the mygene API.
@author: Scott Campit
"""

import mygene
import cobra
import pathlib
import pandas as pd

def getMapfromCOBRAGenes(modelFilePath, inputType='BiGG', outputType='All', Species='human'):
    """

    :param model:
    :param idType:
    :return:
    """
    # Read in COBRA model
    if pathlib.Path(modelFilePath).suffix is '.mat':
        model = cobra.io.load_matlab_model(modelFilePath)
    elif pathlib.Path(modelFilePath).suffix is '.sbml' or '.xml':
        model = cobra.io.read_sbml_model(modelFilePath)
    else:
        model = cobra.io.load_json_model(modelFilePath)

    # Get gene identifiers from COBRA model
    geneNames = []
    for gene in model.genes:
        geneNames.append(gene.id)

    # BiGG uses a combination of Entrez and Ensemble to denote isozymes. Get the real Entrez ID
    if inputType is 'entrezgene':
        geneNames = [x.split('.')[0] for x in geneNames]
    elif inputType is 'BiGG':
        geneNames = [x.split('_')[0] for x in geneNames]
        inputType = 'entrezgene'

    # If not specified, return everything
    if outputType is 'All':
        outputType=['accession.genomic', 'ec',
                    'ensembl.gene', 'entrezgene',
                    'ensembl.transcript',
                    'kegg', 'reactome',
                    'pdb', 'refseq',
                    'reporter', 'symbol',
                    'uniprot']

    # Use MyGene API to query identifiers as Pandas dataframe
    mg = mygene.MyGeneInfo()
    return mg.querymany(geneNames, scopes=inputType,
                        fields=outputType, species=Species,
                        as_dataframe=True)

def getMapFromList(idList, inputType='symbol', outputType='All', Species='human'):
    """

    :param model:
    :param idType:
    :return:
    """

    # If not specified, return everything
    if outputType is 'All':
        outputType=['entrezgene', 'kegg',
                    'ec', 'refseq',
                    'refseq.protein',
                    'HGNC', 'ensembl.gene',
                    'ensembl.protein',
                    'uniprot', 'pdb',
                    'humancyc', 'MIM',
                    'reactome']

    # Use MyGene API to query identifiers as Pandas dataframe
    mg = mygene.MyGeneInfo()
    return mg.querymany(idList, scopes=inputType,
                        fields=outputType, species=Species,
                        as_dataframe=True, df_index=True,
                        returnall=False)

if __name__ == "__main__":

    ## Get identifiers from metabolic model
    # Read in RECON1 from BiGG database
    #if sys.argv[1] is None:
    #    modelFilePath="/home/scampit/Data/CBM/MetabolicModels/RECON1/RECON1.xml"
    #else:
    #    modelFilePath=sys.argv[1]
    #df = getMapfromGenes(modelFilePath, inputType='BiGG', outputType=['symbol', 'ensembl.gene'])

    # Save everything to default genemap.csv in current working directory
    #if sys.argv[2] is None:
    #    df.to_csv('genemap.csv')
    #else:
    #    df.to_csv(sys.argv[2])

    ## Get identifiers from list of gene symbols
    fileName="~/Data/Mappings/MetabolicModelMaps/metabolic_map.xlsx"
    df = pd.read_excel(fileName, sheet_name='Genes')
    genelist = list(df['Gene symbol'])
    query = getMapFromList(genelist)
    print(query)
    query.to_csv('histoneMap.csv')
