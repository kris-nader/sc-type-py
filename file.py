import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import requests
import xml.etree.ElementTree as ET
import scanpy as sc
from sklearn.preprocessing import MinMaxScaler
from collections import defaultdict
import concurrent.futures
import multiprocessing
from functools import partial



def gene_sets_prepare(path_to_db_file, cell_type):
    cell_markers = pd.read_excel(path_to_db_file)
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    #hgnc_table = pd.read_csv("/homes/knader/sctypePy/hgnc_table.txt", names=["Symbol", "Approved_Symbol"], sep=" ", header=1)
    #res = dict(zip(hgnc_table.Symbol, hgnc_table.Approved_Symbol))
    posmarkers = cell_markers['geneSymbolmore1'].dropna().str.split(',', expand=True).stack().reset_index(drop=True)
    negmarkers = cell_markers['geneSymbolmore2'].dropna().str.split(',', expand=True).stack().reset_index(drop=True)
    # Concatenate the two columns and remove 'None' values
    gene_names = pd.concat([posmarkers, negmarkers]).drop_duplicates().reset_index(drop=True)
    gene_names = gene_names[gene_names != 'None']
    gene_names=set(gene_names)
    res=return_gene_symbol(gene_names)
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].apply(lambda row: process_gene_symbols(row, res))
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].apply(lambda row: process_gene_symbols(row, res))
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace("///", ",")
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace("///", ",")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    return {'gs_positive': gs_positive, 'gs_negative': gs_negative}

def process_gene_symbols(gene_symbols, res):
    if pd.isnull(gene_symbols):
        return ""
    markers_all = gene_symbols.split(',')
    markers_all = [marker.strip().upper() for marker in markers_all if marker.strip().upper() not in ['NA', '']]
    markers_all = sorted(markers_all)
    if len(markers_all) > 0:
        markers_all = [res.get(marker) for marker in markers_all]
        markers_all = [symbol for symbol in markers_all if symbol is not None]
        markers_all = list(set(markers_all))
        return ','.join(markers_all)
    else:
        return ""

def check_gene_symbol(gene_symbol):
    url = f'https://rest.genenames.org/search/{gene_symbol}'
    response = requests.get(url)
    
    # Parse the XML data
    root = ET.fromstring(response.text)
    
    # Check if the response has results
    result_elem = root.find("result")
    if result_elem is not None and result_elem.get("numFound") != "0":
        # Get the maximum score value
        max_score = float(result_elem.get("maxScore"))
        
        # Find the doc element with the maximum score
        max_score_doc = None
        for doc in result_elem.findall("doc"):
            score = float(doc.find("float[@name='score']").text)
            if score == max_score:
                max_score_doc = doc
                break
        
        # Retrieve the symbol value
        symbol_elem = max_score_doc.find("str[@name='symbol']")
        if symbol_elem is not None:
            symbol = symbol_elem.text
        else:
            symbol = None
    else:
        symbol = None
    
    return symbol


def process_gene(gene):
    return gene,check_gene_symbol(gene)

def return_gene_symbol(gene_list):
	with multiprocessing.Pool() as pool:
		results = pool.map(partial(process_gene), gene_list)
	# Create a dictionary from the results
	result_dict = dict(results)
	return result_dict



def sctype_score(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
    marker_stat = defaultdict(int, {gene: sum(gene in genes for genes in gs.values()) for gene in set(gene for genes in gs.values() for gene in genes)})
    marker_sensitivity = pd.DataFrame({'gene_': list(marker_stat.keys()), 'score_marker_sensitivity': list(marker_stat.values())})
    # Rescaling the score_marker_sensitivity column
    min_value = marker_sensitivity['score_marker_sensitivity'].min()
    max_value = marker_sensitivity['score_marker_sensitivity'].max()
    # Apply the formula to the column
    marker_sensitivity['score_marker_sensitivity'] = 1 - (marker_sensitivity['score_marker_sensitivity'] - min_value) / (max_value - min_value)
    # Convert gene names to Uppercase
    if gene_names_to_uppercase:
        scRNAseqData.index = scRNAseqData.index.str.upper()
    # Subselect genes only found in data
    names_gs_cp = list(gs.keys())
    names_gs_2_cp = list(gs2.keys())
    gs_ = {key: [gene for gene in scRNAseqData.index if gene in gs[key]] for key in gs}
    gs2_ = {key: [gene for gene in scRNAseqData.index if gene in gs2[key]] for key in gs2}
    gs__ = dict(zip(names_gs_cp, gs_.values()))
    gs2__ = dict(zip(names_gs_2_cp, gs2_.values()))
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity['gene_'].isin(set.union(*map(set, gs__.values())))]
    # Z-scale if not
    if not scaled:
        Z = scale(scRNAseqData.T).T
    else:
        Z = scRNAseqData
    # Multiply by marker sensitivity
    for _, row in cell_markers_genes_score.iterrows():
        Z.loc[row['gene_']] *= row['score_marker_sensitivity']
    marker_genes=list(set().union(*gs__.values(), *gs2__.values()))
    Z = Z.loc[marker_genes]
    # Combine scores
    es = pd.DataFrame(index=gs__.keys(),columns=Z.columns,data=np.zeros((len(gs__), Z.shape[1])))
    print(es)
    for gss_, genes in gs__.items():
        for j in range(Z.shape[1]):
            gs_z = Z.loc[genes, Z.columns[j]]
            gz_2 = Z.loc[gs2__[gss_], Z.columns[j]] * -1 if gs2__ and gss_ in gs2__ else pd.Series(dtype=np.float64)
            sum_t1 = np.sum(gs_z) / np.sqrt(len(gs_z))
            sum_t2 = np.sum(gz_2) / np.sqrt(len(gz_2)) if not gz_2.empty else 0
            if pd.isna(sum_t2):
                sum_t2 = 0
            es.loc[gss_, Z.columns[j]] = sum_t1 + sum_t2
    es = es.dropna(how='all')
    return es

def gene_sets_prepare_corrected(path_to_db_file, cell_type):
    cell_markers = pd.read_excel(path_to_db_file)
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace("///", ",")
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace("///", ",")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    return {'gs_positive': gs_positive, 'gs_negative': gs_negative} 




adata=sc.read_text("Irf7ko_gene_expression_fixed.txt",first_column_names=True)
scaled_matrix_ko = pd.DataFrame(adata.X.T, columns=adata.obs_names, index=adata.var_names)

scRNAseqData=scaled_matrix_ko

#path_to_db_file="/homes/knader/sctypePy/Lung_fibroblast_markers.xlsx"
#cell_type="Lungtissue"

gs_list=gene_sets_prepare(path_to_db_file="/homes/knader/sctypePy/Lung_fibroblast_markers.xlsx",cell_type="Lungtissue")
print(gs_list)

gs_list=gene_sets_prepare_corrected(path_to_db_file="/homes/knader/sctypePy/Lung_fibroblast_markers_corrected.xlsx",cell_type="Lungtissue")


es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
pritn(es.max)
