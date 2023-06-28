import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import requests
import xml.etree.ElementTree as ET


def sctype_score(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
    # Check input matrix
    if not isinstance(scRNAseqData, np.ndarray) or not scRNAseqData.ndim == 2:
        print("scRNAseqData doesn't seem to be a matrix")
    elif np.sum(scRNAseqData.shape) == 0:
        print("The dimension of input scRNAseqData matrix equals to 0. Is it an empty matrix?")
    
    # Marker sensitivity
    marker_stat = pd.Series(np.unique(np.concatenate(gs)), name='gene').value_counts().sort_values(ascending=False)
    marker_sensitivity = pd.DataFrame({
        'score_marker_sensitivity': np.interp(marker_stat, (marker_stat.max(), 1), (0, 1)),
        'gene_': marker_stat.index.str.upper() if gene_names_to_uppercase else marker_stat.index
    })
    
    # Convert gene names to Uppercase
    if gene_names_to_uppercase:
        scRNAseqData.index = scRNAseqData.index.str.upper()
    
    # Subselect genes only found in data
    names_gs_cp = [name for name in gs]
    names_gs_2_cp = [name for name in gs2] if gs2 else None
    
    gs = [
        scRNAseqData.index[scRNAseqData.index.isin(gs_)].tolist()
        for gs_ in gs
    ]
    gs2 = [
        scRNAseqData.index[scRNAseqData.index.isin(gs2_)].tolist()
        for gs2_ in gs2
    ] if gs2 else None
    
    gs = dict(zip(names_gs_cp, gs))
    gs2 = dict(zip(names_gs_2_cp, gs2)) if gs2 else None
    
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity['gene_'].isin(np.unique(np.concatenate(gs)))]
    
    # Z-scale if not
    if not scaled:
        Z = scale(scRNAseqData.T).T
    else:
        Z = scRNAseqData
    
    # Multiply by marker sensitivity
    for _, row in cell_markers_genes_score.iterrows():
        Z.loc[row['gene_']] *= row['score_marker_sensitivity']
    
    # Subselect only with marker genes
    marker_genes = np.unique(np.concatenate([np.unique(np.concatenate(gs)), np.unique(np.concatenate(gs2))]))
    Z = Z.loc[marker_genes]
    
    # Combine scores
    es = pd.DataFrame(
        index=gs.keys(),
        columns=Z.columns,
        data=np.zeros((len(gs), Z.shape[1]))
    )
    for gss_, genes in gs.items():
        for j in range(Z.shape[1]):
            gs_z = Z.loc[genes, Z.columns[j]]
            gz_2 = Z.loc[gs2[gss_], Z.columns[j]] * -1 if gs2 and gss_ in gs2 else pd.Series(dtype=np.float64)
            sum_t1 = np.sum(gs_z) / np.sqrt(len(gs_z))
            sum_t2 = np.sum(gz_2) / np.sqrt(len(gz_2)) if not gz_2.empty else 0
            if pd.isna(sum_t2):
                sum_t2 = 0
            es.loc[gss_, Z.columns[j]] = sum_t1 + sum_t2
    
    es = es.dropna(how='all')
    
    return es




def gene_sets_prepare(path_to_db_file, cell_type):
    cell_markers = pd.read_excel(path_to_db_file)
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    #cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers.loc[:, 'geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    #cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    cell_markers.loc[:, 'geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    

    # read hgnc table from hgnc helper R package
    hgnc_table=pd.read_csv("hgnc_table.txt",names=["Symbol","Approved.Symbol"],sep=" ",header=1)
    res=dict(zip(hgnc_table.Symbol,hgnc_table.Approved_Symbol))
    
    # Correct gene symbols from the given DB (up-genes)
    cell_markers['geneSymbolmore1'] = cell_markers.apply(lambda row: process_gene_symbols(row['geneSymbolmore1']), axis=1)
    
    # Correct gene symbols from the given DB (down-genes)
    cell_markers['geneSymbolmore2'] = cell_markers.apply(lambda row: process_gene_symbols(row['geneSymbolmore2']), axis=1)
    
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace("///", ",")
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace("///", ",")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    
    gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(lambda x: x.str.split(',').explode().unique().tolist()).to_dict()
    
    return {'gs_positive': gs_positive, 'gs_negative': gs_negative}

def process_gene_symbols(gene_symbols):
    markers_all = gene_symbols.split(',')
    markers_all = [marker.strip().upper() for marker in markers_all if marker.strip().upper() not in ['NA', '']]
    markers_all = sorted(markers_all)
    
    if len(markers_all) > 0:
        markers_all = [check_gene_symbol(marker) for marker in markers_all]
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
    # Find the "result" element
    result = root.find("result")
    # Get the maximum score value
    max_score = float(result.get("maxScore"))
    # Find the doc element with the maximum score
    max_score_doc = None
    for doc in result.findall("doc"):
        score = float(doc.find("float[@name='score']").text)
        if score == max_score:
            max_score_doc = doc
            break
    # Retrieve the symbol value
    symbol_element = max_score_doc.find("str[@name='symbol']")
    if symbol_element is not None:
        symbol = symbol_element.text
    else:
        symbol = None
    return symbol


import time
start_time = time.time()
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
print("--- %s seconds ---" % (time.time() - start_time))

gene_symbols_list = []
for index, row in cell_markers.iterrows():
    gene_symbols = row['geneSymbolmore1'].split(',')
    gene_symbols_list.extend(gene_symbols)

gene_symbols_set=set(gene_symbols_list)

hgnc_table=pd.read_csv("hgnc_table.txt",names=["Symbol","Approved.Symbol"],sep=" ",header=1)
res=dict(zip(hgnc_table.Symbol,hgnc_table.Approved_Symbol))

p = {key: res[key] for key in gene_symbol_set if key in res}


