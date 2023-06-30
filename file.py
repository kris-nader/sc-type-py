import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import requests
import xml.etree.ElementTree as ET
import scanpy as sc
from sklearn.preprocessing import MinMaxScaler
from collections import defaultdict


## main


adata = sc.read_h5ad(filename = '/homes/knader/mitro2/merge_seurat_stromal_filtered.h5ad')
scaled_matrix = pd.DataFrame(adata.X.T, columns=adata.obs_names, index=adata.var_names)
scRNAseqData=scaled_matrix
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])

unique_clusters = pd.Series(adata.obs['seurat_clusters']).unique()

# Initialize an empty list to store results
results = []
# Iterate over each unique cluster
for cl in unique_clusters:
    # Subset the AnnData object to the current cluster
    subset_adata = adata[adata.obs['seurat_clusters'] == cl]
    # Calculate row sums for `es.max` matrix using genes present in the current cluster
    es_max_cl = es.max.loc[:, subset_adata.obs_names].sum(axis=1)
    # Create a data frame with cluster, type, scores, and ncells columns
    df = pd.DataFrame({
        'cluster': cl,
        'type': es_max_cl.index,
        'scores': es_max_cl.values,
        'ncells': len(subset_adata)
    })
    # Append the data frame to the results list
    results.append(df)

# Concatenate the results data frames into a single data frame
cL_results = pd.concat(results)


# Group by cluster and select the top row based on the 'scores' column
sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = "Unknown"

# Print the first three columns of sctype_scores
print(sctype_scores.iloc[:, :3])

adata.obs['customclassif_k'] = ""

# Update the customclassif values based on sctype_scores
for cl in sctype_scores['cluster'].unique():
    cl_type = sctype_scores.loc[sctype_scores['cluster'] == cl, 'type']
    adata.obs.loc[adata.obs['seurat_clusters'] == cl, 'customclassif_k'] = cl_type.iloc[0]


## functions

def sctype_score(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
    marker_stat = defaultdict(int, {gene: sum(gene in genes for genes in gs.values()) for gene in set(gene for genes in gs.values() for gene in genes)})
    marker_sensitivity = pd.DataFrame({'gene_': list(marker_stat.keys()), 'score_marker_sensitivity': list(marker_stat.values())})
    # Rescaling the score_marker_sensitivity column
    scaler = MinMaxScaler(feature_range=(0, 1))
    marker_sensitivity['score_marker_sensitivity'] = scaler.fit_transform(marker_sensitivity['score_marker_sensitivity'].values.reshape(-1, 1))
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
    es = pd.DataFrame(
        index=gs__.keys(),
        columns=Z.columns,
        data=np.zeros((len(gs__), Z.shape[1]))
    )
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


def gene_sets_prepare(path_to_db_file, cell_type):
    cell_markers = pd.read_excel(path_to_db_file)
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    cell_markers['geneSymbolmore1'] = cell_markers['geneSymbolmore1'].str.replace(" ", "")
    cell_markers['geneSymbolmore2'] = cell_markers['geneSymbolmore2'].str.replace(" ", "")
    
    hgnc_table = pd.read_csv("hgnc_table.txt", names=["Symbol", "Approved_Symbol"], sep=" ", header=1)
    res = dict(zip(hgnc_table.Symbol, hgnc_table.Approved_Symbol))
    
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



