
# ScType: Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data

**Article**: [https://doi.org/10.1038/s41467-022-28803-w]

<p style="text-align:justify;"> <b>ScType</b> a computational method for automated selection of marker genes based merely on scRNA-seq data. The open-source portal (<a href="//sctype.app">http://sctype.app</a>) provides an interactive web-implementation of the method.</p>

This GitHub covers the implementation of scType ,originally developed in R,  in python. 
<br><br>
![alt text](https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypePlan.png)
<br><br>

## Quick start

```python
import urllib.request

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/file.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)
```

At this point, running functions equivalent to the sctype R implementation becomes straightforward. There are some differences, such as querying HGNC for approved gene symbols. In the R implementation, the <code>HGNChelper::checkGeneSymbols()</code> function is utilized. However, as no suitable equivalent has been identified, the rest.genenames API is employed instead.

```python
# Load the data- this is scaled 
adata = sc.read_text("/Users/naderkri/Desktop/sptype/pbmc_scaled.txt",first_column_names=True)
scRNAseqData = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
gs_list = gene_sets_prepare(path_to_db_file="/Users/naderkri/Downloads/ScTypeDB_full.xlsx",cell_type="Immune system")
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])


```
To validate the consistency between scType-py and scType-R annotations, we'll export ex_max to a txt file. Subsequently, we'll use R to overlay the annotations derived from scType-R alongside those from scType-py.

<br>
<p align="center">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/sctype_py_R.png" width="100%">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/proof.png" width="70%">
</p>
<br>

## Integrating scType in a scanpy workflow
On this page, we will highlight the use of sctype on scRNAseq data. For a tutorial on using sctype in python using spatial transcriptomics data, we refer users to the  <a href="https://github.com/kris-nader/sc-type-py/blob/main/spatial_tutorial.md" target="_blank">following link</a>. 
###  scType for scRNAseq data
```python
pd.set_option("display.precision", 9)
import urllib.request

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)
```

Process the single cell data using scanpy. There are slight differences in the results of the scanpy and Seurat workflows. For this reason, we have followed XXX paper discussing potential alterations to the default functions to achieve similar results to the Seurat workflow.

```python
import scanpy as sc
adata=sc.read_10x_mtx("/Users/naderkri/Downloads/filtered_gene_bc_matrices/hg19/")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize the data
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4,exclude_highly_expressed=False)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")

# Scale and run PCA
sc.pp.scale(adata,max_value=10)
scaled_data = pd.DataFrame(adata.X)
# change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T

decimal_places = 8
modified_matrix = {key: [round(val, decimal_places) for val in values] for key, values in scaled_data.items()}

mod=pd.DataFrame(modified_matrix)
mod.index = adata.var_names
scaled_data=mod

sc.tl.pca(adata, svd_solver="arpack",zero_center=False)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10,use_rep="X_pca")
sc.tl.leiden(adata,resolution=0.8,random_state=0,flavor="igraph",n_iterations=10,directed=False)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'])
```

Run scType with the following functions. These functions may be a bit slower than the original implementation in R. 



```python
scRNAseqData=scaled_data
gs_list=gene_sets_prepare(path_to_db_file="/Users/naderkri/Downloads/ScTypeDB_full.xlsx",cell_type="Immune system")
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])

unique_clusters = adata.obs['leiden'].unique()
# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,adata,es_max) for cluster in unique_clusters])

# Group by cluster and select the top row based on scores
sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

# Set low-confidence clusters to "Unknown"
sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'
adata.obs['sctype_classification'] = ""

# Iterate over unique clusters
for cluster in sctype_scores['cluster'].unique():
    # Filter sctype_scores for the current cluster
    cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
    # Get the type for the current cluster
    cl_type_value = cl_type['type'].iloc[0]
    # Update 'sctype_classification' in pbmc.obs for cells belonging to the current cluster
    adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_classification'] = cl_type_value

# Plot the UMAP with sctype_classification as labels
sc.pl.umap(adata, color='sctype_classification', title='UMAP with sctype_classification')
```






