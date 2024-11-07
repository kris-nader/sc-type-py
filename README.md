
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
import scanpy as sc
import numpy as np
import pandas as pd

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
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

### Load and cluster the data
We load a PBMC 3k example dataset. The raw data can be found <a href='https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz' download>here</a>.
```python
pd.set_option("display.precision", 9)
import urllib.request
import scanpy as sc
import numpy as np
import pandas as pd

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)
```

Process the single cell data using scanpy. There are slight differences in the results of the scanpy and Seurat workflows. For this reason, we have followed XXX paper discussing potential alterations to the default functions to achieve similar results to the Seurat workflow.

```python
import numpy as np
import pandas as pd

np.random.seed(100)
adata=sc.read_10x_mtx("./filtered_gene_bc_matrices/hg19/")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

Normalize, scale and cluster the data.
```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",layer="counts")
adata.raw = adata

# Scale and run PCA
sc.pp.scale(adata,max_value=10)
scaled_data = pd.DataFrame(adata.X)
# change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T


sc.tl.pca(adata,zero_center=False)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10,use_rep="X_pca")
sc.tl.leiden(adata, resolution=0.8)
sc.tl.umap(adata,min_dist=0.3)
sc.pl.umap(adata, color=['leiden'])
```

Now, we can  automatically assign cell types using ScType. For that, we first load 2 additional ScType functions: gene_sets_prepare and sctype_score
These functions may be a bit slower than the original implementation in R. As stated earlier, this is due to the fact that there is no well estabilshed package to accuratly retrieve approved gene symbols. For this reason, we decided to query the HGNC database for approved symbols. 

Users can prepare their gene input cell marker file or use the sctypeDB. The input XLSX must be formatted in the same way as the original scTypeDB. DB file should contain four columns (tissueType - tissue type, cellName - cell type, geneSymbolmore1 - positive marker genes, geneSymbolmore2 - marker genes not expected to be expressed by a cell type)


```python
scRNAseqData=scaled_data
gs_list=gene_sets_prepare(path_to_db_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",cell_type="Immune system")
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])

unique_clusters = adata.obs['leiden'].unique()
# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,adata,es_max,'leiden') for cluster in unique_clusters])

# Group by cluster and select the top row based on scores
sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

# Set low-confidence clusters to "Unknown"
sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

# Iterate over unique clusters
adata.obs['sctype_classification'] = ""
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

<p align="center">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/sctype_py.png" width="70%">
</p>




