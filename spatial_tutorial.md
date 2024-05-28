

# SpType: ScType enables fast and accurate cell type identification from spatial transcriptomics data


**Article**: TBA

In this study, we adapt and showcase the application of scType, renowned for its speed, transparency, and user-friendly interface, to efficiently annotate cell types in **spatial transcriptomics data**.

  
Please refer to the original <a href="https://www.nature.com/articles/s41467-022-28803-w" target="_blank">ScType paper</a>  and <a href="https://github.com/IanevskiAleksandr/sc-type" target="_blank">Github</a> for more information on the method.


<br>

![alt text](https://github.com/kris-nader/sp-type/blob/main/sctype_goes_spatial_fig.png)


<br>
We start with by loading all required functions for this analysis. These are the original scType functions rewritten in python and applied on spatial transcriptomics data. 

```python

import urllib.request
import scanpy as sc
import numpy as np
import pandas as pd
np.random.seed(100)


# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)
```
Load sample data available through Scanpy. Specifically, we will use the "Mouse Brain Sagittal Anterior" dataset generated using V1 Chromium technology. 
In this tutorial, we will follow the preprocessing steps outlined by Scanpy in their basic analysis tutorial see <a href="https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html" target="_blank">scanpy tutorial for more details</a>.
Please note that variations in preprocessing choices can lead to different annotations(SCtransform vs LogNormalize or Leiden vs Louvain). Additionally, the annotations produced can vary significantly between different packages, such as Seurat and Scanpy, as discussed in this <a href="https://www.biorxiv.org/content/10.1101/2024.04.04.588111v2.abstract" target="_blank">paper</a>.

```python
adata = sc.datasets.visium_sge(sample_id="V1_Mouse_Brain_Sagittal_Anterior")
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",layer="counts")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata

# Scale and run PCA
sc.pp.scale(adata,max_value=10)
scaled_data = pd.DataFrame(adata.X)
# change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T


sc.tl.pca(adata,zero_center=False,random_state=0)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10,use_rep="X_pca",random_state=0)
sc.tl.leiden(adata,resolution=0.8,n_iterations=10)

# Visualize clusters using UMAP
sc.tl.umap(adata,min_dist=0.3)
sc.pl.umap(adata, color=['leiden'])
sc.pl.spatial(adata, img_key="hires", color=["leiden"])
```


Now, let's run the sctype functions. This involves several steps: querying the HGNC database for approved gene symbols, calculating sctype scores, aggregating these scores based on cluster information, and overlaying the results onto a lower-dimensional space such as UMAP or t-SNE.
Users are welcome to use their own custom marker datasets. For this example, we will use the default scTypeDB, which contains annotations for various healthy tissues. For more detailed information, please refer to the <a href="https://www.nature.com/articles/s41467-022-28803-w" target="_blank">original paper</a>.

```python

scRNAseqData=scaled_data
gs_list=gene_sets_prepare(path_to_db_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx" ,cell_type="Brain")

es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
unique_clusters = adata.obs['leiden'].unique()
# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,adata,es_max,"leiden") for cluster in unique_clusters])

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
sc.pl.spatial(adata, img_key="hires", color=["sctype_classification"])
```

<br>
<p align="center">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/sptype_py.png" width="70%">
</p>
<br>

