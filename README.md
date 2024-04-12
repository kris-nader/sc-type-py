
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














