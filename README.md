
# ScType: Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data

**Article**: [https://doi.org/10.1038/s41467-022-28803-w]

<p style="text-align:justify;"> <b>ScType</b> a computational method for automated selection of marker genes based merely on scRNA-seq data. The open-source portal (<a href="//sctype.app">http://sctype.app</a>) provides an interactive web-implementation of the method.</p>

This GitHub covers the implementation of scType ,originally developed in R,  in python. 

##
<br><br>

![alt text](https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypePlan.png)

<br><br>

## Quick start
We will show the implementation using interactive python
```python
python3 -i ./sptypepy.py 
```
Once interactive python has been established, we can easily run functions that are equivalent to the sctype R implementation. Some differences include the query HGNC for approved gene symbols. In the R implementation, HGNChelper::checkGeneSymbols() is used. However, as no suitable equivalent has been found, rest.genenames API is used. 
```python
# Load the data- this is scaled 
adata = sc.read_text("/Users/naderkri/Desktop/sptype/pbmc_scaled.txt",first_column_names=True)
scRNAseqData = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
gs_list = gene_sets_prepare(path_to_db_file="/Users/naderkri/Downloads/ScTypeDB_full.xlsx",cell_type="Immune system")
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])


```
For the sake of showing that scType-py and scType-R result in the same annotations, we will export ex_max to a txt file and use R to overlay the scType-R annotations and the scType-py annotations. 

<br>
<p align="center">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/sctype_py_R.png" width="100%">
  <img src="https://github.com/kris-nader/sc-type-py/blob/main/proof.png" width="70%">
</p>
<br>














