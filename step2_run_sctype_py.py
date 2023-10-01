adata=sc.read_text("/Users/naderkri/Desktop/pbmc_scale.txt",first_column_names=True)
scaled_matrix_ko = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

scRNAseqData=scaled_matrix_ko

gs_list=gene_sets_prepare_corrected(path_to_db_file="/Users/naderkri/Desktop/SCTYPE_immune_corrected.xlsx",cell_type="Immune system")
print(gs_list)
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])


df=pd.DataFrame(es.max)
df.to_csv('/Users/naderkri/Desktop/data_pbmc.txt', sep='\t', index=True)

