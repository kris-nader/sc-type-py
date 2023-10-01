## how to make the xlsx file using the R function

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

gs_list = gene_sets_prepare("/Users/naderkri/Downloads/ScTypeDB_full.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

r=read.xlsx("/Users/naderkri/Downloads/ScTypeDB_full.xlsx")

k=data.frame(tissueType=rep("Immune system",length(names(gs_list$gs_positive))),
             cellName=names(gs_list$gs_positive),
             geneSymbolmore1=unlist(lapply(names(gs_list$gs_positive),function(x) {paste(unlist(gs_list$gs_positive[[x]]),collapse = ",")})),
             geneSymbolmore2=unlist(lapply(names(gs_list$gs_negative),function(x) {paste(unlist(gs_list$gs_negative[[x]]),collapse = ",")})))
k$geneSymbolmore2[which(k$geneSymbolmore2=="")]=NA

write.xlsx(k,file="/Users/naderkri/Desktop/SCTYPE_immune_corrected.xlsx")

write.table(pbmc@assays$RNA@scale.data,file="/Users/naderkri/Desktop/pbmc_scale.txt")
