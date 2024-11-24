

library(dplyr)
# 读取化合物信息
DF.TCM=readxl::read_xlsx("../../../../database/TCM/v20240929/TCMSP.xlsx",sheet = 1)
colnames(DF.TCM)=make.names(colnames(DF.TCM))

# 读取CID信息
DF.CID=read.csv("../../../../database/Metab.build.db/output/results/06_20w_name_info_nona_property_unique_separate.csv")
DF.TCM$compound.lower=tolower(DF.TCM$Compound.EN)


# 读取靶基因信息
DF.gene=read.csv("../../../../database/Metab.Gene/output/02_tidy_all.DB_tidy.csv") %>% dplyr::filter(Gene_Value>80)
DF.gene_collapse <- DF.gene %>% group_by(CID) %>% summarise(entrez_gene_symbol = paste(unique(entrez_gene_symbol), collapse = ";"))


# merge CID and gene
DF.CID.Gene=merge(DF.CID,DF.gene_collapse,by="CID",all.x=TRUE)

# merge DF.CID.Gene and TCM
DF.CID.compound.gene=merge(DF.TCM,DF.CID.Gene,by.x="compound.lower",by.y="Synonym.separate.lower",all.x=TRUE)




# 提取中药

DF.TCM_zhebeimu=DF.CID.compound.gene %>% filter(herb_cn_name =="浙贝母")
DF.TCM_xiakucao=DF.CID.compound.gene %>% filter(herb_cn_name =="夏枯草")


write.csv(DF.TCM_zhebeimu,"output/01_TCM_CID_compound_Gene_zhebeimu.csv",row.names = FALSE,quote = T)
write.csv(DF.TCM_xiakucao,"output/01_TCM_CID_compound_Gene_xiakucao.csv",row.names = FALSE,quote = T)




