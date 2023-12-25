library(phyloseq)

rm(list=ls())
############################################phy#################################
TOFU_TSS_kraken2_output_NCBI_SRP151288 <- readRDS("~/data/Projects/github/Li_Project/pipeline/TOFU_16S/Kraken/input/TOFU_TSS_kraken2_output_NCBI_SRP151288.rds")
TOFU_TSS_kraken2_output_NCBI_SRP388727 <- readRDS("~/data/Projects/github/Li_Project/pipeline/TOFU_16S/Kraken/input/TOFU_TSS_kraken2_output_NCBI_SRP388727.rds")
phy<- merge_phyloseq(TOFU_TSS_kraken2_output_NCBI_SRP151288,TOFU_TSS_kraken2_output_NCBI_SRP388727)
############################################Meta#################################
meta1= read.csv("SRP151288.txt")
row.names(meta1)= meta1$Run
meta1$Group=meta1$Sample.Name
meta1$Batch="SRP151288"
meta1=meta1 %>% dplyr::select(Group,Batch)
meta2= read.csv("SRP388727.txt")
row.names(meta2)= meta2$Run
meta2$Group = meta2$isolation_source
meta2$Batch="SRP388727"
meta2=meta2 %>% dplyr::select(Group,Batch)

Meta= rbind(meta1,meta2)
Meta$GROUP=gsub("Control.*","HC",Meta$Group)
Meta$GROUP=gsub("healthy conrtrol","HC",Meta$GROUP)
Meta$GROUP=gsub("thyroid cancer patient","TC",Meta$GROUP)

Meta$GROUP=gsub("Control.*","HC",Meta$GROUP)
Meta$GROUP=gsub("TC.*","TC",Meta$GROUP)
Meta$GROUP=gsub("thyroid nodules patient","TN",Meta$GROUP)

sample_data(phy) = Meta

table(phy@sam_data$Batch,phy@sam_data$GROUP)
saveRDS(phy,"output/Batch2_phy_meta.rds")





