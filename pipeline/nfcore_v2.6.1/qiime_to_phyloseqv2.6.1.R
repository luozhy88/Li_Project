print(" use qiime to phyloseq nf-core ampliseq v2.1.0
      /usr/bin/Rscript /home/zhalia//Desktop/projects/phyloseq_all/scripts/qiime_to_phyloseq.R -f /data/zhalia/nf-core/greenvalley/human/20210820_GV-HC
")

## base on ampliseq v2.1.0 nfcore
library(phyloseq)
#library(mixOmics)
library(viridis)
library(ggsci)
#library(microbiome)
#library(metagMisc)
#library(microbiome)
library(mixOmics)
library(metagMisc)
library(tidyverse)
library(dplyr)
library(DT)
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

# library(metacoder)
library(purrr)
#library(ggplot2)
library(microbiome)
library(ggplot2)
library(ggpubr)
library('grid')
# library(ggsci)
library(qiime2R)
######################   make sure you have a good sample_data for your phyloseq
######################
######################  change the group variable in the following script to match yours
######################   probably you also need to check the group levels?
library(optparse)

print("package finish")
option_list <- list(
  make_option(c("-f", "--first"), type = "character", default = FALSE,
              action = "store", help = "This is first!  phyloseq_work_dir"
  ),
  make_option(c("-s", "--second"), type = "character", default = FALSE,
              action = "store", help = "This is second! read the phyloseq_path"
  ),
  make_option(c("-t", "--third"), type = "character", default = FALSE,
              action = "store", help = "This is third! meta_filename_path"
  )
  # make_option(c("-h", "--help"), type = "logical", default = FALSE,
  #             action = "store_TRUE", help = "This is Help!"
  # )
)


opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))
print(opt$f)
#
#
#
workpath<-opt$f
print(workpath)
#workpath<-"/data2/zhiyu/software/nfcore/v2.1.0/xionghuang.No2"
setwd(workpath)
# browser()
# trace(qza_to_phyloseq,edit=T)
###############

# features <- read_qza(paste0(workpath,"/abundance_table/filtered/table.qza"))$data

# features <- read_q2biom(paste0(workpath, "/qiime2/abundance_tables/feature-table.biom"))
features <-   read.table(paste0(workpath, "/dada2/ASV_table.tsv"), sep = "\t", header = T,row.names = 1)
# phyloseq_rename<-function(x){
#   vec<-c()
#   for (i in 1:length(x)){
#     x[i]<-ifelse(grepl("^[[:digit:]]",x[i]),paste0('X',x[i]),x[i])
#     vec<-append(vec,x[i])
#     # print(x[i])
#   }
#   return(vec)
# }
# colnames(features)<- phyloseq_rename(colnames(features))

# taxonomy <- read_qza(paste0(workpath,"/taxonomy/taxonomy.qza"))$data
#taxonomy <- read_qza(paste0(workpath,"/qiime2/input/taxonomy.qza"))$data
# taxonomy <- read.table(paste0(workpath, "/qiime2/taxonomy/taxonomy.tsv"), sep = "\t", header = T)
taxonomy <- read.table(paste0(workpath, "/dada2/ASV_tax_species.silva_138.tsv"), sep = "\t", header = T,row.names = 1)
# taxonomy <- parse_taxonomy(taxonomy)
taxonomy <- as.matrix(taxonomy)[,1:7]


Taxonomy<-tax_table(taxonomy)
Otu_table<-otu_table(features, taxa_are_rows=T)
# rownames(Taxonomy) %in% rownames(features)



# trace(phyloseq,edit=T)
physeq<-phyloseq::phyloseq( Otu_table  , Taxonomy )

###############

#physeq<-qiime2R::qza_to_phyloseq(features=paste0(workpath,"/abundance_table/filtered/table.qza"), taxonomy=paste0(workpath,"/taxonomy/taxonomy.qza")) #, tree=paste0(workpath,"/phylogenetic_tree/rooted-tree.qza")

otu.table<-otu_table(physeq)

###############
library(phyloseq)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("remotes")
#remotes::install_github("jfq3/RDPutils")
#devtools::install_github("jfq3/RDPutils")
library(RDPutils)
# rep.seqs <- Biostrings::readDNAStringSet(paste0(workpath,"/representative_sequences/filtered/sequences.fasta"), format = "fasta")
rep.seqs <- Biostrings::readDNAStringSet(paste0(workpath,"/dada2/ASV_seqs.fasta"), format = "fasta")

phyloseq_new <- phyloseq::phyloseq(otu.table, rep.seqs)
phyloseq_new
tax_table(phyloseq_new)<-tax_table(physeq)
#phy_tree(phyloseq_new)<-phy_tree(physeq)
phyloseq_new

########### begin to merge tree
# tree <- read_tree("./phylogenetic_tree/tree.nwk")
# phy_tree(phyloseq_new) <- tree


Biostrings::writeXStringSet(rep.seqs,file ="rep.seqs.fasta",   format = "fasta")
phyloseq_new
print(phyloseq_new)
sample_names(phyloseq_new) <- gsub("-", "_",  sample_names(phyloseq_new) )
saveRDS(phyloseq_new, paste0(workpath,"/ampliseq_phyloseq.rds"))
# saveRDS(phyloseq_new, "~/Desktop/ampliseq_phyloseq.rds")







# #########nf-core phyloseq to biom
# library("biomformat")
# library("reltools")
# library(seqinr)
# OTU1 = t(as(otu_table(phyloseq_new), "matrix"))
# otu<-as(otu_table(OTU1,taxa_are_rows = FALSE),"matrix")
# otu_biom<-make_biom(data=otu)
#
# dir.create(paste0(set_work_path,'/output/pircust2/otu_biom'),recursive = TRUE)
# write_biom(otu_biom, paste0(set_work_path,"/output/pircust2/otu_biom/otu_biom.biom"))
# #########nf-corephyloseq to fasta
# #output fasta
# phyloseq_new.seq<-refseq(phyloseq_new)  %>% as.list()
# write.fasta(phyloseq_new.seq , names(phyloseq_new.seq),"uniqueSeqs.fasta")
# # df_<-as.data.frame(otu_table(phyloseq_new))


# print("???????运行完成")
# #########nf-core phyloseq to biom
library("biomformat")
#install.packages("remotes")
#remotes::install_github("DanielSprockett/reltools")
library("reltools")
library(seqinr)
OTU1 = t(as(otu_table(phyloseq_new), "matrix"))
otu<-as(otu_table(OTU1,taxa_are_rows = TRUE),"matrix") %>% t()
head(otu)
otu_biom<-make_biom(data=otu)

# dir.create(paste0(set_work_path,'/output/pircust2/otu_biom'),recursive = TRUE)
write_biom(otu_biom, paste0("otu_biom.biom"))
# #########nf-corephyloseq to fasta
# #output fasta
phyloseq_new.seq<-refseq(phyloseq_new)  %>% as.list()
print("开始写fasta文件")
write.fasta(phyloseq_new.seq , names(phyloseq_new.seq),paste0("uniqueSeqs.fasta"))
# print("完成fasta的书写文件")
# # df_<-as.data.frame(otu_table(phyloseq_new))
