library(LinDA)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
library(dplyr)
library(kableExtra)
rm(list=ls())

# Create a directory named "LinDA"
BATCH="SRP151288"
REF="HC"
notREF="TC"
tax="Species"


path <- glue::glue("../output/H20/",BATCH,"_",paste(REF,notREF,sep  = "_vs_"),".",tax   )
dir.create(path, recursive = TRUE)




# Read the phyloseq object from a file and filter out taxa with zero counts
phy <- readRDS("../../Kraken/input/output/Batch2_phy_meta.rds")
# Meta=meta(phy)
phy



pseq = phyloseq::subset_samples(phy,grepl(BATCH,Batch )  &grepl( paste(REF,notREF,sep  = "|")  ,GROUP ))
pseq
Meta=meta(pseq)

#sample_data(pseq)$age <- scale(sample_data(pseq)$age)
###############################################################
############### empty values for ranks replaced with "undefined"##
###############################################################
# Extract the tax_table from the phyloseq object
taxa <- as.data.frame(tax_table(pseq))

# Replace empty or NA values with "undefined"
taxa[is.na(taxa)] <- "undefined"
taxa[taxa == ""] <- "undefined"

# Convert the modified data frame back to a tax_table format
tax_table_modified <- tax_table(as.matrix(taxa))

# Assign the modified tax_table back to the phyloseq object
tax_table(pseq) <- tax_table_modified
pseq.core <- core(pseq, detection = 0.000001, prevalence = .1)
pseq.core <- tax_glom(pseq.core, taxrank = tax)
pseq.core_Tax= phyloseq_to_df(pseq.core) %>% as.data.frame() 
rownames(pseq.core_Tax)=paste(pseq.core_Tax$OTU,pseq.core_Tax[,tax],sep = "_")
pseq.core_Tax=pseq.core_Tax[,-c(1:9)]%>% t() %>% as.data.frame()
Meta=meta(pseq.core) %>% dplyr::select("GROUP")

pseq.core_df= merge(Meta,pseq.core_Tax,by=0) %>% tibble::column_to_rownames(var="Row.names") %>% as.data.frame()

############################################read ML feature importance###########
ML_feature_importance = read.csv("../output/H20/SRP151288_HC_vs_TC.Species/H20_feature_importance.csv") #%>% head()
ML_feature_importance$Sum=rowSums(ML_feature_importance[,-1])
ML_feature_importance=ML_feature_importance %>% dplyr::arrange(desc(Sum)) %>% head(6)
ML_feature_importance$variable

df_boxplot_input = pseq.core_df[,c("GROUP",ML_feature_importance$variable)]
df_boxplot_input$GROUP = as.factor(df_boxplot_input$GROUP)
colnames(df_boxplot_input)=colnames(df_boxplot_input) %>% make.names()
##three groups
library(ggstatsplot)
library(ggpubr)

plot_list=list()
for (i in colnames(df_boxplot_input)[-1] ){
  # i= "OTU2584_s__[Ruminococcus]_torques"
  print(i)
  my_comparisons <- list( c("HC", "TC"))
  p <- ggboxplot(df_boxplot_input, x = "GROUP", y = i,
                 color = "GROUP", palette =c("#00AFBB",  "#FC4E07"),
                 add = "jitter")+ 
    labs(title = i, x = "",y = "" )+
    theme(axis.title.y=element_blank(),
          plot.title = element_text(size = 8)
          )
  # Add p-values comparing groups
  # Specify the comparisons you want
  p  = p + stat_compare_means(comparisons = my_comparisons) #+ # Add pairwise comparisons p-value
    # stat_compare_means(label.y = 50)                   # Add global p-value
  p
  plot_list[[i]] <- p

}

plot_box1 <- ggarrange(plotlist = plot_list, common.legend = TRUE)#重新排版
plot_box1
ggsave(glue::glue(path,"/",tax,"_boxplot.pdf"),plot_box1,width = 10,height = 8)






