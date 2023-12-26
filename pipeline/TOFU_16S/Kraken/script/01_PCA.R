
rm(list=ls())
library(mixOmics)
phy= readRDS("../input/output/Batch2_phy_meta.rds")

ps_merge<-phy
data.filter <- otu_table(ps_merge) %>% data.frame()
data.filter <- data.filter +1
pca.result <- mixOmics::pca(t(data.filter), logratio = "CLR",scale=T)
Metadata <- ps_merge@sam_data %>% data.frame()
if(!file.exists("../output/PCA")){
  dir.create("../output/PCA",recursive = T)
}


output_sample<-"exploration_all"

pdf("../output/PCA/exploration_all.pdf",width = 12,height = 10)
p=plotIndiv(pca.result,
          group = Metadata$GROUP,
          # cex=Metadata$Batch,
          pch = as.numeric(factor(Metadata$Batch) ),
          pch.levels = Metadata$Batch,
          title = "PCA plot",
          style = "graphics",
          ellipse = T,
          legend = T,
          legend.title = "Category",
          point.lwd = 2,
          abline = T,
          ind.names = Metadata$Batch,
          # point.lwd=1,
          col.per.group = ggsci::pal_lancet("lanonc")(3)[1:nlevels(as.factor(Metadata$GROUP))]
)
dev.off()
