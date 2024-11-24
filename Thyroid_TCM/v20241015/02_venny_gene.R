library(ggVennDiagram)
library(dplyr)
library(ggVennDiagram)
library(glue)
library(ggplot2)


# 1. 甲状腺所有相关基因
GeneCards_SearchResults_Thyroid_tumor <- read.csv("input/GeneCards-SearchResults_Thyroid_tumor.csv")
Disease_gene=GeneCards_SearchResults_Thyroid_tumor$Gene.Symbol

# 2. Load the data

process_tcm_data <- function(Filename, Disease_gene) {
  # 读取数据并处理
  DF.gene <- read.csv(Filename)
  DF.gene <- DF.gene %>% filter(!is.na(entrez_gene_symbol))
  TCM_genes <- unique(unlist(strsplit(DF.gene$entrez_gene_symbol, ";")))
  Type <- gsub(".*_|.csv", "", Filename)
  
  # 计算共同基因并保存
  Tcm_common_gene <- intersect(Disease_gene, TCM_genes)
  Tcm_common_gene_df <- data.frame(Com_Gene = Tcm_common_gene)
  write.csv(Tcm_common_gene_df, 
            glue("output/02_common_gene_{Type}.csv"),
            quote = TRUE, 
            row.names = FALSE)
  
  # 创建Venn图
  TCM_venny <- ggVennDiagram(
    list(Disease_gene = Disease_gene, TCM_genes = TCM_genes),
    edge_size = 0.8,
    edge_lty = 1,
    label_alpha = 0,
    set_color = "darkgreen"
  ) +
    scale_fill_gradient(low = "white", high = "green") +
    ggtitle(glue("{Type}")) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    coord_flip()
  
  # 保存Venn图
  ggsave(
    plot = TCM_venny,
    filename = glue("output/02_Venny_{Type}.png"),
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # 返回一些可能有用的数据
  return(list(
    TCM_genes = TCM_genes,
    Common_genes = Tcm_common_gene,
    Venn_plot = TCM_venny
  ))
}


# zhebeim
Filename <- "output/01_TCM_CID_compound_Gene_ZBM.csv"
result.z <- process_tcm_data(Filename, Disease_gene)

# xaikucao
Filename <- "output/01_TCM_CID_compound_Gene_XKC.csv"
result.x <- process_tcm_data(Filename, Disease_gene)

# all
Filename <- "output/01_TCM_CID_compound_Gene_XYT.csv"
result.p <- process_tcm_data(Filename, Disease_gene)

result.z$Venn_plot
result.x$Venn_plot
result.p$Venn_plot
# 
# 

library(ggplot2)
library(patchwork)

# 假设 result.z$Venn_plot, result.x$Venn_plot, result.p$Venn_plot 是你的三个 Venn 图
p1 <- result.z$Venn_plot
p2 <- result.x$Venn_plot
p3 <- result.p$Venn_plot

# 使用 patchwork 组合三个图，并在左上角添加字母标签
combined_plot <- (p1 | p2 | p3) + plot_annotation(tag_levels = "A")  # 添加 A, B, C 标签

# 显示组合后的图像
print(combined_plot)

ggsave("output/02_combined_venn_plot.png", combined_plot, width = 15, height = 4, dpi = 300)

