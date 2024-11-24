library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
library(dplyr)
library(ggplot2)
library(glue)
plot_go_enrichment <- function(go_enrich, output_file = "GO_enrichment_ordered.png", width = 12, height = 20) {
  library(ggplot2)
  library(dplyr)
  # 准备绘图数据
  plot_data <- go_enrich@result
  plot_data$Count <- as.numeric(plot_data$Count)
  # 定义ONTOLOGY的顺序
  ontology_order <- c("BP", "CC", "MF")
  # 按ONTOLOGY分类并选择每个类型的top30，然后按指定顺序排序
  top30 <- plot_data %>%
    group_by(ONTOLOGY) %>%
    top_n(30, wt = -p.adjust) %>%
    arrange(ONTOLOGY, desc(Count)) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = ontology_order),
           Description = factor(Description, levels = rev(unique(Description))))
  # 创建图形
  p <- ggplot(top30, aes(x = Count, y = Description, fill = ONTOLOGY)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("BP" = "#FC8D62", "CC" = "#66C2A5", "MF" = "#BEBADA")) +
    labs(title = "The Most Enriched GO Terms",
         x = "Gene Number",
         y = "GO term") +
    theme_minimal() +
    theme(
      # axis.text.y = element_text(size = 6),
      axis.text.y = element_text(size = 12, hjust = 1),  # 增加Y轴标签字体大小
      axis.title.y = element_text(size = 12),  # 增加Y轴标题字体大小
      legend.title = element_blank(),
      strip.text.y = element_blank(),  # 隐藏facet标签
      strip.background = element_blank(),  # 移除facet标签背景
      panel.grid.major.y = element_blank()) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(expand = c(0, 0.5)) +  # 调整Y轴标签与图表的距离
    scale_x_continuous(expand = c(0, 0))  # 调整X轴与图表的距离
  if(is.na(height) ){ height=0.15*nrow(top30)+5}
  print(height)
  # 保存图形
  ggsave(output_file, plot = p, width = width, height = height, dpi = 300)
  # 返回ggplot对象，以便进一步修改或显示
  return(p)
}

plot_kegg_enrichment_plot <- function(kegg_enrich, output_file = "enrichment_plot.png", width = 15, height = 25, dpi = 300) {
  # 对数据进行排序和处理
  plot_data <- kegg_enrich@result %>%
    dplyr::filter(p.adjust<0.05) %>%
    arrange(category, subcategory, Description) %>%
    group_by(category, subcategory) %>%
    top_n(10, wt = -p.adjust) %>%
    ungroup() %>%
    mutate(
      category = factor(category, levels = unique(category)),
      subcategory = factor(subcategory, levels = unique(subcategory)),
      Description = factor(Description, levels = rev(unique(Description))),
      category_label = ifelse(!duplicated(category), as.character(category), "")
    )
  # 创建图形
  p <- ggplot(plot_data, aes(x = Count, y = Description)) +
    geom_point(aes(size = -log10(p.adjust), color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    facet_grid(subcategory ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = "富集通路分析",
         x = "基因数量",
         size = "-log10(p.adjust)",
         color = "p.adjust") +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, size = 12, hjust = 0, vjust = 1),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.text.y = element_text(size = 12, hjust = 1),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right"
    )
  # 保存图形
  plot_data
  if(is.na(height) ){ height=0.15*nrow(plot_data)+5}
  print(height)
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  # 返回图形对象
  return(p)
}


analyze_tcm_data <- function(Filename) {
  # 创建输出目录
  dir.create("output", showWarnings = FALSE)
  
  # 读取数据
  DF_gene <- read.csv(Filename)
  genes <- unique(DF_gene$Com_Gene)
  Type <- gsub(".*_|.csv", "", Filename)
  
  # 将基因symbol转换为Entrez ID
  gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  # GO富集分析
  go_enrich <- enrichGO(gene = gene_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  
  # KEGG富集分析
  kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
  
  # 绘制GO富集分析图
  plot_go_enrichment(go_enrich, output_file = glue("output/03_my_go_enrichment_{Type}.png"), width = 15, height = NA)
  
  # 绘制KEGG富集分析图
  plot_kegg_enrichment_plot(kegg_enrich, output_file = glue("output/03_my_kegg_enrichment_plot_{Type}.png"), width = 18, height = NA, dpi = 400)
  
  # 保存富集分析结果
  write.csv(as.data.frame(go_enrich), glue("output/03_go_enrichment_results_{Type}.csv"), row.names = FALSE)
  write.csv(as.data.frame(kegg_enrich), glue("output/03_kegg_enrichment_results_{Type}.csv"), row.names = FALSE)
  
  # 返回结果
  return(list(go_enrich = go_enrich, kegg_enrich = kegg_enrich))
}

# 保留原有的plot_go_enrichment和plot_kegg_enrichment_plot函数

# xiakucao
Filename <- "output/02_common_gene_xiakucao.csv"
results <- analyze_tcm_data(Filename)

# zhebeimu
Filename <- "output/02_common_gene_zhebeimu.csv"
results <- analyze_tcm_data(Filename)





# # GO富集分析结果可视化
# # 点图
# pdf("GO_dotplot.pdf", width = 10, height = 8)
# dotplot(go_enrich, showCategory = 20) + 
#   theme_bw() + 
#   theme(axis.text.y = element_text(size = 8))
# dev.off()
# 
# # 条形图
# pdf("GO_barplot.pdf", width = 10, height = 8)
# barplot(go_enrich, showCategory = 20) + 
#   theme_bw() + 
#   theme(axis.text.y = element_text(size = 8))
# dev.off()
# 
# # 网络图
# pdf("GO_cnetplot.pdf", width = 12, height = 10)
# cnetplot(go_enrich, categorySize = "pvalue", foldChange = gene_entrez$ENTREZID)
# dev.off()
# 
# # KEGG富集分析结果可视化
# # 点图
# pdf("KEGG_dotplot.pdf", width = 10, height = 8)
# dotplot(kegg_enrich, showCategory = 20) + 
#   theme_bw() + 
#   theme(axis.text.y = element_text(size = 8))
# dev.off()
# 
# # 条形图
# pdf("KEGG_barplot.pdf", width = 10, height = 8)
# barplot(kegg_enrich, showCategory = 20) + 
#   theme_bw() + 
#   theme(axis.text.y = element_text(size = 8))
# dev.off()
# 
# # KEGG通路图
# for (pathway in kegg_enrich$ID) {
#   pathview(gene.data = gene_entrez$ENTREZID, 
#            pathway.id = pathway, 
#            species = "hsa", 
#            out.suffix = "pathview")
# }

# 保存富集分析结果


# # 疾病 富集分析
# enrich_do <- enrichDGN(gene = gene_entrez$ENTREZID,
#                        pvalueCutoff = 0.05, 
#                        pAdjustMethod = "BH",
#                        # universe = names(geneList),
#                        minGSSize = 10,
#                        maxGSSize = 500,
#                        qvalueCutoff = 0.2,
#                        readable = FALSE)
