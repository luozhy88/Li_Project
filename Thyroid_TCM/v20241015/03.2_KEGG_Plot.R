
library(dplyr)
rm(list=ls())
###############################################GO#################################
tbl <-
  list.files(path = "output", recursive = T, pattern = "03_go_enrichment_results_.*", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(read.csv)

plot_go_enrichment <- function(go_enrich, my_ontology="BP",output_file = "GO_enrichment_ordered.png", width = 12, height = 20) {
  # go_enrich=tbl
  library(ggplot2)

  # 准备绘图数据
  # browser()
  # my_ontology <- "BP"
  plot_data <- subset(go_enrich, ONTOLOGY == my_ontology)
  plot_data$Count <- as.numeric(plot_data$Count)
  # 定义Type的顺序
  Type_order <- c("ZBM", "XKC", "XYT")
  # 按Type分类并选择每个类型的top30，然后按指定顺序排序
  top30 <- plot_data %>%group_by(Type) %>%top_n(30, wt = -p.adjust) %>% arrange(Type, desc(Count)) %>%
    ungroup() #%>%
  Feature.selected=top30$Description
  
  #筛选表plot_data中的数据，只保留Feature.selected中的数据
  plot_data.selected=plot_data %>% dplyr::filter(Description %in% Feature.selected)
  table(plot_data.selected$Description)
  plot_data.selected2=plot_data.selected %>% mutate(Type = factor(Type, levels = Type_order),
           Description = factor(Description, levels = rev(unique(Description))))
  # 创建图形
  p <- ggplot(plot_data.selected, aes(x = Count, y = Description, fill = Type)) +
    # geom_bar(stat = "identity") +
    geom_bar(stat = "identity", 
             position = "dodge", 
             width = 0.7,            # 调整柱子宽度
             alpha = 0.8) +          # 添加透明度
    scale_fill_manual(values = c("ZBM" = "#FC8D62", "XKC" = "#66C2A5", "XYT" = "#BEBADA")) +
    labs(title = glue::glue("The Most Enriched GO Terms in {my_ontology}"),
         x = "Gene Number",
         y = "") +
    theme_minimal() +
    theme(
      # axis.text.y = element_text(size = 6),
      axis.text.y = element_text(size = 12, hjust = 1),  # 增加Y轴标签字体大小
      axis.title.y = element_text(size = 12),  # 增加Y轴标题字体大小
      legend.title = element_blank(),
      strip.text.y = element_blank(),  # 隐藏facet标签
      strip.background = element_blank(),  # 移除facet标签背景
      panel.grid.major.y = element_blank()) +
    # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(expand = c(0, 0.5),labels = scales::label_wrap(60)) +  # 调整Y轴标签与图表的距离
    scale_x_continuous(expand = c(0, 0))  # 调整X轴与图表的距离
  p
  
  
  if(is.na(height) ){ height=0.15*nrow(plot_data.selected)+5}
  print(height)
  # 保存图形
  ggsave(output_file, plot = p, width = width, height = height, dpi = 300)
  # 返回ggplot对象，以便进一步修改或显示
  return(p)
}
p1=plot_go_enrichment(tbl, my_ontology = "BP", output_file = "output/03.2_GO_enrichment_ordered_BP.png", width = 16, height = 18)
p2=plot_go_enrichment(tbl, my_ontology = "MF", output_file = "output/03.2_GO_enrichment_ordered_MF.png", width = 16, height = 18)
p3=plot_go_enrichment(tbl, my_ontology = "CC", output_file = "output/03.2_GO_enrichment_ordered_CC.png", width = 16, height = 18)

library(ggplot2)
library(patchwork)

# 使用 patchwork 组合三个图，并在左上角添加字母标签
# combined_plot <- (p1 | p2 | p3) + plot_annotation(tag_levels = "A")  # 添加 A, B, C 标签
# 基本方法
combined_plot <- (p1 | p2 | p3) + 
  plot_layout(guides = "collect") +   # 收集并共享图例
  plot_annotation(tag_levels = "A")   # 添加 A, B, C 标签

# 显示组合后的图像
print(combined_plot)

ggsave("output/03.2_combined_GO_plot.png", combined_plot, width = 30, height = 18, dpi = 300)


rm(list=ls())
###############################################KEGG#################################
tbl <-
  list.files(path = "output", recursive = T, pattern = "03_kegg_enrichment_results_.*", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(read.csv)

plot_kegg_enrichment_plot <- function(kegg_enrich, output_file = "enrichment_plot.png", width = 15, height = 25, dpi = 300) {
  # 对数据进行排序和处理
  # kegg_enrich=tbl
  plot_data <- kegg_enrich %>%
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
    # geom_point(aes(size = -log10(p.adjust), color = p.adjust, shape = Type)) +
    geom_point(aes( color = p.adjust, shape = Type)) +
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
  # p
  # 保存图形
  plot_data
  if(is.na(height) ){ height=0.15*nrow(plot_data)+5}
  print(height)
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  # 返回图形对象
  return(p)
}

plot_kegg_enrichment_plot(tbl, output_file = "output/03.2_kegg_enrichment_plot_ALL.png", width = 15, height = 25, dpi = 300)
