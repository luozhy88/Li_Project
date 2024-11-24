
library(dplyr)
library(ggplot2)

df_zhebeimu= read.csv("output/Docking_Vina_zhebeimu/06_Docking_Vina_result_ALL_scores.csv", header = T) %>% 
              dplyr::filter(vina_results<0 & grepl("^X",ligand ) & Model =="Model_1" ) %>% 
              dplyr::mutate(ID=paste0(ligand,"_",receptor),Type="ZBM" ) 

df_xiakucao= read.csv("output/Docking_Vina_xiakucao/06_Docking_Vina_result_ALL_scores.csv", header = T) %>% 
              dplyr::filter(vina_results<0 & grepl("^X",ligand ) & Model =="Model_1" ) %>%
              dplyr::mutate(ID=paste0(ligand,"_",receptor) ,Type="XKC" ) 
df_fufang= read.csv("output/Docking_Vina_fufang/Docking_Vina/06_Docking_Vina_result_ALL_scores.csv", header = T) %>% 
              dplyr::filter(vina_results<0 & grepl("^X",ligand ) & Model =="Model_1" ) %>%
              dplyr::mutate(ID=paste0(ligand,"_",receptor) ,Type="XYT" )

ALL_TCM= rbind(df_zhebeimu, df_xiakucao,df_fufang)


commonID = Reduce(intersect, list(df_zhebeimu$ID, df_xiakucao$ID, df_fufang$ID))


ALL_TCM_filter= ALL_TCM %>% dplyr::filter(ID %in% commonID)

# 使用dplyr包进行分组计算


ALL_TCM_filter <- ALL_TCM_filter %>%
  group_by(ID) %>%
  mutate(Average = mean(vina_results)) %>%
  ungroup()

avg_sorted = ALL_TCM_filter %>% arrange(Average) %>% slice_head(n = 51)

####################################add gene info###############################
Gene_to_protein <- read.csv("input/04_tidy_Gene_PDB.csv", header = T)
head(avg_sorted)
head(Gene_to_protein)
avg_sorted.gene=merge(avg_sorted,Gene_to_protein,by.x="receptor",by.y="identifier",all.x=T)
avg_sorted.gene$ID=gsub("^X","",avg_sorted.gene$ID)
######################################plot#######################################
# 首先合并ID和Gene为新的标签
avg_sorted.gene$label <- paste(avg_sorted.gene$ID, avg_sorted.gene$Gene, sep = " | ")

# 绘制图表
vina_plot <- ggplot(avg_sorted.gene, 
                    aes(x = vina_results, 
                        y = reorder(label, -vina_results), 
                        fill = Type)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           width = 0.7,
           alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Vina Energy Scores of Top 20 Compounds",
    subtitle = "Ranked by binding affinity",
    x = "Binding Energy (kcal/mol)",
    y = "",
    fill = "Compound Type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2")

vina_plot



ggsave("output/07.venny_barplot.png", 
       vina_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

# ggsave("output/07.venny_barplot.png", vina_plot, width = 18, height = 10, units = "cm")






