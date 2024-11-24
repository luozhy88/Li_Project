library(tidyr)
library(dplyr)

process_tcm_data <- function(input_filename, pdb_filename = "input/04_tidy_Gene_PDB.csv") {
  # 读取TCM数据
  df <- read.csv(input_filename)
  df <- df %>% dplyr::filter(!is.na(entrez_gene_symbol))
  
  # 拆分entrez_gene_symbol列
  df_split <- df %>% 
    separate_rows(entrez_gene_symbol, sep = ";") %>%
    filter(entrez_gene_symbol != "")
  
  # 重置行名
  rownames(df_split) <- NULL
  
  # 读取和处理PDB数据
  X04_tidy_Gene_PDB <- read.csv(pdb_filename)
  x04_tidy_Gene_PDB <- X04_tidy_Gene_PDB %>%
    dplyr::group_by(Gene) %>%
    summarise(identifier = paste(unique(identifier), collapse = ";"))
  
  # 合并数据
  df_CID_Gene_PDB <- merge(df_split, x04_tidy_Gene_PDB, 
                           by.x = "entrez_gene_symbol", by.y = "Gene", all.x = TRUE)
  df_CID_Gene_PDB <- df_CID_Gene_PDB %>% dplyr::filter(!is.na(identifier))
  
  # 生成输出文件名并保存结果
  out_name <- gsub("01_", "05_Tidy_", input_filename)
  write.csv(df_CID_Gene_PDB, out_name, row.names = FALSE, quote = TRUE)
  
  return(df_CID_Gene_PDB)
}



# 使用默认的PDB文件名
result <- process_tcm_data("output/01_TCM_CID_compound_Gene_zhebeimu.csv")
result <- process_tcm_data("output/01_TCM_CID_compound_Gene_xiakucao.csv")
result <- process_tcm_data("output/01_TCM_CID_compound_Gene_prescription.csv")
