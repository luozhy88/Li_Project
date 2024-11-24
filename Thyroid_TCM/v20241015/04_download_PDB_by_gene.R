

# # 安装并加载必要的包
# if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
# if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
# 

library(jsonlite)
library(glue)
library(httr)
library(jsonlite)
library(httr)
library(jsonlite)
library(glue)
library(dplyr)

gene_to_pdb <- function(Gene) {
  # 构建查询JSON
  query_json <- glue::glue('
  {{
    "query": {{
      "type": "terminal",
      "label": "full_text",
      "service": "full_text",
      "parameters": {{
        "value": {Gene_json}
      }}
    }},
    "return_type": "entry",
    "request_options": {{
      "paginate": {{
        "start": 0,
        "rows": 25
      }},
      "results_content_type": [
        "experimental"
      ],
      "sort": [
        {{
          "sort_by": "score",
          "direction": "desc"
        }}
      ],
      "scoring_strategy": "combined"
    }}
  }}
  ', Gene_json = jsonlite::toJSON(Gene, auto_unbox = TRUE))
  
  # 发送POST请求到PDB API
  response <- POST(
    url = "https://search.rcsb.org/rcsbsearch/v2/query",
    body = query_json,
    encode = "json",
    add_headers("Content-Type" = "application/json")
  )
  
  # 创建必要的目录
  dir.create("input/Gene_to_PDB_csv", showWarnings = FALSE, recursive = TRUE)
  dir.create("input/PDB", showWarnings = FALSE, recursive = TRUE)
  
  # 检查响应状态
  if (status_code(response) == 200) {
    # 解析JSON响应
    result <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    result_df <- result[["result_set"]]
    result_df$Gene <- Gene
    
    # 过滤和保存结果
    result_df <- result_df %>% dplyr::filter(score > 0.9)
    write.csv(result_df, file = glue::glue("input/Gene_to_PDB_csv/{Gene}_to_PDB.csv"), row.names = FALSE, quote = TRUE)
    
    # 打印结果
    cat("查询到的PDB ID:\n")
    cat(paste(result_df$identifier, collapse = ","), "\n\n")
    
    if (nrow(result_df) > 0 && length(result_df$identifier) > 0) {
      # 创建下载链接
      result_df$download_url <- paste0("https://files.rcsb.org/download/", result_df$identifier, ".pdb")
      
      cat("找到", nrow(result_df), "个PDB结构:\n")
      print(result_df[, c("identifier", "download_url")])
      
      # 下载文件
      for (url in result_df$download_url) {
        dest_file <- file.path("input/PDB", basename(url))
        if (!file.exists(dest_file)) {
          tryCatch({
            download.file(url, destfile = dest_file, mode = "wb")
            cat("下载完成:", basename(url), "\n")
          }, error = function(e) {
            cat("下载失败:", basename(url), "- 错误:", conditionMessage(e), "\n")
          })
        } else {
          cat("文件已存在，跳过下载:", basename(url), "\n")
        }
      }
    } else {
      cat("未找到符合条件的PDB结构。\n")
    }
  } else {
    cat("API请求失败，状态码:", status_code(response), "\n")
  }
}


# 调用函数
gene_to_pdb("TP53")
gene_to_pdb("BRCA1")

# 将Genes生成PDB
tbl <-
  list.files(path = "output", pattern = "^02_common_gene", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(read.csv) %>% dplyr::distinct()

tbl$Com_Gene %>% purrr::walk(gene_to_pdb)

# 合并所有生成PDB的信息
read_and_convert <- function(file) {
  df <- read.csv(file)
  df$identifier <- as.character(df$identifier)  # 将identifier列转换为字符型
  df
}

tbl_pdb <- list.files(path = "input/Gene_to_PDB_csv", pattern = ".*csv", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(read_and_convert) %>% dplyr::distinct()
write.csv(tbl_pdb, "input/04_tidy_Gene_PDB.csv", row.names = FALSE, quote = TRUE)

