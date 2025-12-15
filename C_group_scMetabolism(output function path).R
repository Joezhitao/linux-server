#############################################################
# 指定通路基因提取工具
# 功能：读取指定Excel表中的通路名，导出对应的基因集
# 输入：/home/lin/c_group/pathy.xlsx
# 输出：CSV格式的基因矩阵
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(scMetabolism)
})

# 2. 参数设置 --------------------------------------------------------------
input_excel <- "/home/lin/c_group/pathy.xlsx"        # 输入文件路径
output_csv  <- "/home/lin/c_group/Selected_Pathway_Genes.csv" # 输出文件路径
target_col  <- "function"                            # Excel中存放通路名的列名

# 3. 核心辅助函数：名称清洗 (用于模糊匹配) ---------------------------------
# 将 "Fatty.acid.biosynthesis" 和 "Fatty acid biosynthesis" 统一转换为 "FATTYACIDBIOSYNTHESIS"
clean_name <- function(x) {
  x <- toupper(x)
  x <- gsub("[[:punct:]]", "", x) # 去除标点（包括点号）
  x <- gsub(" ", "", x)           # 去除空格
  return(x)
}

# 4. 主逻辑 ----------------------------------------------------------------

cat(">>> 开始提取基因集...\n")

# (1) 读取您的目标列表
if(!file.exists(input_excel)) stop("找不到输入文件，请检查路径！")
target_df <- read_excel(input_excel)

if(!target_col %in% colnames(target_df)) stop(paste("列名", target_col, "在Excel中不存在！"))
my_pathways <- unique(target_df[[target_col]]) # 去重，防止重复提取
cat(sprintf("  [INFO] 待提取通路数量: %d\n", length(my_pathways)))

# (2) 读取 scMetabolism 的数据库 (GMT文件)
gmtFile <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
if(gmtFile == "") stop("未找到 scMetabolism 包数据，请确认已安装该包。")
raw_lines <- readLines(gmtFile)

# (3) 建立匹配字典
# 我们创建一个列表，键(Key)是清洗后的名字，值(Value)是基因向量
db_list <- list()
db_clean_names <- c()

for(line in raw_lines) {
  parts <- strsplit(line, "\t")[[1]]
  p_name_raw <- parts[1]       # 原始通路名
  p_genes <- parts[-c(1,2)]    # 基因列表 (跳过名字和描述)
  p_genes <- p_genes[p_genes != ""] 
  
  # 存入字典：Key为清洗后的名字
  clean_key <- clean_name(p_name_raw)
  db_list[[clean_key]] <- p_genes
}

# (4) 匹配与提取
result_list <- list()
matched_count <- 0

for(p in my_pathways) {
  # 清洗您的输入名以便去数据库查找
  query_key <- clean_name(p)
  
  if(query_key %in% names(db_list)) {
    # 如果找到了，就用您Excel里的名字作为列名，内容是数据库里的基因
    result_list[[p]] <- db_list[[query_key]]
    matched_count <- matched_count + 1
  } else {
    warning(paste("  [WARN] 未找到通路:", p))
  }
}

# (5) 格式化与导出
if(length(result_list) > 0) {
  # 找出最长的基因列表长度，用于补齐
  max_len <- max(sapply(result_list, length))
  
  # 补齐短的列表，使其等长
  padded_list <- lapply(result_list, function(x) {
    c(x, rep("", max_len - length(x)))
  })
  
  # 转为 Dataframe
  final_df <- data.frame(padded_list, check.names = FALSE)
  
  # 写入 CSV
  write.csv(final_df, output_csv, row.names = FALSE, quote = FALSE)
  
  cat("=================================================\n")
  cat(sprintf("  [SUCCESS] 提取完成！\n"))
  cat(sprintf("  成功匹配: %d / %d\n", matched_count, length(my_pathways)))
  cat(sprintf("  结果已保存至: %s\n", output_csv))
  cat("=================================================\n")
  
} else {
  cat("  [ERROR] 没有匹配到任何通路，请检查Excel中的名字拼写。\n")
}
