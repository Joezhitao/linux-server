library(Seurat)
library(dplyr)
library(readxl)

rm(list = ls())
# 设置工作目录
setwd("/home/lin/c_group/")

# 1. 读取数据
pbmc <- readRDS("/home/lin/c_group/hep_ident.rds")

# 2. 读取需要去除的基因列表
gene_df <- read_excel("/home/lin/c_group/gene.xlsx")
# 确保转换为字符向量，防止因为因子(factor)类型导致错误
genes_to_remove <- as.character(gene_df$gene) 

# --- 核心步骤开始 ---

# 3. 检查哪些待去除的基因实际上存在于你的Seurat对象中
# 这一步很重要，因为如果excel里的基因名和seurat里的不完全匹配（大小写、空格等），可能会报错或删不掉
existing_genes <- rownames(pbmc)

genes_found <- intersect(genes_to_remove, existing_genes)

cat("Excel中提供的基因数量:", length(genes_to_remove), "\n")
cat("Seurat对象中实际匹配到的基因数量:", length(genes_found), "\n")

if(length(genes_found) > 0) {
  # 4. 获取保留下来的基因列表 (总基因 - 要去除的基因)
  genes_to_keep <- setdiff(existing_genes, genes_found)
  
  # 5. 对Seurat对象进行取子集操作
  # subset函数通常用于筛选细胞，筛选基因我们可以直接通过矩阵索引或者使用subset功能的features参数(Seurat V3/V4/V5略有不同，最通用的方法是创建一个新的对象)
  
  # 方法 A (推荐，最稳健): 使用 subset 函数指定 features
  pbmc_subset <- subset(pbmc, features = genes_to_keep)
  
  # 检查结果
  cat("原始基因数:", length(existing_genes), "\n")
  cat("去除后基因数:", nrow(pbmc_subset), "\n")
  
  # 将处理后的对象重新赋值给 pbmc (如果你确定要覆盖)
  pbmc <- pbmc_subset
  
  # 6. 保存结果 (可选，建议保存为新文件名以防覆盖原数据)
  saveRDS(pbmc, "/home/lin/c_group/hep_ident_filtered.rds")
  cat("已保存去除基因后的数据到 hep_ident_filtered.rds\n")
  
} else {
  cat("警告: Excel中的基因在Seurat对象中一个都没找到，请检查基因名格式（如大小写）是否一致。\n")
}
