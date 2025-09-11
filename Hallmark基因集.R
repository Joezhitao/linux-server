library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("/home/lin/B_group")

rm(list = ls())

# 读取已处理的Seurat对象
pbmc <- readRDS('/home/lin/B_group/Hepatocytes.RDS')

# 设置输出文件路径
filepath <- "./"

# 读取下载的GMT文件
read_gmt <- function(file_path) {
  lines <- readLines(file_path)
  gene_sets <- list()
  
  for (line in lines) {
    fields <- unlist(strsplit(line, "\t"))
    gene_set_name <- fields[1]
    gene_set_desc <- fields[2]
    genes <- fields[3:length(fields)]
    gene_sets[[gene_set_name]] <- genes
  }
  
  return(gene_sets)
}

# 读取Hallmark基因集
file_path <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"
all_hallmark_sets <- read_gmt(file_path)

# 需要分析的基因集
needed_pathways <- c(
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

# 由于biomaRt服务器问题，我们使用一个更简单的方法：
# 1. 大多数基因符号在人类和大鼠之间是相似的，只是大小写不同
# 2. 我们将尝试直接匹配（忽略大小写）

# 为每个通路进行分析
for (pathway in needed_pathways) {
  # 获取人类基因
  human_genes <- all_hallmark_sets[[pathway]]
  
  cat(paste("分析", pathway, "...\n"))
  cat(paste("人类基因数量:", length(human_genes), "\n"))
  
  # 获取Seurat对象中的所有基因
  all_genes_in_data <- rownames(pbmc)
  
  # 尝试直接匹配基因名称（忽略大小写）
  matched_genes <- c()
  for (gene in human_genes) {
    # 检查基因名称是否直接存在（精确匹配）
    if (gene %in% all_genes_in_data) {
      matched_genes <- c(matched_genes, gene)
    } else {
      # 尝试忽略大小写匹配
      idx <- which(toupper(all_genes_in_data) == toupper(gene))
      if (length(idx) > 0) {
        matched_genes <- c(matched_genes, all_genes_in_data[idx[1]])
      }
    }
  }
  
  # 去除重复项
  genes_in_data <- unique(matched_genes)
  
  cat(paste("在数据集中找到的基因数量:", length(genes_in_data), "\n"))
  
  if (length(genes_in_data) < 5) {
    cat(paste("警告: 在数据集中找到的基因太少，跳过", pathway, "\n\n"))
    next
  }
  
  # 打印找到的一些基因，以便验证
  cat("找到的部分基因示例: ", paste(head(genes_in_data, 5), collapse=", "), "\n")
  
  # 计算每个细胞中基因集的表达百分比
  gene_name <- paste("肝细胞组间", "_", gsub("HALLMARK_", "", pathway), sep = "")
  pbmc <- PercentageFeatureSet(pbmc, features = genes_in_data, col.name = gene_name)
  
  # 准备组间比较
  factors <- levels(pbmc@meta.data$group)
  combinations <- combn(factors, 2) 
  combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  
  # 设置颜色
  my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')
  
  # 生成小提琴图
  path <- paste(filepath, "组间横向比较_", gene_name, ".png", sep = "")
  p.percentage <- ggviolin(pbmc@meta.data, x = "group", y = gene_name,
                           color = "group", add = 'mean_sd', fill = 'group',
                           add.params = list(color = "black")) + 
    stat_compare_means(comparisons = combinations_list, label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = '', y = paste(gsub("HALLMARK_", "", pathway), "Score")) +
    theme(legend.position = "none")  # 相当于NoLegend()
  
  ggsave(path, width = 8, height = 8)
  
  # 生成UMAP图
  path2 <- paste(filepath, "UMAP_", gene_name, ".png", sep = "")
  p2 <- FeaturePlot(pbmc, gene_name)
  ggsave(path2, width = 6, height = 6)
  
  cat(paste("完成", pathway, "的分析和图表生成\n\n"))
}

# 打印分析完成消息
cat("所有基因集分析完成!\n")
