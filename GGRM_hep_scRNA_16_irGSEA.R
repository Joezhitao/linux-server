# ==============================================================================
# 0. 环境准备与工具函数定义
# ==============================================================================
library(Seurat)
library(irGSEA)
library(tidyverse)
library(ggplot2)

rm(list = ls())
setwd("/home/lin/GGRM_hep/result/")
if(!dir.exists("irGSEA_Result")) dir.create("irGSEA_Result")

# --- 核心：自定义 GMT 读取函数 (沿用你成功的方案) ---
# 这个函数不依赖 clusterProfiler，绝对稳
read_gmt_custom <- function(file_path) {
  if (!file.exists(file_path)) stop("GMT file not found: ", file_path)
  lines <- readLines(file_path)
  gene_sets <- list()
  for (line in lines) {
    # 按制表符分割：第1项是通路名，第2项是描述(跳过)，第3项开始是基因
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) > 2) {
      gene_sets[[fields[1]]] <- fields[3:length(fields)]
    }
  }
  return(gene_sets)
}

# ==============================================================================
# 1. 数据加载与标准化
# ==============================================================================
# 读取数据 -> 更新对象 -> 标准化(如有必要)
pbmc <- readRDS("/home/lin/GGRM_hep/result/cluster_optimization/filtered_both_threshold0.2.rds")
# 优先使用 cell_state，否则降级使用 seurat_clusters
Idents(pbmc) <- if ("cell_state" %in% colnames(pbmc@meta.data)) "cell_state" else "seurat_clusters"
cat("当前分析分组:", levels(Idents(pbmc)), "\n")

# ==============================================================================
# 2. 离线基因集构建 (使用自定义函数)
# ==============================================================================
gmt_path <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"

cat("正在读取 GMT 文件...\n")
gene_sets <- read_gmt_custom(gmt_path)

cat("本地 GMT 加载完毕，共", length(gene_sets), "条通路。\n")
# 检查一下第一条，确保格式正确 (应该是一个字符向量)
print(head(gene_sets[[1]])) 

# ==============================================================================
# 3. irGSEA 核心计算 (Scoring)
# ==============================================================================
cat("开始 irGSEA 评分 (约需几分钟)...\n")

# 重点：custom=TRUE, msigdb=FALSE, geneset传入我们刚读取的列表
pbmc_irGSEA <- irGSEA.score(
  object      = pbmc,
  assay       = "RNA",
  slot        = "data",
  geneset     = gene_sets,  # 传入自定义List
  custom      = TRUE,       # 开启自定义
  msigdb      = FALSE,      # 关闭下载
  method      = c("AUCell", "UCell", "singscore", "ssgsea"),
  species     = "Homo sapiens",
  kcdf        = 'Gaussian',
  ncores      = 20,
  seeds       = 123
)

saveRDS(pbmc_irGSEA, "irGSEA_Result/pbmc_irGSEA_scored.rds")
cat("评分对象已保存。\n")

# ==============================================================================
# 4. 差异整合分析 (Integration & RRA)
# ==============================================================================
cat("正在整合差异结果...\n")
options(future.globals.maxSize = 1000 * 1024^5)

result_dge <- irGSEA.integrate(
  object   = pbmc_irGSEA,
  group.by = "group",
  method   = c("AUCell", "UCell", "singscore", "ssgsea")
)

saveRDS(result_dge, "irGSEA_Result/irGSEA_DGE_result.rds")

# ==============================================================================
# 5. 优雅可视化 (Visualization)
# ==============================================================================
setwd("irGSEA_Result")

# 定义绘图辅助函数
save_plot <- function(plot_obj, filename, w=10, h=8) {
  pdf(filename, width = w, height = h)
  print(plot_obj)
  dev.off()
  cat("已保存:", filename, "\n")
}

# 5.1 概览图表
irGSEA.heatmap(result_dge, method = "RRA", top = 50) %>% save_plot("01_Global_Heatmap.pdf", 12, 12)
irGSEA.bubble(result_dge, method = "RRA", top = 10)  %>% save_plot("02_Bubble_Plot.pdf", 12, 8)
irGSEA.upset(result_dge, method = "RRA")             %>% save_plot("03_Upset_Plot.pdf", 10, 6)
irGSEA.barplot(result_dge, method = c("AUCell", "UCell", "singscore", "ssgsea")) %>% 
  save_plot("04_Barplot.pdf", 10, 6)

# 5.2 Top 通路精细绘图
# 尝试提取第一条显著通路
top_pathway <- tryCatch({
  result_dge$RRA$Name[1]
}, error = function(e) NA)

if (!is.na(top_pathway)) {
  cat("正在绘制 Top Pathway:", top_pathway, "\n")
  
  irGSEA.halfvlnplot(pbmc_irGSEA, method = "AUCell", show.geneset = top_pathway) %>% 
    save_plot(paste0("05_HalfVln_", top_pathway, ".pdf"), 8, 6)
  
  irGSEA.ridgeplot(pbmc_irGSEA, method = "AUCell", show.geneset = top_pathway) %>% 
    save_plot(paste0("06_Ridge_", top_pathway, ".pdf"), 8, 6)
  
  irGSEA.densityheatmap(pbmc_irGSEA, method = "AUCell", show.geneset = top_pathway) %>% 
    save_plot(paste0("07_DensityHeatmap_", top_pathway, ".pdf"), 8, 6)
}

cat("Analysis Completed Successfully! \n")
#密度图
scatterplot <- irGSEA.density.scatterplot(object = pbmc_irGSEA,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION",
                                          reduction = "umap")
scatterplot
