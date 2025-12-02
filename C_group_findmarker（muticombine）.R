# 加载必要的库
library(Seurat)
library(dplyr)
library(openxlsx)

rm(list = ls())

# 设置工作目录
setwd("/home/lin/c_group/")

# 读取之前处理好的单细胞数据
pbmc <- readRDS("/home/lin/c_group/hep_ident.rds")

# 检查可用的标识符
print(levels(pbmc@meta.data$orig.ident))

# 定义时间点和输出目录
timepoints <- c("Cd1", "Cd15", "Cd30")
output_dir <- "/home/lin/c_group/subpopulation_analysis/Hepatocytes/DiffGenes_LR_comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 创建一个空的数据框，用于存储所有比较的结果
all_comparisons <- data.frame()

# 循环处理每个时间点
for (tp in timepoints) {
  # 构建左右样本的标识
  left_ident <- paste0(tp, "L")
  right_ident <- paste0(tp, "R")
  
  cat("分析时间点:", tp, "- 比较", left_ident, "vs", right_ident, "\n")
  
  # 找差异基因
  markers <- FindMarkers(
    object = pbmc, 
    ident.1 = left_ident, 
    ident.2 = right_ident, 
    group.by = "orig.ident",
    test.use = "wilcox", 
    min.pct = 0.1,
    logfc.threshold = 0.25,
    only.pos = FALSE
  )
  
  # 添加基因名称和其他信息
  markers$gene <- rownames(markers)
  markers$comparison <- paste0(tp, "_L_vs_R")
  markers$cluster <- tp  # 将timepoint改名为cluster
  
  # 添加到总表
  all_comparisons <- bind_rows(all_comparisons, markers)
}

# 保存总表
write.xlsx(all_comparisons, file = file.path(output_dir, "All_LR_comparisons.xlsx"), rowNames = FALSE)

# 为每个时间点分别保存结果
for (tp in timepoints) {
  tp_data <- all_comparisons %>% filter(cluster == tp)
  if (nrow(tp_data) > 0) {
    write.xlsx(tp_data, file = file.path(output_dir, paste0(tp, "_LR_comparison.xlsx")), rowNames = FALSE)
  }
}

# 输出完成信息
cat("左右对比分析已完成，结果已保存到", output_dir, "\n")
