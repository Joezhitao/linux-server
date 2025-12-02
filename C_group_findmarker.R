# 加载必要的库
library(Seurat)

# 清空环境
rm(list = ls())

# 读取之前处理好的单细胞数据
pbmc <- readRDS("/home/lin/c_group/subpopulation_analysis/Hepatocytes/Hepatocytes_processed.rds")

# 定义时间点和对应的样本标识
timepoints <- c("Cd1", "Cd15", "Cd30")
sides <- c("L", "R")
pbmc@meta.data$orig.ident
# 创建结果存储目录
output_dir <- "/home/lin/c_group/subpopulation_analysis/Hepatocytes/DiffGenes_LR_comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 循环处理每个时间点
for (tp in timepoints) {
  # 构建左右样本的标识
  left_ident <- paste0(tp, sides[1])
  right_ident <- paste0(tp, sides[2])
  
  cat("分析时间点:", tp, "- 比较", left_ident, "vs", right_ident, "\n")
  
  # 找差异基因
  markers <- FindMarkers(
    object = pbmc, 
    ident.1 = left_ident, 
    ident.2 = right_ident, 
    group.by = "orig.ident",
    test.use = "wilcox", # 使用Wilcoxon秩和检验，更适合单细胞数据
    min.pct = 0.1,
    thresh.use = 0.25,
    only.pos = FALSE
  )
  
  # 添加基因名称列并排序
  markers$gene <- rownames(markers)
  markers <- markers[order(markers$p_val), ]
  
  # 保存结果
  output_file <- file.path(output_dir, paste0(tp, "_L_vs_R_markers.csv"))
  write.csv(markers, output_file, row.names = FALSE)
  
  cat("结果已保存到:", output_file, "\n\n")
}

# 输出完成信息
cat("所有时间点的左右对比分析已完成!\n")

# 保存更新后的Seurat对象（如果需要）
saveRDS(pbmc, "/home/lin/c_group/subpopulation_analysis/Hepatocytes/Hepatocytes_processed.rds")
