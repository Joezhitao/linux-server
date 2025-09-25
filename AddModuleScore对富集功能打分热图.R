# 加载必要的包
library(Seurat)
library(readxl)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(xlsx)

rm(list = ls())
setwd("/home/lin/B_group")

# ============ 参数设置区域 ============
# 数据文件路径
seurat_obj_path <- '/home/lin/B_group/Hepatocytes.RDS'
gene_list_path <- "/home/lin/B_group/gene3.xlsx"
gene_list_sheet <- 1

# 分组设置
group_by <- "seurat_clusters"  # 可选: "seurat_clusters", "group", 或其他meta.data中的列名

# 输出设置
output_prefix <- "Hepatocytes_function"

# 热图颜色设置
use_light_colors <- TRUE  # TRUE使用浅色调，FALSE使用深色调
# ====================================

# 读取Seurat对象
pbmc <- readRDS(seurat_obj_path)

# 确保选择的分组变量存在
if (!group_by %in% colnames(pbmc@meta.data)) {
  stop(paste("错误: 分组变量", group_by, "在Seurat对象的meta.data中不存在"))
}

# 读取包含功能基因列表的Excel文件
gene_list <- read.xlsx(gene_list_path, sheetIndex = gene_list_sheet, header = TRUE)

# 手动定义功能组映射 - 可根据您的数据调整
module_groups <- c(
  "Naïve" = "Differentiation",
  "Activation" = "Differentiation",
  "Exhaustion" = "Differentiation",
  "TCR_Signaling" = "Function",
  "Cytotoxicity" = "Function",
  "Cytokine" = "Function",
  "Chemokine" = "Function",
  "Senescence" = "Function",
  "NFKB_Signaling" = "Function",
  "Stress_response" = "Function",
  "MAPK_Signaling" = "Function",
  "Adhesion" = "Function",
  "IFN_Response" = "Function",
  "Oxidative_phosphorylation" = "Metabolism",
  "Glycolysis" = "Metabolism",
  "Fatty_acid_metabolism" = "Metabolism",
  "Pro_apoptosis" = "Apoptosis"
)

# 显示可用的列名
cat("Excel文件中的可用列名:\n")
print(colnames(gene_list))

# 创建功能模块的基因列表
module_gene_lists <- list()

# 遍历每个模块，提取非NA的基因
for (module in names(module_groups)) {
  if (module %in% colnames(gene_list)) {
    # 过滤掉NA值
    genes <- na.omit(gene_list[[module]])
    # 只保留在Seurat对象中存在的基因
    valid_genes <- genes[genes %in% rownames(pbmc)]
    # 如果有有效基因，则添加到列表中
    if (length(valid_genes) >= 3) {
      module_gene_lists[[module]] <- valid_genes
      print(paste("模块:", module, "有效基因数:", length(valid_genes)))
    } else {
      print(paste("警告: 模块", module, "有效基因数少于3，跳过"))
    }
  } else {
    print(paste("警告: 模块", module, "在Excel文件中不存在"))
  }
}

# 使用AddModuleScore计算每个功能模块的得分
for (i in seq_along(module_gene_lists)) {
  module_name <- names(module_gene_lists)[i]
  gene_set <- module_gene_lists[[i]]
  
  # 创建一个安全的名称，避免特殊字符
  safe_name <- gsub("[^a-zA-Z0-9]", "_", module_name)
  
  pbmc <- AddModuleScore(
    pbmc,
    features = list(gene_set),
    name = paste0(safe_name, "_score"),
    seed = 123
  )
}

# 获取所有模块得分的列名
score_cols <- grep("_score1$", colnames(pbmc@meta.data), value = TRUE)

# 创建模块名称的映射，从安全名称映射回原始名称
module_name_map <- setNames(
  names(module_gene_lists),
  paste0(gsub("[^a-zA-Z0-9]", "_", names(module_gene_lists)), "_score1")
)

# 提取得分矩阵
score_matrix <- pbmc@meta.data[, score_cols]

# 重命名列以匹配原始模块名称
colnames(score_matrix) <- module_name_map[colnames(score_matrix)]

# 按选择的分组变量分组并计算平均得分
cell_type_scores <- data.frame(cell_type = pbmc@meta.data[[group_by]])
cell_type_scores <- cbind(cell_type_scores, score_matrix)

# 计算每个细胞类型的平均得分
avg_scores <- cell_type_scores %>%
  group_by(cell_type) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()

# 将cell_type列设为行名
rownames(avg_scores) <- avg_scores$cell_type
avg_scores$cell_type <- NULL

# 对得分矩阵进行z-score标准化
scale_scores <- t(apply(avg_scores, 2, function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}))

# 确保所有值都是有限的
scale_scores[!is.finite(scale_scores)] <- 0

# 按照功能组对行进行排序
modules_order <- names(module_groups)
modules_order <- modules_order[modules_order %in% rownames(scale_scores)]
scale_scores <- scale_scores[modules_order, ]

# 定义功能组颜色
if (use_light_colors) {
  # 浅色调
  group_colors <- c(
    "Differentiation" = "#FFF2CC", # 浅黄色
    "Function" = "#D0E0E3",        # 浅蓝色
    "Metabolism" = "#E6E6E6",      # 浅灰色
    "Apoptosis" = "#F4CCCC"        # 浅粉色
  )
  
  # 定义热图颜色 - 浅色调
  neg_color <- "#A1C4E6"   # 浅蓝色
  mid_color <- "#F5F5F5"   # 浅灰白色
  pos_color <- "#F8BBD0"   # 浅粉红色
} else {
  # 深色调
  group_colors <- c(
    "Differentiation" = "#FFD966", # 深黄色
    "Function" = "#9BC2E6",        # 深蓝色
    "Metabolism" = "#BFBFBF",      # 深灰色
    "Apoptosis" = "#EA9999"        # 深粉色
  )
  
  # 定义热图颜色 - 深色调
  neg_color <- "#4472C4"   # 深蓝色
  mid_color <- "#FFFFFF"   # 白色
  pos_color <- "#ED7D31"   # 深橙色
}

# 创建功能组注释
module_groups_ordered <- module_groups[modules_order]
ha_row <- rowAnnotation(
  Function_Group = module_groups_ordered,
  col = list(Function_Group = group_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(0.5, "cm")
)

# 创建图例注释
legend_labels <- unique(module_groups_ordered)
ha_legend <- Legend(
  title = "Function_Group",
  labels = legend_labels,
  legend_gp = gpar(fill = group_colors[legend_labels]),
  title_gp = gpar(fontsize = 12),
  labels_gp = gpar(fontsize = 11)
)

# 创建颜色渐变
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c(neg_color, mid_color, pos_color)
)

# 绘制热图
ht <- Heatmap(
  scale_scores,
  name = "Signature score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  left_annotation = ha_row,
  heatmap_legend_param = list(
    title = "Signature score",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 11)
  ),
  rect_gp = gpar(col = "white", lwd = 0.5)
)

# 构建输出文件名
output_filename_base <- paste0(output_prefix, "_by_", group_by)

# 保存热图
pdf(paste0(output_filename_base, "_heatmap.pdf"), width = 10, height = 8)
draw(ht, annotation_legend_list = list(ha_legend))
dev.off()
cat(paste0("PDF热图已保存为: ", output_filename_base, "_heatmap.pdf\n"))

# 另外保存一个PNG版本
png(paste0(output_filename_base, "_heatmap.png"), width = 1200, height = 1000, res = 120)
draw(ht, annotation_legend_list = list(ha_legend))
dev.off()
cat(paste0("PNG热图已保存为: ", output_filename_base, "_heatmap.png\n"))
