library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

setwd("/home/lin/B_group")

rm(list = ls())

# ============ 参数设置区域 ============
# 设置分组变量
group_by <- "seurat_clusters"  # 可以是 "group" 或 "seurat_clusters" 或其他元数据列

# 设置输出路径
output_dir <- "./"

# 是否显示原始分数
show_raw_scores <- TRUE

# 是否为每个类别生成单独热图
generate_category_heatmaps <- TRUE
# ====================================

# 读取已处理的Seurat对象
pbmc <- readRDS('/home/lin/B_group/Hepatocytes.RDS')

# 确保选择的分组变量存在
if (!(group_by %in% colnames(pbmc@meta.data))) {
  stop(paste("Error: Group variable", group_by, "does not exist in Seurat object metadata"))
}

# 读取GMT文件
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

# 添加自定义衰老基因集
senescence_genes <- c(
  "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3",
  "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4",
  "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1",
  "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1",
  "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7",
  "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2",
  "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A",
  "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA",
  "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2",
  "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU",
  "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1",
  "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2",
  "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2"
)
all_hallmark_sets[["HALLMARK_SENESCENCE"]] <- senescence_genes

# 添加自定义DDR基因集
ddr_genes <- c(
  "APLF", "APTX", "ASCC3", "DNTT", "LIG1", "LIG3", "LIG4", "MRE11A", "NBN", "NHEJ1",
  "PARG", "PARP1", "PARP3", "PARPBP", "PNKP", "POLB", "POLL", "POLM", "PRKDC", "RAD50",
  "RNF168", "RNF8", "TP53BP1", "XRCC1", "XRCC2", "XRCC3", "XRCC4", "XRCC5", "XRCC6", "UBE2A",
  "EXO1", "HMGB1", "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PCNA", "PMS1", "PMS2",
  "POLD1", "POLD2", "POLD3", "POLD4", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5", "RPA1",
  "RPA2", "RPA3", "RPA4", "ALKBH1", "ALKBH2", "ALKBH3", "APEX1", "APEX2", "APITD1", "ATM"
)
all_hallmark_sets[["CUSTOM_DDR"]] <- ddr_genes

# 添加增殖基因集
proliferation_genes <- c("CCND1", "CCNE1", "CCNE2", "CCNA2", "CCNB1", "MKI67")
all_hallmark_sets[["CUSTOM_PROLIFERATION"]] <- proliferation_genes

# 定义新的分组方式
damage_stress_response_pathways <- c(
  "HALLMARK_DNA_REPAIR",
  "CUSTOM_DDR",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

cell_state_fate_pathways <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_SENESCENCE",
  "CUSTOM_PROLIFERATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_KRAS_SIGNALING_UP"
)

metabolic_pathways <- c(
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_GLYCOLYSIS"
)

# 合并所有基因集
all_pathways <- c(damage_stress_response_pathways, cell_state_fate_pathways, metabolic_pathways)

# 创建一个数据框来存储每个基因集在每个组中的平均表达值，用于热图
heatmap_data <- data.frame(pathway = character(), 
                           group = character(), 
                           score = numeric(), 
                           category = character(),
                           stringsAsFactors = FALSE)

# 为每个通路进行分析
for (pathway in all_pathways) {
  # 确定该通路属于哪个类别
  if (pathway %in% damage_stress_response_pathways) {
    category <- "Damage & Stress Response"
  } else if (pathway %in% cell_state_fate_pathways) {
    category <- "Cell State & Fate"
  } else {
    category <- "Metabolic Reprogramming"
  }
  
  # 获取人类基因
  human_genes <- all_hallmark_sets[[pathway]]
  
  cat(paste("Analyzing", pathway, "...\n"))
  cat(paste("Number of human genes:", length(human_genes), "\n"))
  
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
  
  cat(paste("Number of genes found in dataset:", length(genes_in_data), "\n"))
  
  if (length(genes_in_data) < 3) {  # 降低阈值以适应增殖基因集
    cat(paste("Warning: Too few genes found in dataset, skipping", pathway, "\n\n"))
    next
  }
  
  # 计算每个细胞中基因集的表达百分比
  gene_name <- paste("Pathway", "_", gsub("HALLMARK_|CUSTOM_", "", pathway), sep = "")
  pbmc <- PercentageFeatureSet(pbmc, features = genes_in_data, col.name = gene_name)
  
  # 计算每个组的平均表达值，用于热图
  group_levels <- unique(pbmc@meta.data[[group_by]])
  for (group_val in group_levels) {
    mean_score <- mean(pbmc@meta.data[pbmc@meta.data[[group_by]] == group_val, gene_name])
    pathway_name <- gsub("HALLMARK_|CUSTOM_", "", pathway)
    heatmap_data <- rbind(heatmap_data, data.frame(pathway = pathway_name, 
                                                   group = as.character(group_val), 
                                                   score = mean_score,
                                                   category = category))
  }
}

# 生成热图
if (nrow(heatmap_data) > 0) {
  # 检查和处理数据中的异常值
  cat("Checking for outliers in heatmap data...\n")
  heatmap_data$score[is.na(heatmap_data$score)] <- 0
  heatmap_data$score[is.infinite(heatmap_data$score)] <- 0
  heatmap_data$score[is.nan(heatmap_data$score)] <- 0
  
  # 获取原始数据中的组顺序
  original_group_order <- unique(pbmc@meta.data[[group_by]])
  
  # 使用简单的方式创建热图矩阵
  pathways <- unique(heatmap_data$pathway)
  groups <- original_group_order
  
  # 创建一个空矩阵
  heatmap_matrix <- matrix(0, nrow = length(pathways), ncol = length(groups))
  rownames(heatmap_matrix) <- pathways
  colnames(heatmap_matrix) <- groups
  
  # 填充矩阵
  for (i in 1:nrow(heatmap_data)) {
    pathway <- heatmap_data$pathway[i]
    group <- heatmap_data$group[i]
    score <- heatmap_data$score[i]
    heatmap_matrix[pathway, group] <- score
  }
  
  # 创建注释数据框
  pathway_categories <- data.frame(
    category = rep("", length(pathways)),
    row.names = pathways,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(pathway_categories)) {
    pathway <- rownames(pathway_categories)[i]
    category_data <- heatmap_data[heatmap_data$pathway == pathway, ]
    if (nrow(category_data) > 0) {
      pathway_categories$category[i] <- category_data$category[1]
    }
  }
  
  # 按照category对行进行排序
  category_order <- c("Damage & Stress Response", "Cell State & Fate", "Metabolic Reprogramming")
  pathway_categories$category <- factor(pathway_categories$category, levels = category_order)
  
  # 对pathway_categories进行排序
  pathway_categories <- pathway_categories[order(pathway_categories$category), , drop = FALSE]
  
  # 对heatmap_matrix按照相同顺序排序
  heatmap_matrix <- heatmap_matrix[rownames(pathway_categories), , drop = FALSE]
  
  # 对数据进行z-score标准化
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
  
  # 检查标准化后的数据，替换任何NA/NaN/Inf值
  heatmap_matrix_scaled[is.na(heatmap_matrix_scaled)] <- 0
  heatmap_matrix_scaled[is.infinite(heatmap_matrix_scaled)] <- 0
  heatmap_matrix_scaled[is.nan(heatmap_matrix_scaled)] <- 0
  
  # 为不同类别设置颜色
  category_colors <- list(category = c("Damage & Stress Response" = "#E41A1C", 
                                       "Cell State & Fate" = "#377EB8", 
                                       "Metabolic Reprogramming" = "#4DAF4A"))
  
  # 设置热图颜色 - 使用blue-white-red配色方案
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # 生成热图
  heatmap_path <- paste0(output_dir, "Pathway_Activity_Heatmap_by_", group_by, ".png")
  png(heatmap_path, width = 12, height = 10, units = "in", res = 300)
  
  pheatmap(
    heatmap_matrix_scaled,
    color = color_palette,  # 使用blue-white-red配色
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 12,
    fontsize_col = 12,
    main = paste0("Pathway Activity in Hepatocytes by ", group_by),
    border_color = NA,
    cellwidth = 40, 
    cellheight = 30,
    annotation_row = pathway_categories,
    annotation_colors = category_colors,
    annotation_names_row = FALSE
  )
  
  dev.off()
  cat(paste0("Heatmap saved to: ", heatmap_path, "\n"))
  
  # 生成未标准化的热图（显示原始分数）
  if (show_raw_scores) {
    heatmap_raw_path <- paste0(output_dir, "Pathway_Activity_Heatmap_Raw_Scores_by_", group_by, ".png")
    png(heatmap_raw_path, width = 12, height = 10, units = "in", res = 300)
    
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(c("blue", "white", "red"))(100),  # 使用blue-white-red配色
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = 12,
      fontsize_col = 12,
      main = paste0("Pathway Activity in Hepatocytes (Raw Scores) by ", group_by),
      border_color = NA,
      cellwidth = 40, 
      cellheight = 30,
      annotation_row = pathway_categories,
      annotation_colors = category_colors,
      annotation_names_row = FALSE,
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 7
    )
    
    dev.off()
    cat(paste0("Raw scores heatmap saved to: ", heatmap_raw_path, "\n"))
  }
  
  # 为三个类别分别生成热图
  if (generate_category_heatmaps) {
    for (category_name in category_order) {
      # 筛选该类别的通路
      category_pathways <- rownames(pathway_categories)[pathway_categories$category == category_name]
      
      if (length(category_pathways) > 0) {
        # 提取该类别的数据
        category_matrix <- heatmap_matrix[category_pathways, , drop = FALSE]
        
        # 对数据进行z-score标准化
        category_matrix_scaled <- t(scale(t(category_matrix)))
        
        # 处理NA/NaN/Inf值
        category_matrix_scaled[is.na(category_matrix_scaled)] <- 0
        category_matrix_scaled[is.infinite(category_matrix_scaled)] <- 0
        category_matrix_scaled[is.nan(category_matrix_scaled)] <- 0
        
        # 生成该类别的热图
        category_name_file <- gsub(" & ", "_and_", category_name)
        category_heatmap_path <- paste0(output_dir, "Pathway_Activity_", category_name_file, "_by_", group_by, ".png")
        png(category_heatmap_path, width = 10, height = 8, units = "in", res = 300)
        
        pheatmap(
          category_matrix_scaled,
          color = colorRampPalette(c("blue", "white", "red"))(100),  # 使用blue-white-red配色
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          fontsize_row = 12,
          fontsize_col = 12,
          main = paste0(category_name, " Pathways in Hepatocytes by ", group_by),
          border_color = NA,
          cellwidth = 40, 
          cellheight = 40
        )
        
        dev.off()
        cat(paste0(category_name, " heatmap saved to: ", category_heatmap_path, "\n"))
      }
    }
  }
}

# 打印分析完成消息
cat("All pathway analysis and heatmap generation completed!\n")
