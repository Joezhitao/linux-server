library(Seurat)
library(pheatmap)
library(readxl)
library(dplyr)

rm(list = ls())
setwd("/home/lin/c_group/")
# ============ 参数设置 ============
# 文件路径
seurat_file <- "/home/lin/c_group/CD4+T_ident.rds"
function_file <- "/home/lin/c_group/function.xlsx"
gmt_file <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"

# 分组设置
group_by <- "orig.ident"
group_order <- c()  # 留空则使用默认顺序

# 选择的Hallmark基因集（留空则不使用）
selected_hallmarks <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_INFLAMMATORY_RESPONSE", 
                        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                        "HALLMARK_IL2_STAT5_SIGNALING",
                        "HALLMARK_TGF_BETA_SIGNALING",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_ALLOGRAFT_REJECTION",
                        "HALLMARK_APOPTOSIS",
                        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "HALLMARK_GLYCOLYSIS")

# 输出设置
output_dir <- "./"
output_prefix <- "CD4"
heatmap_file_zscore <- paste0(output_prefix, "_Pathway_Heatmap_Zscore.pdf")
heatmap_file_raw <- paste0(output_prefix, "_Pathway_Heatmap_Raw.pdf")
show_raw_scores <- TRUE

# 热图设置
heatmap_width <- 10
heatmap_height <- 8
color_low <- "#6495ED"
color_mid <- "white"
color_high <- "#548B54"
# ====================================

# 读取GMT文件函数
read_gmt <- function(file_path) {
  lines <- readLines(file_path)
  gene_sets <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    gene_sets[[fields[1]]] <- fields[3:length(fields)]
  }
  return(gene_sets)
}

# 读取数据
pbmc <- readRDS(seurat_file)
function_genes <- lapply(read_excel(function_file), function(x) x[!is.na(x)])

# 添加Hallmark基因集
if(length(selected_hallmarks) > 0 && file.exists(gmt_file)) {
  hallmark_sets <- read_gmt(gmt_file)
  cat("\nProcessing Hallmark gene sets...\n")
  
  for(set_name in selected_hallmarks) {
    if(set_name %in% names(hallmark_sets)) {
      function_genes[[set_name]] <- hallmark_sets[[set_name]]
      cat(paste0("✓ Found: ", set_name, " (", length(hallmark_sets[[set_name]]), " genes)\n"))
    } else {
      cat(paste0("✗ Not found: ", set_name, "\n"))
    }
  }
  cat("----------------------------------------\n")
}

# 计算基因集评分
scores_list <- list()
cat("\nMatching genes to expression data...\n")

for(name in names(function_genes)) {
  genes <- function_genes[[name]]
  all_genes <- rownames(pbmc)
  
  # 基因匹配
  matched_genes <- genes[genes %in% all_genes]
  if(length(matched_genes) == 0) {
    for(gene in genes) {
      idx <- which(toupper(all_genes) == toupper(gene))
      if(length(idx) > 0) matched_genes <- c(matched_genes, all_genes[idx[1]])
    }
  }
  matched_genes <- unique(matched_genes)
  
  cat(paste0(name, ": ", length(matched_genes), "/", length(genes), " genes matched\n"))
  
  if(length(matched_genes) >= 3) {
    col_name <- paste0("Score_", gsub("HALLMARK_", "", name))
    pbmc <- PercentageFeatureSet(pbmc, features = matched_genes, col.name = col_name)
    
    scores_list[[name]] <- pbmc@meta.data %>%
      group_by(!!sym(group_by)) %>%
      summarise(score = mean(!!sym(col_name))) %>%
      mutate(pathway = name)
  } else {
    cat("  → Skipped (< 3 genes)\n")
  }
}

cat(paste0("\nPathways in heatmap: ", length(scores_list), "\n"))
cat("----------------------------------------\n")

# 创建热图矩阵
score_df <- bind_rows(scores_list)
heatmap_matrix <- score_df %>%
  tidyr::pivot_wider(names_from = all_of(group_by), values_from = score) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# 设置列顺序
if(length(group_order) > 0) {
  heatmap_matrix <- heatmap_matrix[, intersect(group_order, colnames(heatmap_matrix))]
}

# 生成热图
plot_heatmap <- function(mat, filename, scale = TRUE, show_numbers = FALSE) {
  if(scale) {
    mat <- t(scale(t(mat)))
    mat[is.na(mat)] <- 0
  }
  
  pdf(paste0(output_dir, filename), width = heatmap_width, height = heatmap_height)
  pheatmap(mat,
           color = colorRampPalette(c(color_low, color_mid, color_high))(100),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = NA,
           display_numbers = show_numbers,
           number_format = "%.3f",
           fontsize_row = 10,
           fontsize_col = 10)
  dev.off()
}

# 生成热图
plot_heatmap(heatmap_matrix, heatmap_file_zscore, scale = TRUE)
if(show_raw_scores) {
  plot_heatmap(heatmap_matrix, heatmap_file_raw, scale = FALSE, show_numbers = TRUE)
}

# 输出总结
cat("\nAnalysis completed!\n")
cat("Files saved:\n")
cat(paste0("- ", output_dir, heatmap_file_zscore, "\n"))
if(show_raw_scores) cat(paste0("- ", output_dir, heatmap_file_raw, "\n"))
