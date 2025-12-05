library(Seurat)
library(pheatmap)
library(readxl)
library(dplyr)
library(tibble)
library(biomaRt)

rm(list = ls())
setwd("/home/lin/c_group/")

# ============ 1. 参数设置 ============
# 文件路径
seurat_file <- "/home/lin/c_group/hep_ident_filtered.rds"
function_file <- "/home/lin/c_group/function.xlsx"
gmt_file <- "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt"

# 物种设置 (必须准确填写)
# 可选值: "human", "mouse", "rat"
data_species <- "rat"         # 单细胞数据的物种
gene_set_species <- "human"   # function表或GMT的物种

# 分组与筛选
group_by <- "orig.ident"
group_order <- c()  

# 选择的Hallmark基因集（留空 "" 则不使用）
selected_hallmarks <- c("") 

# 输出设置
output_dir <- "./"
output_prefix <- "Hep_Pathway_BiomaRt"
heatmap_file_zscore <- paste0(output_prefix, "_Heatmap_Zscore.pdf")
heatmap_file_raw <- paste0(output_prefix, "_Heatmap_Raw.pdf")
show_raw_scores <- TRUE

# 热图外观
heatmap_width <- 10
heatmap_height <- 8
color_low <- "#6495ED"
color_mid <- "white"
color_high <- "#548B54"
# ====================================

# --- 辅助函数：读取GMT ---
read_gmt <- function(file_path) {
  lines <- readLines(file_path)
  gene_sets <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    gene_sets[[fields[1]]] <- fields[3:length(fields)]
  }
  return(gene_sets)
}

# --- 核心函数：BiomaRt 同源基因转换 ---
get_ortholog_map <- function(from_species, to_species, genes) {
  cat(sprintf("\n[BiomaRt] Starting ortholog conversion: %s -> %s\n", from_species, to_species))
  
  # 1. 定义数据集名称
  species_db <- list(
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl",
    "rat"   = "rnorvegicus_gene_ensembl"
  )
  
  species_attr <- list(
    "human" = "hgnc_symbol",
    "mouse" = "mgi_symbol",
    "rat"   = "rgd_symbol"
  )
  
  if(!from_species %in% names(species_db) || !to_species %in% names(species_db)) {
    stop("Species must be one of: human, mouse, rat")
  }
  
  dataset_from <- species_db[[from_species]]
  dataset_to   <- species_db[[to_species]]
  symbol_from  <- species_attr[[from_species]]
  symbol_to    <- species_attr[[to_species]]
  
  # 2. 连接 Ensembl (使用归档镜像增强稳定性)
  # 使用 2021 或 2022 的归档，通常比 www 稳定
  host_url <- "https://dec2021.archive.ensembl.org" 
  
  cat(sprintf("  Connecting to Ensembl (%s)...\n", host_url))
  
  tryCatch({
    mart_from <- useMart("ensembl", dataset = dataset_from, host = host_url)
    mart_to   <- useMart("ensembl", dataset = dataset_to, host = host_url)
  }, error = function(e) {
    stop("BiomaRt connection failed. Please check your internet or try a different host mirror.\nError: ", e$message)
  })
  
  # 3. 分批查询 (防止超时)
  batch_size <- 500
  batches <- split(genes, ceiling(seq_along(genes) / batch_size))
  results_list <- list()
  
  cat(sprintf("  Querying %d genes in %d batches...\n", length(genes), length(batches)))
  pb <- txtProgressBar(min = 0, max = length(batches), style = 3)
  
  for(i in seq_along(batches)) {
    batch_genes <- batches[[i]]
    
    # getLDS 是跨物种查询的核心函数
    # attributes: 输入物种的基因名
    # attributesL: 输出物种的基因名
    tryCatch({
      res <- getLDS(attributes = symbol_from,
                    filters = symbol_from,
                    values = batch_genes,
                    mart = mart_from,
                    attributesL = symbol_to,
                    martL = mart_to,
                    uniqueRows = TRUE)
      results_list[[i]] <- res
    }, error = function(e) {
      warning(paste("Batch", i, "failed:", e$message))
    })
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 4. 整理结果
  final_map <- do.call(rbind, results_list)
  
  if(is.null(final_map) || nrow(final_map) == 0) {
    stop("Conversion returned no results. Check if species names are correct.")
  }
  
  # 重命名列：第一列是原始基因，第二列是目标基因
  colnames(final_map) <- c("Original_Gene", "Target_Gene")
  
  # 过滤空值和去重
  final_map <- final_map %>%
    filter(Target_Gene != "" & !is.na(Target_Gene)) %>%
    distinct(Original_Gene, Target_Gene)
  
  cat(sprintf("  Done! Mapped %d genes.\n", nrow(final_map)))
  return(final_map)
}

# ================= 主流程 =================

# 1. 读取数据
cat(sprintf("[1/5] Reading Seurat object (%s)...\n", seurat_file))
pbmc <- readRDS(seurat_file)
all_genes_raw <- rownames(pbmc)

# 2. 准备基因映射字典
# 目标：建立一个字典，Key 是 "基因集里的基因(如Human)"，Value 是 "数据里的基因(如Rat)"
gene_dict <- list()

if(data_species == gene_set_species) {
  # 同物种：直接匹配
  cat("[2/5] Species match. Using direct mapping.\n")
  gene_dict <- setNames(all_genes_raw, all_genes_raw) # Key=Value
  
} else {
  # 异物种：使用 BiomaRt 转换
  cat(sprintf("[2/5] Species mismatch (%s vs %s). Using BiomaRt.\n", data_species, gene_set_species))
  
  # 获取转换表: Rat -> Human
  # 我们输入所有 Rat 基因，看它们对应什么 Human 基因
  orth_map <- get_ortholog_map(from_species = data_species, 
                               to_species = gene_set_species, 
                               genes = all_genes_raw)
  
  # 构建字典
  # 字典 Key: Human Gene (来自基因集)
  # 字典 Value: Rat Gene (来自单细胞数据)
  # 注意：可能存在多对一，这里我们为了 AddModuleScore，如果一个Human基因对应多个Rat基因，通常取任意一个或都保留比较复杂
  # 这里采用 Split 列表策略，稍后在循环中处理
  
  # 反向split: 此时 Key 是 Human Gene
  orth_list <- split(orth_map$Original_Gene, orth_map$Target_Gene)
  gene_dict <- orth_list
}

# 3. 读取基因集
cat("[3/5] Loading gene sets...\n")
final_gene_sets <- list()

# (A) 读取 Excel
if(file.exists(function_file)) {
  excel_data <- read_excel(function_file)
  excel_list <- lapply(excel_data, function(x) x[!is.na(x)])
  final_gene_sets <- c(final_gene_sets, excel_list)
}

# (B) 读取 GMT
if(length(selected_hallmarks) > 0 && selected_hallmarks[1] != "" && file.exists(gmt_file)) {
  hallmark_sets <- read_gmt(gmt_file)
  for(set_name in selected_hallmarks) {
    if(set_name %in% names(hallmark_sets)) {
      final_gene_sets[[set_name]] <- hallmark_sets[[set_name]]
    }
  }
}

# 4. 计算评分 (AddModuleScore)
cat("[4/5] Calculating pathway scores...\n")
valid_pathways <- c()

# 初始化结果DF
meta_data <- pbmc@meta.data
results_df <- data.frame(row.names = rownames(meta_data), group = meta_data[[group_by]])

for(pw_name in names(final_gene_sets)) {
  # 基因集里的基因 (例如 Human Genes)
  target_genes <- final_gene_sets[[pw_name]]
  
  # 找到数据中对应的基因 (例如 Rat Genes)
  matched_genes_in_data <- c()
  
  if(data_species == gene_set_species) {
    # 同物种直接交集
    matched_genes_in_data <- intersect(target_genes, all_genes_raw)
  } else {
    # 异物种查字典
    # gene_dict 是一个 list，names 是 Human Gene
    # 查找 target_genes 中哪些在字典里有 key
    found_targets <- intersect(target_genes, names(gene_dict))
    # 取出对应的 Rat Genes (unlist处理多对多)
    matched_genes_in_data <- unlist(gene_dict[found_targets])
    matched_genes_in_data <- unique(matched_genes_in_data)
  }
  
  match_rate <- length(matched_genes_in_data) # 这里分母其实不太好定，因为可能一对多
  cat(sprintf("  > %-30s: Found %d corresponding genes in data\n", pw_name, length(matched_genes_in_data)))
  
  if(length(matched_genes_in_data) >= 3) {
    valid_pathways <- c(valid_pathways, pw_name)
    
    # 生成安全列名
    safe_name <- paste0("Score_", digest::digest(pw_name, algo="crc32"))
    
    # 计算评分
    pbmc <- AddModuleScore(
      object = pbmc,
      features = list(matched_genes_in_data),
      name = safe_name,
      seed = 42
    )
    
    res_col <- paste0(safe_name, "1")
    results_df[[pw_name]] <- pbmc@meta.data[[res_col]]
    
  } else {
    cat("    [SKIPPED] < 3 genes matched\n")
  }
}

# 5. 绘图
cat("[5/5] Generating heatmaps...\n")

if(length(valid_pathways) > 0) {
  
  heatmap_matrix <- results_df %>%
    group_by(group) %>%
    summarise(across(all_of(valid_pathways), mean, .names = "{col}")) %>%
    column_to_rownames("group") %>%
    as.matrix() %>%
    t()
  
  if(length(group_order) > 0) {
    idx <- intersect(group_order, colnames(heatmap_matrix))
    if(length(idx) > 0) heatmap_matrix <- heatmap_matrix[, idx, drop=FALSE]
  }
  
  do_heatmap <- function(mat, filename, do_scale) {
    plot_mat <- mat
    if(do_scale) {
      plot_mat <- t(scale(t(plot_mat)))
      plot_mat[is.na(plot_mat)] <- 0
      plot_mat[plot_mat > 2.5] <- 2.5
      plot_mat[plot_mat < -2.5] <- -2.5
    }
    
    pdf(paste0(output_dir, filename), width = heatmap_width, height = heatmap_height)
    pheatmap(plot_mat,
             color = colorRampPalette(c(color_low, color_mid, color_high))(100),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             border_color = "white",
             display_numbers = !do_scale,
             number_format = "%.2f",
             main = ifelse(do_scale, "Pathway Activity (Z-score)", "Pathway Activity (Raw)"))
    dev.off()
  }
  
  do_heatmap(heatmap_matrix, heatmap_file_zscore, TRUE)
  if(show_raw_scores) do_heatmap(heatmap_matrix, heatmap_file_raw, FALSE)
  
  cat("\nSuccess! Heatmaps saved to:", output_dir, "\n")
  
} else {
  cat("\n[ERROR] No valid pathways found.\n")
}
