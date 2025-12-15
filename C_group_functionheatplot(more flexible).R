#############################################################
# 单细胞通路热图参数化脚本 (v3: 指定通路筛选版)
# 功能：多格式基因集读取(xlsx/csv/gmt)、跨物种转化、指定通路筛选、评分、热图
# 亮点：支持 selected_hallmarks 参数，仅分析感兴趣的通路
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载包与环境设置 ------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readxl)
  library(pheatmap)
  library(biomaRt)
  library(digest) # 用于生成安全列名
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 辅助功能函数 ----------------------------------------------------------

# --- 函数 A: 智能读取基因集 (支持 xlsx, csv, gmt) ---
read_geneset_file <- function(file_path) {
  if (is.null(file_path) || file_path == "") return(list())
  if (!file.exists(file_path)) { warning(paste("文件不存在:", file_path)); return(list()) }
  
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext == "gmt") {
    lines <- readLines(file_path)
    gene_sets <- list()
    for (line in lines) {
      fields <- strsplit(line, "\t")[[1]]
      gene_sets[[fields[1]]] <- fields[3:length(fields)]
    }
    return(gene_sets)
    
  } else if (ext == "xlsx") {
    df <- read_excel(file_path)
    return(lapply(df, function(x) x[!is.na(x)]))
    
  } else if (ext == "csv") {
    df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    return(lapply(df, function(x) x[!is.na(x) & x != ""]))
    
  } else {
    warning("不支持的文件格式，请使用 .gmt, .xlsx 或 .csv")
    return(list())
  }
}

# --- 函数 B: BiomaRt 跨物种转换核心 ---
get_ortholog_map <- function(genes, from_sp, to_sp) {
  if (from_sp == to_sp) return(NULL) # 同物种无需转换
  
  cat(sprintf("  [BiomaRt] 启动转换: %s -> %s\n", from_sp, to_sp))
  
  # 数据库配置
  db_cfg <- list(
    human = c("hsapiens_gene_ensembl", "hgnc_symbol"),
    mouse = c("mmusculus_gene_ensembl", "mgi_symbol"),
    rat   = c("rnorvegicus_gene_ensembl", "rgd_symbol")
  )
  
  if (!all(c(from_sp, to_sp) %in% names(db_cfg))) stop("仅支持 human, mouse, rat")
  
  host_url <- "https://dec2021.archive.ensembl.org" # 使用归档镜像更稳定
  
  tryCatch({
    mart_from <- useMart("ensembl", dataset = db_cfg[[from_sp]][1], host = host_url)
    mart_to   <- useMart("ensembl", dataset = db_cfg[[to_sp]][1], host = host_url)
  }, error = function(e) stop("BiomaRt连接失败，请检查网络或更换Host"))
  
  # 分批查询
  batch_size <- 500
  batches <- split(genes, ceiling(seq_along(genes) / batch_size))
  res_list <- list()
  
  cat("  [BiomaRt] 正在查询同源基因...\n")
  for (i in seq_along(batches)) {
    try({
      res <- getLDS(attributes = db_cfg[[from_sp]][2], filters = db_cfg[[from_sp]][2],
                    values = batches[[i]], mart = mart_from,
                    attributesL = db_cfg[[to_sp]][2], martL = mart_to, uniqueRows = TRUE)
      res_list[[i]] <- res
    }, silent = TRUE)
  }
  
  final_map <- do.call(rbind, res_list)
  if (is.null(final_map)) stop("未找到匹配基因，请检查物种名称是否正确")
  
  colnames(final_map) <- c("Original", "Target")
  # 返回字典列表: Key=Target(Human), Value=Original(Rat)
  return(split(final_map$Original, final_map$Target))
}

# 3. 主分析函数 ------------------------------------------------------------

run_pathway_heatmap <- function(
    seurat_obj,
    species_data = "rat",       # 数据物种
    species_set = "human",      # 基因集物种
    geneset_file = NULL,        # 路径 (xlsx/csv)
    gmt_file = NULL,            # 路径 (gmt)
    selected_pathways = NULL,   # 【新增】指定要分析的通路列表，NULL则分析所有
    subset_mode = FALSE,        # 是否亚群分析
    subset_col = NULL,          # 亚群分组列
    subset_idents = NULL,       # 亚群ID
    group_by = "orig.ident",    # 热图展示的分组
    output_dir = "./heatmap_results",
    file_prefix = "Pathway_Activity"
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始通路热图分析流程\n")
  
  # --- 1. 数据预处理 ---
  work_obj <- seurat_obj
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("亚群模式下需指定 subset_col 和 subset_idents")
    cat(sprintf("  [INFO] 模式: 亚群分析 (筛选 %s: %s)\n", subset_col, paste(subset_idents, collapse=",")))
    Idents(work_obj) <- subset_col
    work_obj <- subset(work_obj, idents = subset_idents)
  } else {
    cat("  [INFO] 模式: 全局分析\n")
  }
  
  # --- 2. 读取并合并基因集 ---
  all_sets <- c(read_geneset_file(geneset_file), read_geneset_file(gmt_file))
  if (length(all_sets) == 0) stop("未读取到任何基因集，请检查文件路径")
  
  # --- 2.5 筛选指定通路 (新增核心逻辑) ---
  if (!is.null(selected_pathways) && length(selected_pathways) > 0) {
    cat(sprintf("  [FILTER] 正在筛选 %d 个指定通路...\n", length(selected_pathways)))
    # 查找匹配的通路 (区分大小写，严格匹配)
    matched_names <- intersect(names(all_sets), selected_pathways)
    
    if (length(matched_names) == 0) {
      stop("错误：指定的 selected_pathways 在文件/GMT中一个都没找到！请检查拼写。")
    }
    
    if (length(matched_names) < length(selected_pathways)) {
      missing <- setdiff(selected_pathways, names(all_sets))
      cat(sprintf("  [WARN] 以下 %d 个通路未找到，将被忽略:\n    %s\n", length(missing), paste(missing, collapse=", ")))
    }
    
    all_sets <- all_sets[matched_names]
    cat(sprintf("  [INFO] 最终保留 %d 个通路进行分析\n", length(all_sets)))
  } else {
    cat(sprintf("  [INFO] 未指定筛选，共分析 %d 个通路\n", length(all_sets)))
  }
  
  # --- 3. 基因映射准备 ---
  all_genes <- rownames(work_obj)
  gene_dict <- list()
  
  if (species_data != species_set) {
    gene_dict <- get_ortholog_map(all_genes, species_data, species_set)
    cat(sprintf("  [INFO] 跨物种映射完成，字典包含 %d 个目标基因\n", length(gene_dict)))
  }
  
  # --- 4. 计算评分 (AddModuleScore) ---
  cat("  [INFO] 正在计算通路评分...\n")
  valid_pws <- c()
  
  for (pw_name in names(all_sets)) {
    target_genes <- all_sets[[pw_name]]
    
    # 获取数据中存在的基因
    if (species_data == species_set) {
      matched <- intersect(target_genes, all_genes)
    } else {
      # 异物种：查字典 (Target -> Original)
      matched <- unlist(gene_dict[intersect(target_genes, names(gene_dict))])
    }
    matched <- unique(matched)
    
    if (length(matched) >= 3) {
      valid_pws <- c(valid_pws, pw_name)
      safe_name <- paste0("S_", digest(pw_name, algo="crc32")) # 生成短且安全的列名
      
      work_obj <- AddModuleScore(work_obj, features = list(matched), name = safe_name, seed = 1)
      # Seurat会自动加 "1" 后缀，需重命名回通路名
      work_obj@meta.data[[pw_name]] <- work_obj@meta.data[[paste0(safe_name, "1")]]
    } else {
      cat(sprintf("    - 跳过通路 %s (匹配基因 < 3)\n", pw_name))
    }
  }
  
  if (length(valid_pws) == 0) {
    warning("没有通路的匹配基因数超过3个，无法绘图")
    return(NULL)
  }
  
  # --- 5. 生成矩阵并绘图 ---
  cat(sprintf("  [INFO] 正在绘制 %d 个通路的平均热图...\n", length(valid_pws)))
  
  # 聚合矩阵
  res_df <- work_obj@meta.data %>%
    group_by(!!sym(group_by)) %>%
    summarise(across(all_of(valid_pws), mean)) %>%
    column_to_rownames(group_by)
  
  mat_raw <- t(as.matrix(res_df))
  
  # 内部绘图函数
  plot_and_save <- function(mat, suffix, scale_row) {
    fname <- file.path(output_dir, paste0(file_prefix, suffix))
    
    # 数据修整
    plot_mat <- mat
    if (scale_row) {
      plot_mat <- t(scale(t(plot_mat))) # 按行标准化
      plot_mat[is.na(plot_mat)] <- 0
      plot_mat[plot_mat > 2.5] <- 2.5; plot_mat[plot_mat < -2.5] <- -2.5
    }
    
    pdf(fname, width = 10, height = max(5, length(valid_pws) * 0.35)) # 略微调高高度系数
    pheatmap(plot_mat,
             color = colorRampPalette(c("#4575B4", "white", "#7CCD7C"))(100),
             cluster_rows = TRUE, cluster_cols = FALSE,
             border_color = "white",
             display_numbers = !scale_row, number_format = "%.2f",
             main = paste("Pathway Activity", ifelse(scale_row, "(Z-score)", "(Raw)")))
    dev.off()
  }
  
  plot_and_save(mat_raw, "_Zscore.pdf", TRUE)
  plot_and_save(mat_raw, "_Raw.pdf", FALSE)
  
  # 导出数据
  write.csv(res_df, file.path(output_dir, paste0(file_prefix, "_Score_Matrix.csv")))
  
  cat("=================================================\n")
  cat(sprintf("分析完成！结果保存在: %s\n", output_dir))
  return(work_obj)
}


# 4. 参数定义与运行示例 ----------------------------------------------------

# --- 定义感兴趣的 Hallmark 基因集 (留空则使用全部) ---
selected_hallmarks <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V2", 
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_PEROXISOME",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_XENOBIOTIC_METABOLISM"
)

# 读取 Seurat 对象 (仅需运行一次)
seurat_obj <- readRDS("/home/lin/c_group/hep.rds")
levels(seurat_obj@meta.data$orig.ident)
# ==============================================
# 场景 A: 全局分析 (使用 selected_hallmarks 筛选)
# ==============================================
run_pathway_heatmap(
  seurat_obj = seurat_obj,
  species_data = "rat",
  species_set = "human",
  geneset_file = "/home/lin/c_group/function.xlsx",  # 可以同时读取 xlsx
  gmt_file = "", # 和 gmt
  selected_pathways = NULL,             # <--- 这里传入筛选列表
  group_by = "orig.ident",
  subset_mode = FALSE,
  output_dir = "/home/lin/c_group/New Folder/pathway_heatmap_results",
  file_prefix = "Global_Selected_Hallmarks"
)

# ==============================================
# 场景 B: 亚群分析 (使用 selected_hallmarks 筛选)
# ==============================================
run_pathway_heatmap(
  seurat_obj = seurat_obj,
  species_data = "rat",
  species_set = "human",
  geneset_file = "", # 不需要额外的 xlsx
  gmt_file = "/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt",
  selected_pathways = selected_hallmarks,             # <--- 这里传入筛选列表
  group_by = "cellcluster",
  subset_mode = TRUE,
  subset_col = "orig.ident",
  subset_idents = "Con",
  output_dir = "./pathway_heatmap_results",
  file_prefix = "Stressed_Subgroup_Selected_Hallmarks"
)
#/home/lin/B_group/h.all.v2025.1.Hs.symbols.gmt