#############################################################
# 单细胞代谢分析参数化脚本 (v7: TopN参数化版)
# 功能：基因转换、代谢评分、亚群分析、桑葚图、基因集导出
# 更新：增加 top_n 参数，可自定义展示通路的数量，若超过总数则自动调整
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(biomaRt)
  library(scMetabolism)
  library(networkD3)
  library(viridis)
  library(AUCell)
  library(GSEABase)
  library(stringr)
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 核心功能函数定义 ------------------------------------------------------

# 内部函数：基因转换
convert_genes_internal <- function(obj, species) {
  expr_matrix <- tryCatch({
    as.matrix(GetAssayData(obj, layer = "counts"))
  }, error = function(e) {
    as.matrix(GetAssayData(obj, slot = "counts")) 
  })
  
  genes <- rownames(expr_matrix)
  species <- tolower(species)
  
  if (species == "human") {
    cat("  [INFO] 检测为人类数据，仅标准化基因名为大写...\n")
    rownames(expr_matrix) <- toupper(rownames(expr_matrix))
    new_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = obj@meta.data)
    new_obj@reductions <- obj@reductions
    Idents(new_obj) <- Idents(obj)
    return(new_obj)
  }
  
  cat(sprintf("  [INFO] 开始基因同源转换 (%s -> Human)...\n", species))
  if (species == "rat") {
    dataset_name <- "rnorvegicus_gene_ensembl"; symbol_attr <- "rgd_symbol"
  } else if (species == "mouse") {
    dataset_name <- "mmusculus_gene_ensembl"; symbol_attr <- "mgi_symbol"
  } else { stop("不支持的物种") }
  
  tryCatch({
    host_url <- "https://dec2021.archive.ensembl.org"
    input_mart <- useMart("ensembl", dataset = dataset_name, host = host_url)
    human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
  }, error = function(e) { stop("Ensembl 连接失败，请检查网络或更换 Host") })
  
  batch_size <- 1000
  batches <- split(genes, ceiling(seq_along(genes) / batch_size))
  gene_map_list <- list()
  
  for (i in seq_along(batches)) {
    try({
      res <- getLDS(attributes = symbol_attr, filters = symbol_attr, values = batches[[i]], 
                    mart = input_mart, attributesL = "hgnc_symbol", martL = human_mart, uniqueRows = TRUE)
      gene_map_list[[i]] <- res
    }, silent = TRUE)
  }
  
  gene_map <- do.call(rbind, gene_map_list)
  colnames(gene_map) <- c("Original_Gene", "Human_Gene")
  gene_map <- gene_map[gene_map$Human_Gene != "" & !duplicated(gene_map$Human_Gene), ]
  gene_map <- gene_map[gene_map$Original_Gene %in% rownames(expr_matrix), ]
  
  cat(sprintf("  [INFO] 成功转换基因: %d / %d\n", nrow(gene_map), length(genes)))
  
  sub_matrix <- expr_matrix[gene_map$Original_Gene, ]
  rownames(sub_matrix) <- gene_map$Human_Gene
  
  new_obj <- CreateSeuratObject(counts = sub_matrix, meta.data = obj@meta.data)
  Idents(new_obj) <- Idents(obj)
  return(new_obj)
}

# 内部函数：代谢评分
run_metabolism_internal <- function(obj) {
  cat("  [INFO] 正在计算 AUCell 代谢评分...\n")
  countexp <- as.matrix(GetAssayData(obj, layer = "counts"))
  gmtFile <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  geneSets <- getGmt(gmtFile)
  
  cells_rankings <- AUCell_buildRankings(countexp, plotStats = F, verbose = FALSE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, verbose = FALSE)
  
  auc_matrix <- getAUC(cells_AUC)
  auc_matrix <- auc_matrix[, colnames(obj), drop = FALSE]
  
  obj[["METABOLISM"]] <- CreateAssayObject(data = as.matrix(auc_matrix))
  return(obj)
}

# 辅助函数：清洗通路名称
clean_pathway_name <- function(x) {
  x <- toupper(x) 
  x <- gsub("[[:punct:]]", "", x) 
  x <- gsub(" ", "", x) 
  return(x)
}

# 3. 主函数：代谢分析流程 ------------------------------------------------------

run_metabolism_analysis <- function(
    seurat_obj,                  
    species = "rat",             
    group_by = "seurat_clusters",
    subset_mode = FALSE,         
    subset_col = NULL,           
    subset_idents = NULL,        
    output_dir = "./metabolism_result", 
    file_prefix = "Analysis",
    top_n = 20  # <--- 新增参数：默认为20
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始单细胞代谢分析流程\n")
  
  # --- 1. 数据预处理与筛选 ---
  final_obj <- seurat_obj
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("错误: 请设置 subset 参数！")
    cat(sprintf("模式: 亚群分析 (筛选 %s == %s)\n", subset_col, paste(subset_idents, collapse=",")))
    Idents(seurat_obj) <- subset_col
    final_obj <- subset(seurat_obj, idents = subset_idents)
  } else {
    cat("模式: 全局分析\n")
  }
  
  if (!group_by %in% colnames(final_obj@meta.data)) stop(paste("分组列不存在:", group_by))
  
  # --- 2. 基因转换 ---
  human_obj <- convert_genes_internal(final_obj, species)
  
  # --- 3. 代谢评分 ---
  human_obj <- run_metabolism_internal(human_obj)
  
  # --- 4. 绘制桑葚图 ---
  cat("  [INFO] 正在计算通路差异并绘制桑葚图...\n")
  met_data <- GetAssayData(human_obj, assay = "METABOLISM", layer = "data")
  idents <- human_obj@meta.data[[group_by]]
  
  avg_df <- data.frame(t(met_data))
  avg_df$group <- idents
  avg_res <- avg_df %>% group_by(group) %>% summarise(across(everything(), mean, .names = "{col}"))
  avg_mat <- t(avg_res[,-1])
  colnames(avg_mat) <- avg_res$group
  
  # --- 核心修改：筛选 Top N 通路逻辑 ---
  pathway_sd <- apply(avg_mat, 1, sd)
  
  # 计算实际可用的通路数量
  total_available <- length(pathway_sd)
  real_top_n <- min(top_n, total_available) # 取较小值，防止越界
  
  cat(sprintf("  [INFO] 请求 Top %d，实际可用 %d，将展示 Top %d 通路...\n", top_n, total_available, real_top_n))
  
  top_pathways <- names(sort(pathway_sd, decreasing = TRUE))[1:real_top_n]
  
  plot_data <- as.data.frame(avg_mat[top_pathways, , drop = FALSE])
  plot_data$Pathway <- rownames(plot_data)
  
  links <- plot_data %>%
    pivot_longer(cols = -Pathway, names_to = "TargetGroup", values_to = "Score") %>%
    filter(Score > quantile(Score, 0.25)) %>%
    mutate(Value = Score * 100)
  
  nodes <- data.frame(name = unique(c(links$Pathway, links$TargetGroup)))
  nodes$group <- ifelse(nodes$name %in% links$Pathway, "Pathway", "Group")
  links$IDsource <- match(links$Pathway, nodes$name) - 1
  links$IDtarget <- match(links$TargetGroup, nodes$name) - 1
  
  my_color <- 'd3.scaleOrdinal() .domain(["Pathway", "Group"]) .range(["#69b3a2", "#404080"])'
  
  # 动态调整高度，如果通路很多，增加图片高度
  plot_height <- max(600, real_top_n * 25) 
  
  sankey <- sankeyNetwork(
    Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
    Value = "Value", NodeID = "name", NodeGroup = "group", LinkGroup = "Pathway",
    colourScale = my_color, units = "Score", fontSize = 14, nodeWidth = 30,
    nodePadding = 10, sinksRight = FALSE, height = plot_height, width = 900
  )
  
  # --- 5. 保存结果 ---
  cat("  [INFO] 正在导出结果文件...\n")
  
  out_html <- file.path(output_dir, paste0(file_prefix, "_Sankey_Top", real_top_n, ".html"))
  out_rds <- file.path(output_dir, paste0(file_prefix, "_Metabolism.rds"))
  out_csv_all <- file.path(output_dir, paste0(file_prefix, "_Score_Mean.csv"))
  out_csv_genes <- file.path(output_dir, paste0(file_prefix, "_Pathway_Genes_Top", real_top_n, ".csv"))
  
  saveNetwork(sankey, file = out_html)
  saveRDS(human_obj, out_rds)
  
  all_pathways_df <- as.data.frame(avg_mat)
  all_pathways_df$Pathway <- rownames(all_pathways_df)
  all_pathways_df <- all_pathways_df %>% dplyr::select(Pathway, everything()) 
  write.csv(all_pathways_df, out_csv_all, row.names = FALSE)
  
  # (4) 导出基因集
  cat(sprintf("  [INFO] 正在解析并导出 Top %d 通路基因...\n", real_top_n))
  
  gmtFile <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  raw_lines <- readLines(gmtFile)
  gmt_list <- list()
  
  clean_top_dict <- setNames(top_pathways, clean_pathway_name(top_pathways))
  
  for(line in raw_lines) {
    parts <- strsplit(line, "\t")[[1]]
    p_name_raw <- parts[1]
    p_genes <- parts[-c(1,2)] 
    p_genes <- p_genes[p_genes != ""] 
    
    p_name_clean <- clean_pathway_name(p_name_raw)
    
    if(p_name_clean %in% names(clean_top_dict)) {
      original_name <- clean_top_dict[[p_name_clean]]
      gmt_list[[original_name]] <- p_genes
    }
  }
  
  if(length(gmt_list) > 0) {
    max_len <- max(sapply(gmt_list, length))
    pad_vec <- function(x, n) { c(x, rep("", n - length(x))) }
    padded_list <- lapply(gmt_list, pad_vec, n = max_len)
    genes_df <- data.frame(padded_list, check.names = FALSE)
    write.csv(genes_df, out_csv_genes, row.names = FALSE, quote = FALSE)
    cat(sprintf("  [SUCCESS] 成功导出 %d 条通路的基因信息！\n", length(gmt_list)))
  } else {
    warning("依然没有匹配到通路，请检查 scMetabolism 版本或 GMT 文件内容。")
  }
  
  cat("=================================================\n")
  cat(sprintf("分析完成！结果保存在: %s\n", output_dir))
  cat(sprintf("1. 桑葚图: %s\n", basename(out_html)))
  cat(sprintf("2. 评分均值表: %s\n", basename(out_csv_all)))
  cat(sprintf("3. Top%d通路基因集: %s\n", real_top_n, basename(out_csv_genes)))
  cat(sprintf("4. Seurat对象: %s\n", basename(out_rds)))
  
  return(human_obj)
}

# 4. 使用示例 ---------------------------------------------------------------

# 请确保这里读取的是您的 Seurat 对象路径
seurat_obj <- readRDS("/home/lin/c_group/hep.rds")
levels(seurat_obj@meta.data$cellcluster)
# ==============================================
# 场景 A: 全局分析 (自定义 Top 30)
# ==============================================
# 如果总通路只有 25 条，脚本会自动只展示 25 条
res_global <- run_metabolism_analysis(
  seurat_obj = seurat_obj,
  species = "rat",                  
  group_by = "cellcluster",
  subset_mode = FALSE,
  top_n = 40,                         # <--- 示例：改为查看前30条
  output_dir = "/home/lin/c_group/metabolism_results",
  file_prefix = "Rat_Global_Analysis"
)

# ==============================================
# 场景 B: 亚群分析 (使用默认 Top 20)
# ==============================================
res_subset <- run_metabolism_analysis(
  seurat_obj = seurat_obj,
  species = "rat",
  subset_mode = TRUE,               
  subset_col = "cellcluster",   
  subset_idents = "Regenerative-Unit-Hep", 
  group_by = "orig.ident",
  top_n = 20,                         # <--- 示例：显式设置为20（或不写该行，默认即为20）
  output_dir = "/home/lin/c_group/metabolism_results/stress",
  file_prefix = "Rat_Regenerative_Subgroup"
)
