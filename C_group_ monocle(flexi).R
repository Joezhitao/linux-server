#############################################################
# Monocle2 拟时序分析参数化脚本 (修正版 v6)
# 功能：Seurat转Monocle2、轨迹构建、BEAM分析、热图
# 修复点：
#   1. 按照 orderCells 文档，增加 reverse_traj 参数
#   2. 若发现拟时序方向反了（起点变终点），只需将 reverse_traj = TRUE 即可
# 日期：2023-10-27
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 ----------------------------------------------------------
suppressPackageStartupMessages({
  library(monocle)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(patchwork)
  library(ggplot2)
  library(readxl)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)

# 2. 核心功能函数定义 ------------------------------------------------------

run_monocle2_analysis <- function(
    seurat_obj,                  
    group_by = "seurat_clusters", # 轨迹图着色依据
    subset_mode = FALSE,          # 是否开启亚群模式
    subset_col = NULL,            # 亚群筛选列
    subset_idents = NULL,         # 亚群筛选ID
    
    # --- [关键修改] 反转轨迹参数 ---
    reverse_traj = NULL,          # 是否反转拟时序 (TRUE/FALSE/NULL)
    # 如果跑出来的图发现起点是反的，设为 TRUE 即可
    
    gene_file = NULL,             # 指定基因Excel文件路径
    min_expr_prop = 0.05,         # 排序基因过滤阈值
    output_dir = "./monocle_result", 
    file_prefix = "Traj",
    cores = 4,
    
    # --- 高级分析开关 ---
    run_pseudotime_heatmap = TRUE, # 是否画随时间变化的热图
    run_branch_heatmap = TRUE,     # 是否画分支热图 (BEAM)
    heatmap_top_n = 50             # 热图展示基因数
) {
  
  # --- 0. 初始化 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat(sprintf("开始 Monocle2 分析: %s\n", file_prefix))
  
  # --- 1. 数据筛选 ---
  work_obj <- seurat_obj
  if (subset_mode) {
    if (is.null(subset_col) || is.null(subset_idents)) stop("错误: 亚群模式下请设置 subset_col 和 subset_idents！")
    cat(sprintf("  [INFO] 模式: 亚群分析 (筛选 %s == %s)\n", subset_col, paste(subset_idents, collapse=",")))
    Idents(work_obj) <- subset_col
    work_obj <- subset(work_obj, idents = subset_idents)
  } else {
    cat("  [INFO] 模式: 全局分析\n")
  }
  
  if (!group_by %in% colnames(work_obj@meta.data)) stop(paste("分组列不存在:", group_by))
  
  # --- 2. 构建 CDS 对象 ---
  cat("  [INFO] 构建 CellDataSet...\n")
  # 兼容 Seurat V4/V5 写法
  expr_matrix <- tryCatch({
    as(as.matrix(GetAssayData(work_obj, layer = "counts")), "sparseMatrix")
  }, error = function(e) {
    as(as.matrix(GetAssayData(work_obj, slot = "counts")), "sparseMatrix")
  })
  
  feature_ann <- data.frame(gene_id = rownames(expr_matrix), gene_short_name = rownames(expr_matrix))
  rownames(feature_ann) <- rownames(expr_matrix)
  pd <- new("AnnotatedDataFrame", data = work_obj@meta.data)
  fd <- new("AnnotatedDataFrame", data = feature_ann)
  
  cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  
  # --- 3. 预处理 ---
  cat("  [INFO] 估算 SizeFactors 和 Dispersions...\n")
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, cores = cores)
  cds <- detectGenes(cds, min_expr = 0.1)
  
  cat(sprintf("  [INFO] 筛选排序基因 (阈值: %s%% 的细胞表达)...\n", min_expr_prop * 100))
  expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= min_expr_prop * ncol(cds)))
  cds <- setOrderingFilter(cds, expressed_genes)
  
  # --- 4. 降维与轨迹 ---
  cat("  [INFO] DDRTree 降维 (这可能需要几分钟)...\n")
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
  
  # --- 5. 排序细胞 (应用 reverse 参数) ---
  # 根据您的要求，这里直接在 orderCells 中使用 reverse 参数
  if (is.null(reverse_traj)) {
    cat("  [INFO] 正在排序细胞 (Reverse = NULL, 自动默认)...\n")
    cds <- orderCells(cds)
  } else {
    cat(sprintf("  [INFO] 正在排序细胞 (Reverse = %s)...\n", reverse_traj))
    # 这里直接调用 reverse 参数
    cds <- orderCells(cds, reverse = reverse_traj)
  }
  
  # --- 6. 绘制核心轨迹图 (三联图) ---
  cat("  [INFO] 绘制轨迹图 (Cluster | Pseudotime | State)...\n")
  
  p1 <- plot_cell_trajectory(cds, color_by = group_by) + 
    theme(legend.position = "top") + 
    ggtitle(paste0("By ", group_by))
  
  p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime") + 
    theme(legend.position = "top") + 
    ggtitle("By Pseudotime")
  
  p3 <- plot_cell_trajectory(cds, color_by = "State") + 
    theme(legend.position = "top") + 
    ggtitle("By State")
  
  combined_plot <- p1 | p2 | p3
  
  ggsave(file.path(output_dir, paste0(file_prefix, "_Trajectory_Triple.png")), combined_plot, width = 18, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(file_prefix, "_Trajectory_Triple.pdf")), combined_plot, width = 18, height = 6)
  
  # -------------------------------------------------------------------------
  # 功能 A: 拟时序热图
  # -------------------------------------------------------------------------
  if (run_pseudotime_heatmap) {
    cat("  [INFO] 计算拟时序差异基因...\n")
    diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                          fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                          cores = cores)
    
    sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
    
    if (length(sig_gene_names) > 0) {
      out_csv_time <- file.path(output_dir, paste0(file_prefix, "_Pseudotime_All_Sig_Genes.csv"))
      write.csv(diff_test_res[sig_gene_names,], out_csv_time)
      
      top_genes <- head(sig_gene_names[order(diff_test_res[sig_gene_names,]$qval)], heatmap_top_n)
      
      out_ph_png <- file.path(output_dir, paste0(file_prefix, "_Heatmap_Pseudotime.png"))
      png(out_ph_png, width = 8, height = 10, units = "in", res = 300)
      plot_pseudotime_heatmap(cds[top_genes,], num_clusters = 4, cores = 1, show_rownames = TRUE)
      dev.off()
    }
  }
  
  # -------------------------------------------------------------------------
  # 功能 B: 分支热图 (BEAM)
  # -------------------------------------------------------------------------
  has_branches <- length(unique(pData(cds)$State)) > 1
  
  if (run_branch_heatmap && has_branches) {
    cat("  [INFO] 进行 BEAM 分支分析...\n")
    
    # 强制指定 progenitor_method 修复报错
    BEAM_res <- BEAM(cds[expressed_genes,], 
                     branch_point = 1, 
                     cores = cores,
                     progenitor_method = "duplicate") 
    
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
    sig_beam_genes <- row.names(subset(BEAM_res, qval < 0.01))
    
    if (length(sig_beam_genes) > 0) {
      write.csv(BEAM_res[sig_beam_genes,], file.path(output_dir, paste0(file_prefix, "_Branch_Dependent_Genes.csv")))
      
      top_beam_genes <- head(sig_beam_genes, heatmap_top_n)
      out_bh_png <- file.path(output_dir, paste0(file_prefix, "_Heatmap_Branched.png"))
      png(out_bh_png, width = 8, height = 10, units = "in", res = 300)
      plot_genes_branched_heatmap(cds[top_beam_genes,],
                                  branch_point = 1,
                                  num_clusters = 4,
                                  cores = 1,
                                  use_gene_short_name = TRUE,
                                  show_rownames = TRUE)
      dev.off()
    }
  } else if (run_branch_heatmap && !has_branches) {
    cat("  [WARN] 轨迹中未检测到分支，跳过 BEAM 分析。\n")
  }
  
  saveRDS(cds, file = file.path(output_dir, paste0(file_prefix, "_cds.rds")))
  cat("=================================================\n")
  return(cds)
}


# 3. 使用示例 ---------------------------------------------------------------

seu_obj <- readRDS("/home/lin/c_group/hep.rds")

# ==============================================
# 场景 A: 默认情况 (不反转)
# ==============================================
cds_a <- run_monocle2_analysis(
  seurat_obj = seu_obj,
  group_by = "cellcluster",
  
  reverse_traj = NULL,           # <--- 默认不反转
  
  subset_mode = FALSE,             
  output_dir = "/home/lin/c_group/monocle/result_A",
  file_prefix = "Normal_Run"
)
#levels(seu_obj@meta.data$cellcluster)
# ==============================================
# 场景 B: 反转起点 (Reverse = TRUE)
# 假如场景A跑出来的 Pseudotime 颜色是反的 (起点是深蓝，终点是浅蓝)
# 用这个参数可以把它们掉个头
# ==============================================
cds_b <- run_monocle2_analysis(
  seurat_obj = seu_obj,
  group_by = "seurat_clusters",
  
  reverse_traj = TRUE,           # <--- 设置为 TRUE 强制反转
  
  subset_mode = TRUE,              
  subset_col = "cellcluster",
  subset_idents = c("Core-Metabolic-Hep"),
  output_dir = "/home/lin/c_group/monocle/result_B_Reversed",
  file_prefix = "Reversed_Run"
)
