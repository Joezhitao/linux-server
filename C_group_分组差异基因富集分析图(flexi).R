#############################################################
# 单细胞成对差异分析 (L vs R) 与 GO 富集绘图脚本 (v6.1 自定义尺寸版)
# 功能：L vs R 差异 -> GO富集 -> 气泡图
# 特点：
# 1. 布局：单图分面 (facet)，共用左侧Y轴。
# 2. 等宽：强制三列等宽。
# 3. 排序：以【Cd1】组的Count从大到小为基准锁定Y轴顺序。
# 4. 内容：只展示左叶上调 (Up_in_Left)。
# 5. 输出：支持自定义 PDF 长宽和 DPI。
# 作者：Sci专家
#############################################################

rm(list = ls())
gc()

# 1. 加载必要的包 -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(openxlsx)
  library(forcats)
  library(clusterProfiler)
  
  if(!requireNamespace("org.Rn.eg.db", quietly=TRUE)) {
    library(org.Rn.eg.db)
  } else {
    library(org.Rn.eg.db)
  }
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 核心功能函数 -----------------------------------------------------------

calculate_enrich_factor <- function(df) {
  df$GeneRatio_Val <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  df$BgRatio_Val <- sapply(strsplit(df$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  df$EnrichFactor <- df$GeneRatio_Val / df$BgRatio_Val
  return(df)
}

# 3. 主流程 -----------------------------------------------------------------

run_L_vs_R_go_analysis <- function(
    seurat_file = "/home/lin/c_group/hep.rds",
    timepoints = c("Cd1", "Cd15", "Cd30"), 
    species = "rat",
    logfc_threshold = 0.25,
    p_val_cutoff = 0.05,
    top_n_pathways = 30,
    output_dir = "/home/lin/c_group/GO_LR_Comparison",
    # --- 新增绘图参数 ---
    pdf_width = 12,       # PDF 宽度 (英寸)
    pdf_height = NULL,    # PDF 高度 (英寸)，设为 NULL 则根据通路数量自动计算
    plot_dpi = 300        # 图片分辨率
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始分析...展示 Top", top_n_pathways, "条通路\n")
  
  # 读取数据
  if (is.character(seurat_file)) {
    obj <- readRDS(seurat_file)
  } else {
    obj <- seurat_file
  }
  
  # --- 1. 差异分析 ---
  cat("  [INFO] 执行成对差异分析...\n")
  all_markers <- data.frame()
  Idents(obj) <- obj@meta.data$orig.ident
  
  for (tp in timepoints) {
    left_ident <- paste0(tp, "L")
    right_ident <- paste0(tp, "R")
    
    available_idents <- levels(Idents(obj))
    if (!left_ident %in% available_idents || !right_ident %in% available_idents) next
    
    cat(sprintf("    -> 比较: %s vs %s\n", left_ident, right_ident))
    markers <- FindMarkers(object = obj, ident.1 = left_ident, ident.2 = right_ident, 
                           group.by = "orig.ident", test.use = "wilcox", 
                           min.pct = 0.1, logfc.threshold = logfc_threshold, 
                           only.pos = FALSE, verbose = FALSE)
    
    if (nrow(markers) > 0) {
      markers$gene <- rownames(markers)
      markers$cluster <- tp
      markers$direction <- ifelse(markers$avg_log2FC > 0, "Up_in_Left", "Up_in_Right")
      all_markers <- rbind(all_markers, markers)
    }
  }
  
  sig_markers <- all_markers %>% filter(p_val_adj < p_val_cutoff)
  write.xlsx(all_markers, file.path(output_dir, "All_L_vs_R_DiffGenes.xlsx"), rowNames = F)
  
  # --- 2. GO 富集 (只看 Up_in_Left) ---
  target_direction <- "Up_in_Left"
  cat(sprintf("  [INFO] GO 富集 (方向: %s)...\n", target_direction))
  
  org_db <- org.Rn.eg.db
  all_go_results <- list()
  groups_to_plot <- intersect(timepoints, unique(sig_markers$cluster))
  
  for (grp in groups_to_plot) {
    genes <- sig_markers %>% filter(cluster == grp & direction == target_direction) %>% pull(gene)
    if(length(genes) < 5) next
    
    gene_conv <- tryCatch({ bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db) }, error = function(e) NULL)
    if (is.null(gene_conv)) next
    
    ego <- enrichGO(gene = gene_conv$ENTREZID, OrgDb = org_db, ont = "BP", 
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
    
    if (!is.null(ego) && nrow(ego) > 0) {
      res <- as.data.frame(ego)
      res$Group <- grp
      all_go_results[[grp]] <- res
    }
  }
  
  if (length(all_go_results) == 0) stop("未发现显著结果。")
  final_df <- do.call(rbind, all_go_results)
  final_df <- calculate_enrich_factor(final_df)
  write.csv(final_df, file.path(output_dir, paste0("GO_Enrichment_", target_direction, ".csv")), row.names = F)
  
  # --- 3. 绘图数据准备 ---
  cat("  [INFO] 整理绘图数据...\n")
  
  # 筛选数据：每个组只留 Top N
  plot_data <- final_df %>%
    group_by(Group) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n_pathways) %>%
    ungroup()
  
  # 强制 Group 顺序：Cd1 -> Cd15 -> Cd30
  plot_data$Group <- factor(plot_data$Group, levels = timepoints)
  
  # 处理 p值
  plot_data$neg_log10_p <- -log10(plot_data$p.adjust)
  max_val <- max(plot_data$neg_log10_p[is.finite(plot_data$neg_log10_p)], na.rm=T)
  plot_data$neg_log10_p[is.infinite(plot_data$neg_log10_p)] <- max_val * 1.1
  
  # === 【关键排序逻辑：以 Cd1 为基准】 ===
  # 1. 找出 Cd1 组的数据
  cd1_data <- plot_data %>% filter(Group == timepoints[1])
  
  # 2. 找出 Cd1 中 Count 的顺序 (Count 从小到大排，这样画图时大的在上面)
  if(nrow(cd1_data) > 0) {
    # 如果有重复的Description，取Count最大的那个
    cd1_order <- cd1_data %>% 
      arrange(Count) %>% 
      pull(Description) %>% 
      unique()
    
    # 3. 获取其他组特有的 Description (不在Cd1里的)
    other_desc <- setdiff(unique(plot_data$Description), cd1_order)
    
    # 4. 合并顺序：Cd1的排前面，其他的排后面
    final_levels <- c(other_desc, cd1_order)
    
    # 5. 应用因子顺序
    plot_data$Description <- factor(plot_data$Description, levels = final_levels)
  } else {
    # 如果没Cd1数据，就全局按Count排
    plot_data$Description <- fct_reorder(plot_data$Description, plot_data$Count)
  }
  
  # --- 4. 绘图 ---
  cat("  [INFO] 正在绘图...\n")
  
  p <- ggplot(plot_data, aes(x = EnrichFactor, y = Description)) +
    geom_point(aes(size = Count, color = neg_log10_p)) +
    
    # 配色：深蓝 -> 红
    scale_color_gradient(low = "darkblue", high = "red") +
    
    # 气泡大小：整体调大
    scale_size_continuous(range = c(4, 10)) +
    
    labs(
      color = expression(-log[10](p.adjust)),
      size = "Gene Number",
      y = NULL,
      x = "Enrichment Factor",
      title = paste("L vs R Comparison:", target_direction)
    ) +
    
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      strip.text = element_text(size = 14, face = "bold"), 
      strip.background = element_rect(fill = "#E5E5E5"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    
    # 文字自动换行
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    
    # 【核心：等宽设置】
    facet_grid(~ Group, scales = "free_x", space = "fixed")
  
  # --- 5. 保存设置 ---
  filename <- paste0("Dotplot_L_vs_R_", target_direction, "_Cd1_Sorted.pdf")
  
  # 决定高度：如果用户没给(NULL)，就自动算；如果给了，就用给定的
  if (is.null(pdf_height)) {
    # 自动计算逻辑：基础高度 + 每个通路给一点空间
    final_height <- max(6, length(unique(plot_data$Description)) * 0.4) 
    cat(sprintf("  [Auto-Height] 自动计算高度为: %.2f inch\n", final_height))
  } else {
    final_height <- pdf_height
    cat(sprintf("  [Manual-Height] 使用自定义高度: %.2f inch\n", final_height))
  }
  
  ggsave(file.path(output_dir, filename), p, 
         width = pdf_width, 
         height = final_height, 
         dpi = plot_dpi)
  
  cat(sprintf("  [完成] 结果保存在: %s\n", output_dir))
  return(p)
}

# 4. 运行示例 ---------------------------------------------------------------

# 读取对象 (请替换为实际路径)
seurat_obj <- readRDS("/home/lin/c_group/hep.rds")

# 运行 (这里设置了长宽和DPI)
p_final <- run_L_vs_R_go_analysis(
  seurat_file = seurat_obj,
  timepoints = c("Cd1", "Cd15", "Cd30"), 
  top_n_pathways = 20,  
  output_dir = "/home/lin/c_group/New Folder/GO_Results",
  
  # --- 设置 PDF 尺寸 ---
  pdf_width = 14,    # 宽度设为 14 英寸
  pdf_height = 10,   # 高度设为 10 英寸 (若不想指定，写 NULL)
  plot_dpi = 600     # DPI 设为 600
)

print(p_final)
