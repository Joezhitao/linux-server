#############################################################
# 单细胞多分组 GO 富集分析与可视化脚本 (v3.1 宽度修正版)
# 功能：自动读取分组顺序 -> 找差异基因 -> GO富集 -> 仿Nature气泡图
# 特点：修复了某些组因富集因子过大导致其他组被挤压的问题
# 作者：Sci专家
# 日期：2023-10-27
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
  library(org.Rn.eg.db)
  library(clusterProfiler)
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 核心工具函数 -----------------------------------------------------------

# 自动匹配物种数据库
get_species_db <- function(species) {
  species_map <- list(
    "human" = "org.Hs.eg.db", "hs" = "org.Hs.eg.db",
    "mouse" = "org.Mm.eg.db", "mm" = "org.Mm.eg.db",
    "rat"   = "org.Rn.eg.db", "rn" = "org.Rn.eg.db"
  )
  
  pkg <- species_map[[tolower(species)]]
  if (is.null(pkg)) stop("仅支持 human, mouse, rat")
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("请安装物种包:", pkg))
  }
  library(pkg, character.only = TRUE)
  return(get(pkg))
}

# 计算富集因子 (Enrichment Factor)
calculate_enrich_factor <- function(df) {
  # GeneRatio "10/100" -> 0.1
  df$GeneRatio_Val <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  # BgRatio "500/10000" -> 0.05
  df$BgRatio_Val <- sapply(strsplit(df$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  # Factor = 0.1 / 0.05 = 2
  df$EnrichFactor <- df$GeneRatio_Val / df$BgRatio_Val
  return(df)
}

# 3. 主分析函数 -------------------------------------------------------------

run_auto_go_analysis <- function(
    seurat_data,                 # Seurat 对象或路径
    species = "rat",             # 物种
    group_by = "orig.ident",     # 分组列
    subset_mode = FALSE,         # 是否亚群分析
    subset_col = NULL,           
    subset_idents = NULL,        
    logfc_threshold = 0.25,      # 差异基因阈值
    top_n_pathways = 5,          # 展示通路数
    output_dir = "./GO_result", 
    file_prefix = "GO_Analysis"
) {
  
  # --- 0. 准备工作 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat("开始自动化 GO 富集分析流程\n")
  
  # 读取数据
  if (is.character(seurat_data)) {
    cat(sprintf("  [INFO] 读取 Seurat 文件: %s \n", seurat_data))
    obj <- readRDS(seurat_data)
  } else {
    obj <- seurat_data
  }
  
  # --- 1. 数据筛选与分组处理 ---
  final_obj <- obj
  
  # 亚群过滤模式
  if (subset_mode) {
    cat(sprintf("  [INFO] 筛选亚群: %s == %s\n", subset_col, paste(subset_idents, collapse=",")))
    Idents(final_obj) <- subset_col
    final_obj <- subset(final_obj, idents = subset_idents)
  }
  
  # 检查分组列
  if (!group_by %in% colnames(final_obj@meta.data)) stop(paste("元数据中找不到:", group_by))
  
  # *** 核心修改：自动获取分组顺序 ***
  Idents(final_obj) <- final_obj@meta.data[[group_by]]
  current_levels <- levels(Idents(final_obj))
  if (is.null(current_levels)) {
    current_levels <- sort(unique(Idents(final_obj))) 
  }
  cat(sprintf("  [INFO] 检测到分组顺序: %s\n", paste(current_levels, collapse = " -> ")))
  
  # --- 2. 差异分析 (FindAllMarkers) ---
  cat("  [INFO] 正在计算各分组的高表达基因 (FindAllMarkers)...\n")
  markers <- FindAllMarkers(final_obj, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = logfc_threshold, verbose = FALSE)
  
  sig_markers <- markers %>% filter(p_val_adj < 0.05)
  cat(sprintf("  [INFO] 找到 %d 个显著基因\n", nrow(sig_markers)))
  
  # --- 3. GO 富集 (循环每个组) ---
  org_db <- get_species_db(species)
  all_go_results <- list()
  
  # 只分析数据中存在的组
  groups_to_test <- current_levels[current_levels %in% unique(sig_markers$cluster)]
  
  cat("  [INFO] 开始 GO (BP) 富集...\n")
  for (grp in groups_to_test) {
    genes <- sig_markers %>% filter(cluster == grp) %>% pull(gene)
    
    # 转换基因 ID
    gene_conv <- tryCatch({
      bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
    }, error = function(e) NULL)
    
    if (is.null(gene_conv) || nrow(gene_conv) < 3) next
    
    # 运行 GO
    ego <- enrichGO(gene = gene_conv$ENTREZID,
                    OrgDb = org_db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
    
    if (!is.null(ego) && nrow(ego) > 0) {
      res <- as.data.frame(ego)
      res$Group <- grp
      all_go_results[[grp]] <- res
    }
  }
  
  if (length(all_go_results) == 0) stop("未发现显著的富集通路。")
  final_df <- do.call(rbind, all_go_results)
  
  # --- 4. 绘图数据准备 ---
  cat("  [INFO] 整理绘图数据...\n")
  final_df <- calculate_enrich_factor(final_df)
  
  # 筛选 Top N
  plot_data <- final_df %>%
    group_by(Group) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n_pathways) %>%
    ungroup()
  
  # *** 关键：应用分组顺序 ***
  plot_data$Group <- factor(plot_data$Group, levels = current_levels)
  
  # 处理 p值显示 (避免 Inf)
  plot_data$neg_log10_p <- -log10(plot_data$p.adjust)
  max_val <- max(plot_data$neg_log10_p[is.finite(plot_data$neg_log10_p)], na.rm=T)
  plot_data$neg_log10_p[is.infinite(plot_data$neg_log10_p)] <- max_val * 1.1
  
  # 导出表格
  write.csv(final_df, file.path(output_dir, paste0(file_prefix, "_All.csv")), row.names = F)
  
  # --- 5. 绘制气泡图 ---
  cat("  [INFO] 生成分面气泡图 (已修正宽度)...\n")
  
  p <- ggplot(plot_data, aes(x = EnrichFactor, y = Description)) +
    geom_point(aes(size = Count, color = neg_log10_p)) +
    scale_color_gradient(low = "green", high = "red") +
    
    labs(
      color = expression(-log[10](p.adjust)), # 数学表达式图例
      size = "Gene Number",
      y = NULL,
      x = "Enrichment Factor"
    ) +
    
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"), # 字体稍微调小一点，防重叠
      strip.text = element_text(size = 14, face = "bold"), 
      strip.background = element_rect(fill = "#E5E5E5"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines") # 增加一点分面间距
    ) +
    
    # 自动换行
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    
    # *** 关键修改：去掉了 space="free_x"，强制等宽 ***
    facet_grid(~ Group, scales = "free_x")
  
  # 保存
  # 自动计算宽度，每多一个组增加一些宽度
  plot_width <- 4 + 1.8 * length(groups_to_test)
  ggsave(file.path(output_dir, paste0(file_prefix, "_Dotplot.pdf")), p, 
         width = plot_width, height = 8)
  
  cat(sprintf("  [完成] 结果已保存至: %s\n", output_dir))
  return(p)
}

# 4. 运行部分 ---------------------------------------------------------------

seurat_obj <- readRDS("/home/lin/c_group/hep.rds")

# ==============================================
# 场景 A: 全局分析 (使用 orig.ident 的默认顺序)
# ==============================================
p1 <- run_auto_go_analysis(
  seurat_data = seurat_obj,
  species = "rat",                
  group_by = "orig.ident",        
  top_n_pathways = 5,             
  output_dir = "/home/lin/c_group/New Folder/GO_Results",
  file_prefix = "Global_Hep_Analysis"
)

# ==============================================
# 场景 B: 亚群分析
# ==============================================
# p2 <- run_auto_go_analysis(
#   seurat_data = seurat_obj,
#   species = "rat",
#   subset_mode = TRUE,
#   subset_col = "cell_type",      
#   subset_idents = "Hepatocytes", 
#   group_by = "orig.ident",       
#   top_n_pathways = 5,
#   output_dir = "/home/lin/c_group/GO_Results",
#   file_prefix = "Subgroup_Hep_Analysis"
# )

print(p1)
