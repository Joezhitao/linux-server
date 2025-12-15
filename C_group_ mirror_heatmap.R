#############################################################
# 单细胞成对差异分析 (L vs R) 四格镜像热图 (v9.0 Mirror版)
# 功能：
# 1. 针对每天 (Cd1, Cd15...) 生成一张独立的图。
# 2. 布局为“四格图”效果：
#    - 左列 (Cd1L): 上蓝下红
#    - 右列 (Cd1R): 上红下蓝
# 3. 完美解决宽度不一致问题，两列等宽。
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

run_go_enrichment <- function(genes, org_db = org.Rn.eg.db) {
  if(length(genes) < 3) return(NULL) 
  
  gene_conv <- tryCatch({ 
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db) 
  }, error = function(e) NULL)
  
  if (is.null(gene_conv)) return(NULL)
  
  ego <- enrichGO(gene = gene_conv$ENTREZID, OrgDb = org_db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    return(as.data.frame(ego))
  } else {
    return(NULL)
  }
}

# 3. 主流程函数 -------------------------------------------------------------

run_mirror_heatmap_analysis <- function(
    seurat_file = "/home/lin/c_group/hep.rds",
    timepoints = c("Cd1", "Cd15", "Cd30"), 
    logfc_threshold = 0.25,
    p_val_cutoff = 0.05,
    top_n_pathways = 15,  # 上下各展示多少条通路
    output_dir = "/home/lin/c_group/GO_Mirror_Heatmap" 
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 读取数据
  if (is.character(seurat_file)) {
    obj <- readRDS(seurat_file)
  } else {
    obj <- seurat_file
  }
  
  Idents(obj) <- obj@meta.data$orig.ident
  
  # --- 循环处理每一个时间点 ---
  for (tp in timepoints) {
    cat(sprintf("\n================ Processing: %s ================\n", tp))
    
    left_ident <- paste0(tp, "L")
    right_ident <- paste0(tp, "R")
    
    available_idents <- levels(Idents(obj))
    if (!left_ident %in% available_idents || !right_ident %in% available_idents) next
    
    # 1. 差异分析 (Cd1L vs Cd1R)
    cat(sprintf("  -> Diff Analysis: %s vs %s\n", left_ident, right_ident))
    markers <- FindMarkers(object = obj, ident.1 = left_ident, ident.2 = right_ident, 
                           group.by = "orig.ident", test.use = "wilcox", 
                           min.pct = 0.1, logfc.threshold = logfc_threshold, 
                           only.pos = FALSE, verbose = FALSE)
    
    if (nrow(markers) == 0) next
    
    markers$gene <- rownames(markers)
    
    # 基因集 A: 左高右低 (Up in L)
    genes_L_High <- markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC > 0) %>% pull(gene)
    # 基因集 B: 右高左低 (Up in R)
    genes_R_High <- markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC < 0) %>% pull(gene)
    
    cat(sprintf("     Genes L_High: %d | Genes R_High: %d\n", length(genes_L_High), length(genes_R_High)))
    
    # 2. 运行 GO 富集
    # 我们只运行“高表达”的那一边，然后通过镜像逻辑构建另一边的数据
    
    # --- 通路集 1: 左侧特异性通路 (L红, R蓝) ---
    go_L_High <- run_go_enrichment(genes_L_High)
    df_L_pathways <- data.frame()
    
    if (!is.null(go_L_High)) {
      top_L <- go_L_High %>% arrange(p.adjust) %>% slice_head(n = top_n_pathways)
      
      # 构建镜像数据：这些通路在 L 中是 Up (Score > 0)，在 R 中是 Down (Score < 0)
      for (i in 1:nrow(top_L)) {
        score <- -log10(top_L$p.adjust[i])
        desc <- top_L$Description[i]
        
        # 添加 Cd1L 数据 (红色)
        df_L_pathways <- rbind(df_L_pathways, data.frame(
          Description = desc, Group = left_ident, Score = score, Type = "High in Left (L>R)"
        ))
        # 添加 Cd1R 数据 (蓝色)
        df_L_pathways <- rbind(df_L_pathways, data.frame(
          Description = desc, Group = right_ident, Score = -score, Type = "High in Left (L>R)"
        ))
      }
    }
    
    # --- 通路集 2: 右侧特异性通路 (R红, L蓝) ---
    go_R_High <- run_go_enrichment(genes_R_High)
    df_R_pathways <- data.frame()
    
    if (!is.null(go_R_High)) {
      top_R <- go_R_High %>% arrange(p.adjust) %>% slice_head(n = top_n_pathways)
      
      # 构建镜像数据：这些通路在 R 中是 Up (Score > 0)，在 L 中是 Down (Score < 0)
      for (i in 1:nrow(top_R)) {
        score <- -log10(top_R$p.adjust[i])
        desc <- top_R$Description[i]
        
        # 添加 Cd1L 数据 (蓝色 - 因为这是右边高的通路)
        df_R_pathways <- rbind(df_R_pathways, data.frame(
          Description = desc, Group = left_ident, Score = -score, Type = "High in Right (R>L)"
        ))
        # 添加 Cd1R 数据 (红色)
        df_R_pathways <- rbind(df_R_pathways, data.frame(
          Description = desc, Group = right_ident, Score = score, Type = "High in Right (R>L)"
        ))
      }
    }
    
    # 3. 合并数据
    final_df <- rbind(df_L_pathways, df_R_pathways)
    if (nrow(final_df) == 0) next
    
    # 4. 关键：设置因子顺序以控制四格图布局
    
    # X轴顺序：Cd1L 在左，Cd1R 在右
    final_df$Group <- factor(final_df$Group, levels = c(left_ident, right_ident))
    
    # 面板顺序：
    # 按照您的要求：左边上面是蓝的 -> 意味着上面那个面板必须是 "High in Right" 的通路
    # 因为 "High in Right" 的通路在左侧(Cd1L) 是下调(蓝)的。
    # 所以我们将 "High in Right" 设为第一个 level
    final_df$Type <- factor(final_df$Type, levels = c("High in Right (R>L)", "High in Left (L>R)"))
    
    # Y轴排序：在每个 Type 内部，按 Score 的绝对值排序，让颜色深的在一起
    final_df <- final_df %>%
      group_by(Type) %>%
      mutate(AbsScore = abs(Score)) %>%
      arrange(Type, AbsScore) %>%
      ungroup()
    
    final_df$Description <- factor(final_df$Description, levels = unique(final_df$Description))
    
    # 5. 绘图 (Heatmap)
    cat("  -> Plotting Mirror Heatmap...\n")
    
    p <- ggplot(final_df, aes(x = Group, y = Description, fill = Score)) +
      # 画热图方块，color="white" 加个白边框更好看
      geom_tile(color = "white", size = 0.2) +
      
      # 颜色：红蓝对抗
      # 负值(蓝)=Down, 0(白), 正值(红)=Up
      scale_fill_gradient2(
        low = "#3288BD", 
        mid = "white", 
        high = "#D53E4F", 
        midpoint = 0,
        name = expression(Signed ~ -log[10](P.adj))
      ) +
      
      # 分面：上下两部分
      # scales="free_y" 保证上面的格子只显示上面的通路，下面显示下面的
      facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
      
      labs(
        title = paste0("Pathway Dysregulation: ", tp, " (Left vs Right)"),
        subtitle = "Red = Up-regulated | Blue = Down-regulated",
        x = NULL, y = NULL
      ) +
      
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 14, face = "bold", color = "black"), # X轴(Cd1L/Cd1R)字体大一点
        axis.text.y = element_text(size = 11, color = "black"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), # 热图不要网格线
        strip.background = element_rect(fill = "grey90"), # 分面标题背景
        strip.text = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      
      # 自动换行防止通路名太长
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    
    # 6. 保存
    # 高度自适应：通路越多图越高
    n_pathways <- length(unique(final_df$Description))
    plot_height <- max(5, n_pathways * 0.3 + 2)
    
    out_file <- file.path(output_dir, paste0("Mirror_Heatmap_", tp, ".pdf"))
    ggsave(out_file, p, width = 10, height = plot_height)
    
    cat(sprintf("  [DONE] Saved: %s\n", out_file))
  }
  
  cat("\n所有分析完成！\n")
}

# 4. 运行 -------------------------------------------------------------------

seurat_obj <- readRDS("/home/lin/c_group/hep.rds")

run_mirror_heatmap_analysis(
  seurat_file = seurat_obj,
  timepoints = c("Cd1", "Cd15", "Cd30"), 
  top_n_pathways = 20, # 上下各20条，总共40条
  output_dir = "/home/lin/c_group/New_Folder/Mirror_Heatmap_Results"
)
