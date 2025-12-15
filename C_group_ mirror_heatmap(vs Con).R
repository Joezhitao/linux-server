#############################################################
# 单细胞差异分析 (Exp vs Con) 四格镜像热图与表格脚本 (v11.0 Table版)
# 功能：
# 1. 循环处理 6 个实验组 vs Con。
# 2. 生成 6 张四格镜像热图 (PDF)。
# 3. 生成 6 个对应的详细数据表格 (CSV)。
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

run_Exp_vs_Con_heatmap_and_table <- function(
    seurat_file = "/home/lin/c_group/hep.rds",
    exp_groups = c("Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R"), 
    control_group = "Con",
    logfc_threshold = 0.25,
    p_val_cutoff = 0.05,
    top_n_pathways = 20,  
    output_dir = "/home/lin/c_group/GO_Vs_Con_Mirror_Results" 
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 读取数据
  if (is.character(seurat_file)) {
    obj <- readRDS(seurat_file)
  } else {
    obj <- seurat_file
  }
  Idents(obj) <- obj@meta.data$orig.ident
  
  if (!control_group %in% levels(Idents(obj))) {
    stop(paste("Control group", control_group, "not found!"))
  }
  
  # --- 循环处理每一个实验组 ---
  for (exp_grp in exp_groups) {
    cat(sprintf("\n================ Processing: %s vs %s ================\n", exp_grp, control_group))
    
    if (!exp_grp %in% levels(Idents(obj))) next
    
    # 1. 差异分析
    cat("  -> Diff Analysis...\n")
    markers <- FindMarkers(object = obj, ident.1 = exp_grp, ident.2 = control_group, 
                           group.by = "orig.ident", test.use = "wilcox", 
                           min.pct = 0.1, logfc.threshold = logfc_threshold, 
                           only.pos = FALSE, verbose = FALSE)
    
    if (nrow(markers) == 0) next
    markers$gene <- rownames(markers)
    
    genes_Exp_High <- markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC > 0) %>% pull(gene)
    genes_Con_High <- markers %>% filter(p_val_adj < p_val_cutoff & avg_log2FC < 0) %>% pull(gene)
    
    # 2. GO 富集 & 数据构建
    
    # --- Part 1: Exp High (High in Exp) ---
    go_Exp <- run_go_enrichment(genes_Exp_High)
    df_Exp_pathways <- data.frame()
    
    if (!is.null(go_Exp)) {
      top_Exp <- go_Exp %>% arrange(p.adjust) %>% slice_head(n = top_n_pathways)
      # 为了生成表格，我们要保留更多原始信息
      for (i in 1:nrow(top_Exp)) {
        base_score <- -log10(top_Exp$p.adjust[i])
        desc <- top_Exp$Description[i]
        orig_p <- top_Exp$p.adjust[i]
        gene_ratio <- top_Exp$GeneRatio[i]
        
        # Exp列 (红)
        df_Exp_pathways <- rbind(df_Exp_pathways, data.frame(
          Description = desc, Group = exp_grp, Signed_Score = base_score, 
          Raw_Log10P = base_score, P_adj = orig_p, GeneRatio = gene_ratio, Type = paste0("High in ", exp_grp)
        ))
        # Con列 (蓝)
        df_Exp_pathways <- rbind(df_Exp_pathways, data.frame(
          Description = desc, Group = control_group, Signed_Score = -base_score, 
          Raw_Log10P = base_score, P_adj = orig_p, GeneRatio = gene_ratio, Type = paste0("High in ", exp_grp)
        ))
      }
    }
    
    # --- Part 2: Con High (High in Con) ---
    go_Con <- run_go_enrichment(genes_Con_High)
    df_Con_pathways <- data.frame()
    
    if (!is.null(go_Con)) {
      top_Con <- go_Con %>% arrange(p.adjust) %>% slice_head(n = top_n_pathways)
      for (i in 1:nrow(top_Con)) {
        base_score <- -log10(top_Con$p.adjust[i])
        desc <- top_Con$Description[i]
        orig_p <- top_Con$p.adjust[i]
        gene_ratio <- top_Con$GeneRatio[i]
        
        # Exp列 (蓝)
        df_Con_pathways <- rbind(df_Con_pathways, data.frame(
          Description = desc, Group = exp_grp, Signed_Score = -base_score, 
          Raw_Log10P = base_score, P_adj = orig_p, GeneRatio = gene_ratio, Type = paste0("High in ", control_group)
        ))
        # Con列 (红)
        df_Con_pathways <- rbind(df_Con_pathways, data.frame(
          Description = desc, Group = control_group, Signed_Score = base_score, 
          Raw_Log10P = base_score, P_adj = orig_p, GeneRatio = gene_ratio, Type = paste0("High in ", control_group)
        ))
      }
    }
    
    # 3. 合并
    final_df <- rbind(df_Exp_pathways, df_Con_pathways)
    if (nrow(final_df) == 0) next
    
    # --- 保存表格 (新增部分) ---
    table_filename <- paste0("GO_Table_", exp_grp, "_vs_", control_group, ".csv")
    write.csv(final_df, file.path(output_dir, table_filename), row.names = FALSE)
    cat(sprintf("  [DONE] Table Saved: %s\n", table_filename))
    
    # 4. 绘图准备
    final_df$Group <- factor(final_df$Group, levels = c(exp_grp, control_group))
    panel_levels <- c(paste0("High in ", control_group), paste0("High in ", exp_grp))
    final_df$Type <- factor(final_df$Type, levels = panel_levels)
    
    final_df <- final_df %>%
      group_by(Type) %>%
      mutate(AbsScore = abs(Signed_Score)) %>%
      arrange(Type, AbsScore) %>%
      ungroup()
    
    final_df$Description <- factor(final_df$Description, levels = unique(final_df$Description))
    
    # 5. 绘图
    p <- ggplot(final_df, aes(x = Group, y = Description, fill = Signed_Score)) +
      geom_tile(color = "white", size = 0.2) +
      scale_fill_gradient2(
        low = "#3288BD", mid = "white", high = "#D53E4F", midpoint = 0,
        name = expression(Signed ~ -log[10](P.adj))
      ) +
      facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
      labs(
        title = paste0(exp_grp, " vs ", control_group),
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    
    # 6. 保存图片
    plot_height <- max(5, length(unique(final_df$Description)) * 0.3 + 2)
    plot_filename <- paste0("Mirror_Heatmap_", exp_grp, "_vs_", control_group, ".pdf")
    ggsave(file.path(output_dir, plot_filename), p, width = 10, height = plot_height)
    cat(sprintf("  [DONE] Plot Saved: %s\n", plot_filename))
  }
  
  cat("\n所有 6 组表格和图片均已生成完毕！\n")
}

# 4. 运行 -------------------------------------------------------------------

seurat_obj <- readRDS("/home/lin/c_group/hep.rds")

run_Exp_vs_Con_heatmap_and_table(
  seurat_file = seurat_obj,
  exp_groups = c("Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R"), 
  control_group = "Con",
  output_dir = "/home/lin/c_group/New_Folder/Vs_Con_Results"
)
