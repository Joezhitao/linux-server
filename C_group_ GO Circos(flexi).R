#############################################################
# 单细胞 GO Circos 圈图 (v6.4 浅蓝配色版)
# 功能：Seurat 差异 -> GOChord -> 深度定制图例
# 特点：
# 1. 颜色优化：logFC 渐变从【浅蓝 (cornflowerblue)】到【红色】。
# 2. 去除中间透明色，适合只有正值 logFC 的情况。
# 3. 保持强制大写逻辑，杜绝报错。
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
  library(GOplot)
  library(clusterProfiler)
  
  if(!requireNamespace("org.Rn.eg.db", quietly=TRUE)) library(org.Rn.eg.db)
})

options(future.globals.maxSize = 100 * 1024^3)

# 2. 主分析流程 -------------------------------------------------------------

run_go_circos <- function(
    seurat_data, 
    species = "rat",
    scene = "B",               
    cluster_col = "seurat_clusters", 
    target_cluster = NULL,           
    ident_1 = NULL,            
    ident_2 = NULL,            
    group_by_col = "orig.ident", 
    top_n_pathways = 5,        
    logfc_threshold = 0.25,
    p_val_cutoff = 0.05,
    
    # --- 尺寸与字体控制 ---
    pdf_width = 13,       
    pdf_height = 10,      
    dpi = 300,            
    gene_font_size = 4,   
    
    output_dir = "./GO_Circos_Result"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat(sprintf("开始 GO Circos 分析 | 场景: %s\n", scene))
  
  # --- 1. 数据准备 ---
  if (is.character(seurat_data)) obj <- readRDS(seurat_data) else obj <- seurat_data
  
  markers <- NULL
  title_prefix <- ""
  
  if (scene == "A") {
    Idents(obj) <- obj@meta.data[[cluster_col]]
    markers <- FindMarkers(obj, ident.1 = target_cluster, only.pos = TRUE, 
                           min.pct = 0.25, logfc.threshold = logfc_threshold, verbose = FALSE)
    title_prefix <- paste0(target_cluster, "_Cluster")
  } else {
    Idents(obj) <- obj@meta.data[[group_by_col]]
    markers <- FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2, only.pos = TRUE,
                           min.pct = 0.1, logfc.threshold = logfc_threshold, verbose = FALSE)
    title_prefix <- paste0(ident_1, "_vs_", ident_2)
  }
  
  sig_markers <- markers %>% filter(p_val_adj < p_val_cutoff)
  sig_markers$Symbol <- rownames(sig_markers)
  if (nrow(sig_markers) < 5) stop("显著基因太少，无法画图")
  
  # --- 2. ID 转换 ---
  cat("  [INFO] ID 转换...\n")
  org_pkg <- switch(species, "rat" = "org.Rn.eg.db", "mouse" = "org.Mm.eg.db", "human" = "org.Hs.eg.db")
  library(org_pkg, character.only = TRUE)
  gene_map <- bitr(sig_markers$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = get(org_pkg))
  
  if (is.null(gene_map)) stop("ID 转换失败")
  sig_markers <- merge(sig_markers, gene_map, by.x = "Symbol", by.y = "SYMBOL")
  
  # --- 3. GO 富集 ---
  cat("  [INFO] 执行 GO 分析...\n")
  ego <- enrichGO(gene = sig_markers$ENTREZID, OrgDb = get(org_pkg), ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(ego) == 0) stop("未富集到 GO 通路")
  
  # --- 4. 数据整理 ---
  cat("  [INFO] 构建绘图数据...\n")
  go_full_df <- as.data.frame(ego)
  go_top <- go_full_df %>% arrange(p.adjust) %>% head(top_n_pathways)
  
  genelist <- data.frame(ID = toupper(sig_markers$Symbol), logFC = sig_markers$avg_log2FC)
  genelist <- genelist %>% group_by(ID) %>% summarize(logFC = mean(logFC)) %>% as.data.frame()
  
  clean_genes <- function(s) { paste(toupper(unlist(strsplit(s, "/"))), collapse = ", ") }
  
  david <- data.frame(
    Category = "BP",
    ID = go_top$ID,
    Term = go_top$Description,
    Genes = sapply(go_top$geneID, clean_genes),
    adj_pval = go_top$p.adjust
  )
  
  circ <- circle_dat(david, genelist)
  chord <- chord_dat(data = circ, genes = genelist, process = david$Term)
  
  # --- 5. 绘图 ---
  cat("  [INFO] 正在绘图...\n")
  
  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  use_cols <- all_cols[1:length(david$Term)]
  names(use_cols) <- david$Term 
  
  # 5.1 基础图
  p <- GOChord(chord, 
               space = 0.02,           
               gene.order = 'logFC',   
               gene.space = 0.25,      
               gene.size = gene_font_size,  
               ribbon.col = use_cols   
  )
  
  # 5.2 提取 Limits
  gb <- ggplot_build(p)
  idx <- which(vapply(gb$plot$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1)))
  sc <- gb$plot$scales$scales[[idx]]
  lims <- tryCatch({ sc$limits %||% range(sc$range$range) }, error = function(e) c(-1, 1))
  brks <- pretty(lims, n = 5)  
  
  # 5.3 深度定制图例
  p_final <- p + 
    guides(shape = "none", 
           size  = guide_legend(
             title = "GO Terms", order = 1,
             ncol  = 1, byrow = TRUE,
             override.aes = list(shape = 22, fill = use_cols, size = 8), 
             theme = theme(legend.text = element_text(color="black", size=10))
           )
    ) +
    # === 【颜色修改点】 ===
    scale_fill_gradient(
      name   = "logFC",
      limits = lims,
      breaks = brks,
      labels = brks,
      
      # 这里把原来的 "blue" 改为了 "cornflowerblue" (矢车菊蓝/浅蓝)
      # 这种蓝色比较柔和，不会黑乎乎的
      low = "cornflowerblue", 
      high = "red",
      
      guide = guide_colorbar(
        title = "logFC",
        order = 2, 
        title.position = "top", title.hjust = 0.5,
        barheight = unit(5, "cm"),
        barwidth  = unit(0.5, "cm"))
    ) +
    labs(title = paste("GO Circos:", title_prefix)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position  = "right",
      legend.direction = "vertical",   
      legend.box       = "vertical",
      legend.margin    = margin(l = 10) 
    )
  
  # 6. 保存
  filename <- paste0("Circos_", title_prefix, ".pdf")
  ggsave(file.path(output_dir, filename), p_final, 
         width = pdf_width, 
         height = pdf_height, 
         dpi = dpi)
  
  cat(sprintf("  [完成] 结果保存在: %s\n", output_dir))
  return(p_final)
}

# 4. 运行示例 ---------------------------------------------------------------

seurat_obj <- readRDS("/home/lin/c_group/hep.rds") 

p2 <- run_go_circos(
  seurat_data = seurat_obj,
  species = "rat",
  scene = "B",
  ident_1 = "Cd1L",                
  ident_2 = "Cd1R",                
  top_n_pathways = 5,
  
  pdf_width = 11.5,      
  pdf_height = 10,   
  dpi = 300,
  gene_font_size = 5,  
  
  output_dir = "/home/lin/c_group/Circos_Results"
)

print(p2)
