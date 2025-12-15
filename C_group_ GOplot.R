#############################################################
# 单细胞 GO Circos 圈图 (v6.0 深度定制图例版)
# 功能：Seurat 差异 -> 数据清洗 -> GOChord -> 深度重构图例
# 特点：
# 1. 完美复刻您提供的右侧垂直图例方案 (方块Terms + 竖条LogFC)。
# 2. 修复 sprintf 报错。
# 3. 保持强制大写逻辑，确保数据一定能画出来。
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
    pdf_width = 14,       # 宽一点给右侧图例
    pdf_height = 10,      
    dpi = 300,            
    gene_font_size = 4,   
    
    output_dir = "./GO_Circos_Result"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("=================================================\n")
  cat(sprintf("开始 GO Circos 分析 | 场景: %s\n", scene))
  
  # --- 1. 数据准备 (保持不变) ---
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
  
  # --- 4. 数据整理 (强制大写匹配，防止报错) ---
  cat("  [INFO] 构建绘图数据...\n")
  go_full_df <- as.data.frame(ego)
  go_top <- go_full_df %>% arrange(p.adjust) %>% head(top_n_pathways)
  
  # 强制大写，确保匹配
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
  
  # --- 5. 绘图 (应用您的定制代码) ---
  cat("  [INFO] 正在绘图与重构图例...\n")
  
  # 定义颜色 (确保和 david$Term 数量一致)
  # 这里我们先定义足够多的颜色，然后按需截取
  all_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  use_cols <- all_cols[1:length(david$Term)]
  names(use_cols) <- david$Term # 命名向量，保证一一对应
  
  # 5.1 生成基础图
  p <- GOChord(chord, 
               space = 0.02,           
               gene.order = 'logFC',   
               gene.space = 0.25,      
               gene.size = gene_font_size,  
               ribbon.col = use_cols   
  )
  
  # 5.2 获取图形信息 (参考您的代码)
  # 这一步是为了获取 logFC 的范围 (limits)
  gb    <- ggplot_build(p)
  
  # 自动寻找含有 fill 属性的 scale (对应 logFC)
  idx <- which(vapply(gb$plot$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1)))
  sc  <- gb$plot$scales$scales[[idx]]
  
  # 取 limits（若为空就用色标的范围）
  # 使用 tryCatch 防止提取失败导致报错
  lims <- tryCatch({ sc$limits %||% range(sc$range$range) }, error = function(e) c(-1, 1))
  brks <- pretty(lims, n = 5)  
  
  # 5.3 深度定制图例 (完全应用您的参考代码逻辑)
  p_final <- p + 
    # 隐藏多余的形状
    guides(shape = "none", 
           # 关键点：劫持 size 图例来制作 GO Terms 图例
           # 因为 GOplot 原图里并没有显式的 fill 图例给 ribbon，我们需要自己造一个
           size  = guide_legend(
             title = "GO Terms", order = 1,
             ncol  = 1, byrow = TRUE,
             # 这里使用我们预定义的 use_cols，确保颜色绝对正确
             override.aes = list(shape = 22, fill = use_cols, size = 8), 
             theme = theme(legend.text = element_text(color="black", size=10))
           )
    ) +
    
    # 重构 logFC 颜色条
    scale_fill_gradient2(
      name   = "logFC",
      limits = lims,
      breaks = brks,
      labels = brks,      
      low = "cornflowerblue", mid = "white", high = "brown1",
      guide = guide_colorbar(
        title = "logFC",
        order = 2, 
        title.position = "top", title.hjust = 0.5,
        barheight = unit(5, "cm"),
        barwidth  = unit(0.5, "cm"))
    ) +
    
    # 主题调整：放置在右侧
    labs(title = paste("GO Circos:", title_prefix)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position  = "right",
      legend.direction = "vertical",   
      legend.box       = "vertical",
      legend.margin    = margin(l = 10) # 左侧留点空隙
    )
  
  # 6. 保存
  filename <- paste0("Circos_", title_prefix, ".pdf")
  ggsave(file.path(output_dir, filename), p_final, 
         width = pdf_width, 
         height = pdf_height, 
         dpi = dpi)
  
  # 修复 sprintf 报错：使用 %s 接收任何类型，或者确保是数字
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
  
  # --- 尺寸与字体 ---
  pdf_width = 12,      # 宽度设置大一点，因为图例在右边很占地方
  pdf_height = 10,   
  dpi = 300,
  gene_font_size = 5,  
  
  output_dir = "/home/lin/c_group/Circos_Results"
)

print(p2)
