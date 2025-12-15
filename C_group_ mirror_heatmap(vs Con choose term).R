#############################################################
# 单细胞特定通路汇总热图脚本 (v12.1 Fix_GO_Mapping版)
# 修复：org.Rn.eg.db 不包含 TERM 列的问题，引入 GO.db
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
  library(clusterProfiler)
  
  # 注释包
  if(!requireNamespace("org.Rn.eg.db", quietly=TRUE)) {
    BiocManager::install("org.Rn.eg.db")
  }
  library(org.Rn.eg.db)
  
  # 必须加载 GO.db 用于转换通路名
  if(!requireNamespace("GO.db", quietly=TRUE)) {
    BiocManager::install("GO.db")
  }
  library(GO.db)
})

# 2. 定义感兴趣的通路列表 ---------------------------------------------------
# 这里是您指定的通路，保持原样
target_pathways_list <- list(
  "Cd1L" = c("carboxylic acid biosynthetic process", "organic acid biosynthetic process", "fatty acid metabolic process"),
  "Cd1R" = c("sulfur compound metabolic process", "steroid metabolic process", "neutral lipid metabolic process"),
  "Cd15L" = c("cellular detoxification", "immunoglobulin mediated immune response", "sterol biosynthetic process"),
  "Cd15R" = c("adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", "humoral immune response", "B cell mediated immunity"),
  "Cd30L" = c("acyl-CoA metabolic process", "fatty acid metabolic process", "purine nucleotide metabolic process"), 
  "Cd30R" = c("cellular response to interleukin-1", "antigen processing and presentation of exogenous peptide antigen", "complement activation, classical pathway")
)

# 提取所有唯一的通路名称作为Y轴
all_target_terms <- unique(unlist(target_pathways_list))

# 3. 核心计算函数 (修正版) ---------------------------------------------------

calculate_specific_pathway_stats <- function(obj, ident.1, ident.2, pathway_terms) {
  
  # 1. 差异分析
  # 降低阈值以确保热图有颜色填充
  markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, 
                         group.by = "orig.ident", test.use = "wilcox", 
                         logfc.threshold = 0.1, min.pct = 0.05, verbose = FALSE)
  
  if(nrow(markers) == 0) return(NULL)
  markers$gene <- rownames(markers)
  
  # 2. 准备 GO 映射数据
  
  # A. 从 GO.db 获取所有 BP (Biological Process) 的 Term -> GOID 映射
  # keys=keys(GO.db) 获取所有GOID，然后选出 TERM 和 ONTOLOGY
  go_map <- suppressMessages(AnnotationDbi::select(GO.db, 
                                                   keys = keys(GO.db), 
                                                   columns = c("TERM", "ONTOLOGY"), 
                                                   keytype = "GOID"))
  # 只保留 BP (生物学过程)，加快速度并防止重名
  go_map <- go_map %>% filter(ONTOLOGY == "BP")
  
  # B. 从 org.Rn.eg.db 获取 Gene -> GOID 映射
  gene_map <- suppressMessages(AnnotationDbi::select(org.Rn.eg.db, 
                                                     keys = keys(org.Rn.eg.db, keytype="ENTREZID"), 
                                                     columns = c("GOALL", "SYMBOL"), 
                                                     keytype = "ENTREZID"))
  
  results <- data.frame()
  
  for (term in pathway_terms) {
    # 3. 匹配 GO ID (从名字找ID)
    
    # 清理输入字符 (处理可能的特殊连字符)
    term_clean <- str_replace_all(term, "−", "-") 
    
    # 在 GO.db 中查找
    matched_entry <- go_map %>% filter(tolower(TERM) == tolower(term_clean))
    
    if (nrow(matched_entry) == 0) {
      cat(sprintf("  [WARN] Pathway not found in GO.db: %s\n", term))
      # 填0防止报错，但在图中是白色
      results <- rbind(results, data.frame(Description = term, Comparison = ident.1, Score = 0))
      next
    }
    
    target_go_id <- matched_entry$GOID[1] # 如果有重复取第一个
    
    # 4. 匹配基因 (从ID找基因)
    pathway_genes <- gene_map %>% 
      filter(GOALL == target_go_id) %>% 
      pull(SYMBOL) %>% 
      unique()
    
    # 5. 计算分数
    matched_markers <- markers %>% filter(gene %in% pathway_genes)
    
    if (nrow(matched_markers) < 2) {
      final_score <- 0
    } else {
      # 计算逻辑：平均LogFC方向 * 显著性强度
      # P值极小值处理，防止 Inf
      p_vals <- matched_markers$p_val_adj
      p_vals[p_vals == 0] <- 1e-300
      
      avg_logfc <- mean(matched_markers$avg_log2FC)
      # 使用中位数P值代表通路整体显著程度
      sig_score <- -log10(median(p_vals))
      
      final_score <- sign(avg_logfc) * sig_score
    }
    
    results <- rbind(results, data.frame(
      Description = term,
      Comparison = ident.1, 
      Score = final_score
    ))
  }
  
  return(results)
}

# 4. 主流程 -----------------------------------------------------------------

run_selected_pathway_heatmap <- function(
    seurat_file = "/home/lin/c_group/hep.rds",
    output_dir = "/home/lin/c_group/Summary_Heatmap"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (is.character(seurat_file)) {
    obj <- readRDS(seurat_file)
  } else {
    obj <- seurat_file
  }
  Idents(obj) <- obj@meta.data$orig.ident
  
  comparisons <- c("Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R")
  control_grp <- "Con"
  
  # 检查 Control 组
  if (!control_grp %in% levels(Idents(obj))) {
    stop("Error: 'Con' group not found in Seurat object.")
  }
  
  all_results <- data.frame()
  
  cat("开始计算指定通路的富集分数...\n")
  
  for (comp in comparisons) {
    cat(sprintf("  -> Processing %s vs %s ...\n", comp, control_grp))
    
    if (!comp %in% levels(Idents(obj))) {
      cat(sprintf("     [SKIP] Group %s not found.\n", comp))
      next
    }
    
    res <- calculate_specific_pathway_stats(obj, ident.1 = comp, ident.2 = control_grp, 
                                            pathway_terms = all_target_terms)
    
    if (!is.null(res) && nrow(res) > 0) {
      all_results <- rbind(all_results, res)
    }
  }
  
  # 5. 绘图准备 ---------------------------------------------------------------
  
  if (nrow(all_results) == 0) stop("No results generated. Check gene names or groups.")
  
  # X轴顺序
  all_results$Comparison <- factor(all_results$Comparison, levels = comparisons)
  
  # Y轴顺序 (反转列表，让第一组排在最上面，或者按ggplot默认从下往上)
  # 这里我们让Cd1L的通路排在最下面，Cd30R的排在最上面
  final_y_levels <- unique(unlist(target_pathways_list)) # 保持列表顺序
  all_results$Description <- factor(all_results$Description, levels = final_y_levels)
  
  # 颜色截断
  limit_val <- 10 # 这里的P值往往比较显著，设为10或者15看起来对比度好
  all_results$PlotScore <- all_results$Score
  all_results$PlotScore[all_results$PlotScore > limit_val] <- limit_val
  all_results$PlotScore[all_results$PlotScore < -limit_val] <- -limit_val
  
  # 6. 绘图 -------------------------------------------------------------------
  cat("  -> Plotting Summary Heatmap...\n")
  
  p <- ggplot(all_results, aes(x = Comparison, y = Description, fill = PlotScore)) +
    geom_tile(color = "white", size = 0.5) +
    
    scale_fill_gradient2(
      low = "#3288BD", 
      mid = "white", 
      high = "#D53E4F", 
      midpoint = 0,
      name = "Score",
      limits = c(-limit_val, limit_val),
      oob = scales::squish
    ) +
    
    labs(
      title = "Selected Pathway Enrichment (Exp vs Con)",
      subtitle = "Red = Up in Exp (CdX) | Blue = Up in Con",
      x = NULL, y = NULL
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
  
  out_file <- file.path(output_dir, "Summary_Selected_Pathways_Heatmap_Fixed.pdf")
  ggsave(out_file, p, width = 9, height = 11)
  
  cat(sprintf("[DONE] Result saved to: %s\n", out_file))
}

# 运行
seurat_obj <- readRDS("/home/lin/c_group/hep.rds")
run_selected_pathway_heatmap(seurat_obj)
