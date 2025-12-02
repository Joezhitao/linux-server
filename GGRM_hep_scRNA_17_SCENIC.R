#!/usr/bin/env Rscript

# ==============================================================================
# 1. 环境配置与包加载
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(GENIE3)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

rm(list = ls())

# 设置工作目录与并行核心数
work_dir <- "/home/lin/GGRM_hep/result/SCENIC_analysis"
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)
dir.create("int", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
nCores_set <- 20 

print(paste("Start Time:", Sys.time()))

# ==============================================================================
# 2. 数据加载与预处理
# ==============================================================================
print(">>> Loading Data...")
seurat_obj <- readRDS("/home/lin/GGRM_hep/result/cluster_optimization/filtered_both_threshold0.2.rds")
exprMat <- as.matrix(GetAssayData(seurat_obj, layer = "counts"))
cellInfo <- seurat_obj@meta.data[, c("orig.ident", "group", "cell_state")]

# ==============================================================================
# 3. 加载motif注释和初始化SCENIC
# ==============================================================================
print(">>> Loading motif annotations and initializing SCENIC...")

# 加载motif注释
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# 使用您指定的v1版本数据库文件
dbDir <- "/home/lin/cisTarget_db/SCENIC" 
dbs <- c(
  '500bp' = 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather',
  '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'
)

# 初始化SCENIC
scenicOptions <- initializeScenic(org = "hgnc", 
                                  nCores = nCores_set,
                                  datasetTitle = "Liver_Organoid_Resistance",
                                  dbDir = dbDir,
                                  dbs = dbs)

# 保存元数据
saveRDS(cellInfo, file="int/cellInfo.Rds")

# ==============================================================================
# 4. 基因过滤与共表达网络推断 (GENIE3)
# ==============================================================================
print(">>> Step 1: Gene Filtering & Co-expression Network (GENIE3)...")

# 过滤低表达基因
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

# 计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)

# 运行GENIE3
exprMat_log <- log2(exprMat_filtered + 1)
runGenie3(exprMat_log, scenicOptions)

# ==============================================================================
# 5. 构建Regulon与AUCell评分
# ==============================================================================
print(">>> Step 2 & 3: Build Regulons & Score Cells...")

# 推断共表达模块
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# 推断转录调控网络
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

# 对细胞进行评分
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, log2(exprMat + 1))

# 二值化活性状态
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# 保存结果
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# ==============================================================================
# 6. 下游分析与可视化
# ==============================================================================
print(">>> Step 4: Visualization & Downstream Analysis...")

# 提取AUC矩阵
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
auc_matrix <- getAUC(regulonAUC)

# 比较不同组别的Regulon活性
avg_activity_group <- sapply(split(rownames(cellInfo), cellInfo$group),
                             function(cells) rowMeans(auc_matrix[, cells]))
scaled_activity_group <- t(scale(t(avg_activity_group), center = TRUE, scale = TRUE))

# 绘制热图 - 按组别
pdf("output/Heatmap_Regulon_Activity_by_Group.pdf", width = 8, height = 12)
Heatmap(scaled_activity_group, 
        name = "Regulon Activity\n(Z-score)",
        cluster_rows = TRUE, 
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_title = "Regulon Activity: Resistant vs Non-resistant")
dev.off()

# 比较不同细胞状态的Regulon活性
avg_activity_state <- sapply(split(rownames(cellInfo), cellInfo$cell_state),
                             function(cells) rowMeans(auc_matrix[, cells]))
scaled_activity_state <- t(scale(t(avg_activity_state), center = TRUE, scale = TRUE))

pdf("output/Heatmap_Regulon_Activity_by_CellState.pdf", width = 10, height = 14)
Heatmap(scaled_activity_state, 
        name = "Regulon Activity\n(Z-score)",
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_title = "Regulon Activity across Cell States")
dev.off()

# 计算RSS并绘制点图
rss_group <- calcRSS(AUC=auc_matrix, cellAnnotation=cellInfo[colnames(auc_matrix), "group"])
rss_group_data <- rss_group$RSS

dot_plot_data <- data.frame()
for(group in colnames(scaled_activity_group)) {
  group_data <- data.frame(
    TF = rownames(scaled_activity_group),
    Group = group,
    Zscore = scaled_activity_group[,group],
    RSS = rss_group_data[rownames(scaled_activity_group), group]
  )
  dot_plot_data <- rbind(dot_plot_data, group_data)
}

library(ggplot2)
pdf("output/DotPlot_TF_Activity_by_Group.pdf", width = 10, height = 8)
ggplot(dot_plot_data, aes(x = Group, y = TF, color = Zscore, size = RSS)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Transcription Factor Activity by Group",
       x = "", y = "", color = "Z-score", size = "RSS")
dev.off()

print(paste("All Done! End Time:", Sys.time()))
