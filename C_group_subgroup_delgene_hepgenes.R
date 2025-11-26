# 加载必要的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)

# 清空环境
rm(list = ls())

# 设置工作目录
setwd("/home/lin/c_group/")

# 读取单细胞数据
pbmc <- readRDS("/home/lin/c_group/ident.rds")

# 查看聚类水平
print(levels(pbmc@meta.data$seurat_clusters))

# 定义要分析的细胞类型
celltypes <- c("HSCs", "LESCs")

# 创建要排除的肝细胞相关基因列表
gene <- read_excel("wrgene.xlsx", sheet = "Sheet1")
hepatocyte_genes <- gene$gene
# 为每个细胞类型创建输出目录
dir.create("subpopulation_analysis_filtered_complete", showWarnings = FALSE)

# 对每个细胞类型进行亚群分析
for (celltype in celltypes) {
  # 创建细胞类型特定目录
  cell_dir <- file.path("subpopulation_analysis_filtered_complete", celltype)
  dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(paste0("\nProcessing ", celltype, "...\n"))
  
  # 提取特定细胞类型
  tryCatch({
    # 检查该细胞类型是否存在于聚类中
    if(!celltype %in% levels(pbmc@meta.data$seurat_clusters)) {
      cat("Error: '", celltype, "' not found in seurat_clusters\n")
      cat("Available clusters are:", paste(levels(pbmc@meta.data$seurat_clusters), collapse=", "), "\n")
      next
    }
    
    subset_data <- subset(pbmc, seurat_clusters == celltype)
    
    # 检查子集是否为空
    if(ncol(subset_data) == 0) {
      cat("Error: No cells found for", celltype, "\n")
      next
    }
    
    cat(paste0("Found ", ncol(subset_data), " cells for ", celltype, "\n"))
    
    # 从表达矩阵中删除肝细胞相关基因
    all_genes <- rownames(subset_data)
    genes_to_keep <- setdiff(all_genes, hepatocyte_genes)
    cat(paste0("Removing ", length(all_genes) - length(genes_to_keep), 
               " hepatocyte-related genes from the expression matrix\n"))
    
    # 创建一个新的Seurat对象，只包含非肝细胞基因
    subset_data_filtered <- subset_data[genes_to_keep, ]
    
    cat(paste0("New expression matrix has ", nrow(subset_data_filtered), " genes\n"))
    
    # 数据预处理和降维
    subset_data_filtered <- subset_data_filtered %>%
      NormalizeData() %>%
      ScaleData(vars.to.regress = c('mt_percent')) %>%
      FindVariableFeatures(nfeatures = 2000) %>%
      RunPCA() %>%
      FindNeighbors(reduction = 'pca', dims = 1:50) %>%
      FindClusters(resolution = 0.4) %>%
      RunUMAP(reduction = "pca", dims = 1:50)
    
    # 关键修复：在FindAllMarkers之前先合并图层
    subset_data_filtered <- JoinLayers(subset_data_filtered)
    
    # 寻找差异表达基因
    markers <- FindAllMarkers(
      object = subset_data_filtered,
      only.pos = TRUE,
      logfc.threshold = 0.25
    )
    
    # 保存差异表达基因结果
    if(!is.null(markers) && nrow(markers) > 0) {
      write.csv(
        markers,
        file = file.path(cell_dir, paste0(celltype, "_markers_filtered_complete.csv")),
        quote = FALSE,
        row.names = FALSE
      )
      cat("Saved markers to", file.path(cell_dir, paste0(celltype, "_markers_filtered_complete.csv")), "\n")
    } else {
      cat("No markers found for", celltype, "\n")
    }
    
    # 保存UMAP可视化
    p1 <- DimPlot(
      subset_data_filtered, 
      reduction = "umap", 
      label = TRUE, 
      group.by = "seurat_clusters", 
      pt.size = 0.4
    ) + 
      ggtitle(paste("UMAP of", celltype, "(Completely Filtered)"))
    
    ggsave(
      file.path(cell_dir, paste0(celltype, "_umap_filtered_complete.pdf")),
      p1, 
      width = 10, 
      height = 8, 
      dpi = 300
    )
    cat("Saved UMAP plot to", file.path(cell_dir, paste0(celltype, "_umap_filtered_complete.pdf")), "\n")
    
    # 保存处理后的Seurat对象
    saveRDS(
      subset_data_filtered,
      file = file.path(cell_dir, paste0(celltype, "_processed_filtered_complete.rds"))
    )
    cat("Saved processed Seurat object to", file.path(cell_dir, paste0(celltype, "_processed_filtered_complete.rds")), "\n")
    
    # 输出处理完成信息
    cat(paste0("Finished processing ", celltype, "\n"))
    
  }, error = function(e) {
    cat("Error occurred while processing", celltype, ":\n")
    cat(as.character(e), "\n")
  })
}

# 汇总所有亚群分析结果
cat("\nAll subpopulation analyses completed!\n")
