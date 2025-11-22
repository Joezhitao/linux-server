# 加载必要的R包
library(Seurat)
library(dplyr)

rm(list = ls())
# 设置工作目录（根据实际情况调整）
setwd("/home/lin/c_group/")

# 读取之前处理好的单细胞数据
pbmc <- readRDS("/home/lin/c_group/optimization_results/umap/filtered_all_clusters_threshold0.2.rds")
pbmc <- JoinLayers(pbmc)
combined_markers <- FindAllMarkers(object = pbmc, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)

all <- combined_markers %>% group_by(cluster)

write.csv(all, 
          file = "/home/lin/c_group/combined_all_markers(resolution=0.3).csv", 
          quote = FALSE, 
          row.names = FALSE)
pdf()
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE,group.by = "seurat_clusters", pt.size=0.4) + ggtitle("By umap")

ggsave("/home/lin/c_group/optimization_results/umap/filtered_all_clusters_threshold0.2.pdf", 
       p1, width = 12, height = 6, dpi = 300)
