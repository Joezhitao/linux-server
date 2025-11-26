library(Seurat)

# 读取已处理的Seurat对象
pbmc <- readRDS("/home/lin/c_group/T_ident.rds")
levels(pbmc@meta.data$seurat_clusters)

# 筛选特定的细胞类型
clusters <- c("CD4+ Memory T cells (with gut-homing potential)", "CD4+ Naive/Central Memory T cells", 
              "CD4+ Regulatory T cells")
pbmc <- pbmc[, pbmc@meta.data$seurat_clusters %in% clusters]
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

saveRDS(pbmc,"/home/lin/c_group/CD4+T_ident.rds")
