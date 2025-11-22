library(celda)
library(Seurat)
library(cowplot)

rm(list = ls())

scobj <- readRDS("/home/lin/c_group/C_group_harmony.rds")

#可视化污染情况
#seurat V5在GetAssayData代码前添加这段代码
scobj = JoinLayers(scobj)
counts <- scobj@assays$RNA@layers$counts
# 计算
# 需要把结果储存在新变量
decontX_results <- decontX(counts) 
# 你可以使用str()查看结果的形式
# RNA污染的计算结果储存在：decontX_results$contamination
# 他是一个数值型向量，长度与细胞数量一致，顺序也与矩阵的colnames一致
# 我们可以直接把他写入metadata
scobj$Contamination =decontX_results$contamination
head(scobj@meta.data) # 可以查看一下结果

# 我们对他进行可视化
library(ggplot2)
FeaturePlot(scobj, 
            features = 'Contamination', 
            raster=FALSE       # 细胞过多时候需要加这个参数
) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')

## 在这里看到，对比前文最开始的状态，我们经过了严格的质控之后，散在的细胞已经少了很多，但是仍然有亚群之间不清晰的界限
## 我们可视化Contamination，惊奇地发现这些边缘的毛躁就是Contamination较高地细胞

low_con_scobj = scobj[,scobj$Contamination < 0.3] #保留小于等于0.3的细胞

# 计算每个样本细胞删除的比例
sample_counts <- table(scobj$orig.ident)
print(paste0("细胞删除比例: ", round((1 - sum(table(low_con_scobj$orig.ident))/sum(sample_counts))*100, 2), "%, 各样本删除比例: ", 
             paste(sapply(names(sample_counts), function(s) paste0(s, ": ", round((1 - table(low_con_scobj$orig.ident)[s]/sample_counts[s])*100, 2), "%")), collapse=", ")))


## 查看去污染后的结果
p1 <- DimPlot(scobj, reduction = "umap", group.by = "orig.ident", pt.size=0.4) + ggtitle("before decoutX")
p2 <- DimPlot(low_con_scobj, reduction = "umap", group.by = "orig.ident", pt.size=0.4, label=T, repel=T) + ggtitle("after decoutX")

plot_grid(p1, p2, ncol = 2)

# 保存RDS文件
saveRDS(low_con_scobj, file = "/home/lin/c_group/C_group_decontX(0.3)_11m14d.rds")
