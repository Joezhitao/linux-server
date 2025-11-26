library(Seurat)
library(clustree)
library(tidydr)

rm(list = ls())

setwd("/home/lin/GGRM_hep/result/")

pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_decontX_11m26d.rds")
levels(pbmc@meta.data$orig.ident)
################
## 1. 不同分辨率聚类树可视化
################

res_values <- seq(0.1, 1.2, by=0.1)
for(res in res_values){
  pbmc <- FindClusters(pbmc, resolution = res)
}
clustree(pbmc)
pbmc <- FindClusters(pbmc, resolution = 0.6)

DimPlot(pbmc, reduction = "tsne", label = TRUE, label.size = 5) + NoLegend()

################
## 2. SingleR细胞鉴定
################
library(SingleR)
library(ggsci)
library(dplyr) #管道符
##获得单细胞表达feature_barcode矩阵
query <- pbmc@assays$RNA$counts
class(query)
#注释细胞名
#Test 输入表达矩阵
#Ref 下载好的参考数据
#ref <- celldex::HumanPrimaryCellAtlasData()
#save(ref,file = '/Users/joe/Desktop/GG_liver_scRNA_seq/HumanPrimaryCellAtlasData.Rdata')
#直接加载之前下好的鉴定参考数据
load('/home/lin/GGRM_hep/HumanPrimaryCellAtlasData.Rdata')
#seurat V5在GetAssayData代码前添加这段代码
ref$label.fine
pbmc = JoinLayers(pbmc)
pbmc_for_SingleR <- GetAssayData(pbmc, layer="data") # 标准化后数据
clusters <- pbmc@meta.data$seurat_clusters
pred <- SingleR(test = pbmc_for_SingleR, 
                ref = ref, # 参考数据
                labels = ref$label.fine, #标签列 当'cluster='参数被指定时，method="cluster"参数可以缺省
                clusters = clusters, # 为cluster注释celltype
                assay.type.test = "logcounts", 
                assay.type.ref = "logcounts")

head(pred)
table(pred$labels)

#保存SingleR注释结果
celltype <- data.frame(ClusterID = rownames(pred), 
                       celltype = pred$labels, 
                       stringsAsFactors = F)
write.csv(celltype,"/home/lin/GGRM_hep/result/SingleR_anno.csv")

# 添加celltype到metadata
pbmc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

df=pbmc@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=pbmc@meta.data$celltype)
names(pbmc@reductions)
head(df)
tsne_singleR <- ggplot(df, aes(tSNE_1, tSNE_2, color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = color) +
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 14, family = "myFont1"),
        legend.key.size = unit(0.6, 'cm'),
        text = element_text(family = "myFont1")) +
  stat_ellipse(aes(x = tSNE_1, y = tSNE_2, fill = celltype),
               geom = "polygon",
               linetype = 0,
               alpha = 0,
               show.legend = FALSE,
               level = 0.93) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = color)

tsne_singleR1 <- tsne_singleR + 
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())
pdf("singleR.pdf", width = 7, height = 6)
print(tsne_singleR1)
dev.off()
pbmc@meta.data$celltype <- as.factor(pbmc@meta.data$celltype)
levels(pbmc@meta.data$celltype) <- c("Epithelial_cells", "Epithelial_cells", "Hepatocytes", "MSC")
levels(pbmc@meta.data$celltype)
# 保存RDS文件
saveRDS(pbmc, file = "/home/lin/GGRM_hep/result/GGRM_hep_singleR_11m26d.rds")

DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident", label = TRUE, label.size = 5) + NoLegend()
pbmc@meta.data$orig.ident <- as.factor(pbmc@meta.data$orig.ident)
pbmc@meta.data$group <- pbmc@meta.data$orig.ident
DimPlot(pbmc, reduction = "tsne", group.by = "group", label = TRUE, label.size = 5) + NoLegend()
