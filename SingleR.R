library(Seurat)
library(celldex)
library(SingleR)
library(scrapper)
library(dplyr)
library(ggsci)
library(tidydr)

install.packages("tidydr")

rm(list = ls())
################
## 1. 加载数据
################

pbmc <- readRDS("/home/lin/CRC/CRC_scRNA_harmony.rds")
##获得单细胞表达feature_barcode矩阵
query <- pbmc@assays$RNA$counts
class(query)
#注释细胞名
#Test 输入表达矩阵
#Ref 下载好的参考数据
ref <- celldex::HumanPrimaryCellAtlasData()
#save(ref,file = '/home/lin/CRC/HumanPrimaryCellAtlasData.Rdata')
#seurat V5在GetAssayData代码前添加这段代码
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
write.csv(celltype,"/home/lin/CRC/SingleR_anno.csv")

# 添加celltype到metadata
pbmc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

df=pbmc@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=pbmc@meta.data$celltype)
names(pbmc@reductions)

umap_singleR <- ggplot(df, aes(umap_1, umap_2, color = celltype)) +
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
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = celltype),
               geom = "polygon",
               linetype = 0,
               alpha = 0,
               show.legend = FALSE,
               level = 0.93) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = color)

umap_singleR1 <- umap_singleR + 
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())
print(umap_singleR1)
