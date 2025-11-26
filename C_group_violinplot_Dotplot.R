library(Seurat)

setwd("/home/lin/c_group/")

rm(list = ls())

pbmc <- readRDS("/home/lin/c_group/HSCs_ident.rds")
pbmc@meta.data$seurat_clusters
# 设置分组变量的顺序
desired_order <- c("Ad1R", "Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R")
pbmc@meta.data$orig.ident <- factor(pbmc@meta.data$orig.ident, levels = desired_order)

VlnPlot(object = pbmc, features = 'Tgfb2')

gene <- c("Il10")

plot <- DotPlot(pbmc, features = gene,group.by = 'orig.ident',dot.scale = 30)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
  scale_y_discrete(limits = rev(levels(pbmc@meta.data$orig.ident)))

print(plot)
