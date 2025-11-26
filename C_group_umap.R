library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(showtext)
library(showtext)
library(sysfonts)

# 加载合并的数据集
pbmc <- readRDS("/home/lin/c_group/ident.rds")
levels(pbmc@reductions$umap)
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

df=pbmc@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(seurat_clusters=pbmc@meta.data$seurat_clusters)

head(df)

p_tsne1 <- ggplot(df, aes(umap_1, umap_2, color = seurat_clusters)) +
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
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = seurat_clusters),
               geom = "polygon",
               linetype = 0,
               alpha = 0,
               show.legend = FALSE,
               level = 0.93) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = color)

p_tsne2 <- p_tsne1 + 
  theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())

png("E:/B组数据备份(4.29)/单细胞结果/UMAP图美化/umap_1.2.png", width = 12, height = 8, units = "in", res = 1000)
print(p_tsne2)
dev.off()