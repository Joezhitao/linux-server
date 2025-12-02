library(CytoTRACE2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggpubr)

rm(list = ls())

pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")

sce2 <- pbmc
sce2@meta.data$seurat_clusters
sce2@meta.data$CB <- rownames(sce2@meta.data)
sample_CB <- sce2@meta.data %>% 
  group_by(seurat_clusters) %>% 
  sample_frac(0.3)
sce3 <- subset(sce2,CB %in% sample_CB$CB) 
sce3
# An object of class Seurat 

#######输入seurat 对象###########
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = sce3@meta.data$seurat_clusters) %>% 
  set_rownames(., colnames(sce3))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

png("cytotrace2_4.png", width = 12, height = 12, units = "in", res = 1000)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
dev.off()

png("cytotrace2.png", width = 5, height = 5, units = "in", res = 1000)
FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1.5) + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)
dev.off()

library(ggpubr)
p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x="seurat_clusters", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="seurat_clusters",#填充
                palette = "npg",
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") #图例放右边 
###指定组比较
factors <- c(levels(sce2@meta.data$seurat_clusters))

combinations <- combn(factors, 2) 
my_comparisons <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))

pdf("/home/lin/GGRM_hep/result//cytotrace2_5.pdf", width = 8, height = 5)
p1+stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")
dev.off()