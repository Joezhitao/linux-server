library(Seurat)
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(dplyr)

rm(list = ls())

setwd("/home/lin/GGRM_hep/result/")
pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_ident_11m27d.rds")

# 确保细胞类型和分组是因子，并按想要顺序排列
pbmc@meta.data$seurat_clusters <- factor(pbmc@meta.data$celltype,
                                         levels = c("Epithelial_cells", "Hepatocytes", "MSC")
)
pbmc@meta.data$orig.ident <- factor(pbmc@meta.data$orig.ident,
                                    levels = c("Org1", "Org2", "Org3", "Org4")
)
levels(pbmc@meta.data$group)

# 统计细胞类型和分组的数量
cell_count <- table(pbmc@meta.data$celltype, pbmc@meta.data$group)

# 按分组计算比例
cell_ratio <- prop.table(cell_count, margin = 2) %>%
  as.data.frame() %>%
  set_names(c("celltype", "group", "ratio"))
write_csv(cell_ratio,"/home/lin/c_group/ratio.csv")
#自定义配色
colour <- c(
  "Epithelial_cells" = "#e67e22", 
  "Hepatocytes" = "#2980b9",
  "MSC" = "#16a085"
)
cell_ratio$celltype <- factor(cell_ratio$celltype, levels = names(colour))

#绘制冲击图 (Alluvial plot)
# 基础ggplot对象
pp <- ggplot(cell_ratio, aes(x = group, y = ratio, fill = celltype,
                             stratum = celltype, alluvium = celltype)) +
  scale_fill_manual(values = colour) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()

# 冲击图（直线连线效果）
p5 <- pp +
  geom_col(width = 0.6, color = NA, linewidth = 0.5) +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0, color = 'white', linewidth = 0.5) +
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0, fill = NA, color = 'white', linewidth = 0.5) +
  labs(x = 'Group', y = 'Cell Ratio', fill = 'Cell Type') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

# 显示图
pdf("ratio.pdf", width = 4, height = 6)
print(p5)
dev.off()