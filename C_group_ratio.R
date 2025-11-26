library(Seurat)
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(dplyr)

setwd("/home/lin/c_group/")
pbmc <- readRDS("/home/lin/c_group/ident.rds")

# 确保细胞类型和分组是因子，并按想要顺序排列
pbmc@meta.data$seurat_clusters <- factor(pbmc@meta.data$seurat_clusters,
                                         levels = c("Hepatocytes", "T cells", "NK/NKT cells", "LESCs", "B cells",
                                                    "Neutrophils", "Kupffer cells", "HSCs",
                                                    "Inflammatory Monocyte-Derived Macrophages", "Dendritic Cells",
                                                    "Cycling cells", "Hepatic Progenitor Cells")
)
pbmc@meta.data$orig.ident <- factor(pbmc@meta.data$orig.ident,
                                    levels = c("Ad1R", "Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R")
)
levels(pbmc@meta.data$orig.ident)

# 统计细胞类型和分组的数量
cell_count <- table(pbmc@meta.data$seurat_clusters, pbmc@meta.data$orig.ident)

# 按分组计算比例
cell_ratio <- prop.table(cell_count, margin = 2) %>%
  as.data.frame() %>%
  set_names(c("celltype", "group", "ratio"))
write_csv(cell_ratio,"/home/lin/c_group/ratio.csv")
#自定义配色
colour <- c(
  "Hepatocytes" = "#e67e22", 
  "T cells" = "#2980b9",
  "NK/NKT cells" = "#16a085", 
  "LESCs" = "#636e72",
  "B cells" = "#fbc531",
  "Neutrophils" = "#00b894",
  "Kupffer cells" = "#00cec9",
  "HSCs" = "#d35400",
  "Inflammatory Monocyte-Derived Macrophages" = "#8e44ad",
  "Dendritic Cells" = "#2ecc71",
  "Cycling cells" = "#e84393",
  "Hepatic Progenitor Cells" = "#fdcb6e"
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
print(p5)
