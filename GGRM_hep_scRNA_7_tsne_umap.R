library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(showtext)
library(sysfonts)

rm(list = ls())

# 函数：创建降维可视化图
create_dim_plot <- function(seurat_obj, 
                            group_by = "seurat_clusters", 
                            reduction = "tsne",
                            output_file = NULL,
                            width = 7, 
                            height = 6) {
  
  # 设置颜色
  color <- c(pal_d3("category20")(20),
             pal_d3("category20b")(20),
             pal_d3("category20c")(20),
             pal_d3("category10")(10))
  
  # 提取降维坐标和分组信息
  df <- seurat_obj@reductions[[reduction]]@cell.embeddings %>%
    as.data.frame() %>%
    cbind(group = seurat_obj@meta.data[[group_by]])
  
  # 获取坐标轴名称
  dim_names <- colnames(df)[1:2]
  
  # 创建基础图
  p <- ggplot(df, aes_string(x = dim_names[1], y = dim_names[2], color = "group")) +
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
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.6, 'cm')) +
    stat_ellipse(aes_string(x = dim_names[1], y = dim_names[2], fill = "group"),
                 geom = "polygon",
                 linetype = 0,
                 alpha = 0,
                 show.legend = FALSE,
                 level = 0.93) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    scale_color_manual(values = color)
  
  # 添加tidydr主题
  p <- p + 
    theme_dr(xlength = 0.22, ylength = 0.22,
             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) +
    theme(panel.grid = element_blank())
  
  # 保存图像
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = width, height = height)
  }
  
  return(p)
}

# 主程序
main <- function() {
  # 加载数据
  pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")
  
  # 创建并保存图像
  p <- create_dim_plot(
    seurat_obj = pbmc,
    group_by = "seurat_clusters",  # 可以改为任何元数据列
    reduction = "tsne",            # 可以改为"umap"
    output_file = "/home/lin/GGRM_hep/result/Ep_tsne_clusters.pdf"
  )
  
  # 显示图像
  print(p)
}

# 运行主程序
main()
