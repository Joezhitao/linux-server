rm(list = ls())
# 加载包
library("Seurat")
library("Nebulosa")
library("ggplot2") # 加载ggplot2用于ggsave

# 设置工作目录和读取数据
setwd("/home/lin/c_group/subpopulation_analysis")

# 读取数据
pbmc <- readRDS("/home/lin/c_group/hep.rds")

gene_list <- c(
  # 列表 1
  "Cyp2c13", "Cyp2a1", "Vegfa", "Igfbp1", "Cyp27a1", 
  "Cyp39a1", "Cyp4a2", "Cyp2c11", "Srebf1", "Fasn", 
  "Hint1", "Cox7c", "Apoc1", "Fabp1", "Pck1", 
  "Ass1", "Ag1", "Uox", "Mfsd2a", "Hmgcr", 
  "Hnf4a", "Acly", "G6pc", "Igf1", "Cps1",
  
  # 列表 2 (Wnt/Zonation markers 等)
  "Axin2", "Rspo3", "Lgr5", "Cyp2f2", "Oit3", 
  "Cyp1a2", "Sds", "Hal", "Adh1", "Gulo", 
  "Pxmp2"
)

# 获取Seurat对象中所有的基因名，用于检查是否存在
all_genes <- rownames(pbmc)

# 开始循环
for (g in gene_list) {
  
  # 1. 检查基因是否存在于数据中
  if (g %in% all_genes) {
    
    message(paste("正在绘制:", g))
    
    # 2. 绘制单个基因的密度图
    p <- plot_density(
      pbmc, 
      features = g, 
      reduction = "umap", 
      size = 0.8  # 保持你设置的点大小
    )
    
    # 修改标题（可选，让标题更简洁）
    p <- p + ggtitle(label = paste0(g, " Density"))
    
    # 3. 保存图片
    # 调整了width和height，因为是单张图，不需要像m3那样宽
    ggsave(
      plot = p,
      filename = paste0(g, "_density.png"), # 文件名例如: Me1_density.png
      width = 6,      
      height = 5, 
      dpi = 600,
      bg = "white"
    )
    
  } else {
    message(paste("警告: 基因", g, "未在 Seurat 对象中找到，已跳过。"))
  }
}

message("所有绘图已完成！")
