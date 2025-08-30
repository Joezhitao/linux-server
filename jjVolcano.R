library(Seurat) 
library(tidyverse)
library(scRNAtoolVis)
library(schard)  # 用于h5ad到Seurat的转换
library(xlsx)

setwd("/home/lin/gfp")

rm(list = ls())

# 读取已处理的Seurat对象
pbmc <- schard::h5ad2seurat('gfp.h5ad')
pbmc@meta.data$cell_type <- as.factor(pbmc@meta.data$cell_type)
levels(pbmc@meta.data$cell_type)

# 创建新的分组变量
pbmc@meta.data$merged_groups <- as.character(pbmc@meta.data$orig.ident)

# 将"RT"和"RT+HT"合并为一个新组，例如命名为"RT_combined"
pbmc@meta.data$merged_groups[pbmc@meta.data$merged_groups %in% c("RT", "RT+HT")] <- "RT"

# 转换为因子并设置因子水平顺序（如果需要特定顺序）
pbmc@meta.data$merged_groups <- factor(pbmc@meta.data$merged_groups, 
                                       levels = c("Con", "HT", "RT"))
# 检查新的分组
table(pbmc@meta.data$cell_type)

group <- levels(pbmc@meta.data$merged_groups)

#group去除Con
group <- group[group != "Con"]

cluster <- levels(pbmc@meta.data$cell_type)
# 首先设置当前的识别标签为 cell_type
Idents(pbmc) <- "cell_type"

for (j in group) {
  combined_table <- data.frame()
  
  for (i in cluster) {
    # 创建一个临时的 Seurat 对象，只包含当前细胞类型的细胞
    temp_pbmc <- subset(pbmc, idents = i)
    
    # 在这个临时对象上设置 merged_groups 作为识别标签
    Idents(temp_pbmc) <- "merged_groups"
    
    # 检查是否包含所需的组别
    available_groups <- levels(factor(temp_pbmc@meta.data$merged_groups))
    if(!(j %in% available_groups) || !("Con" %in% available_groups)) {
      message(paste("跳过", i, "细胞类型，因为缺少", j, "或 Con 组的细胞"))
      next
    }
    
    # 检查每个组的细胞数量是否足够
    cell_counts <- table(Idents(temp_pbmc))
    if(cell_counts[j] < 3 || cell_counts["Con"] < 3) {
      message(paste("跳过", i, "细胞类型，因为", j, "组或 Con 组的细胞数量少于3个"))
      next
    }
    
    # 运行 FindMarkers
    tryCatch({
      pbmc.markers <- FindMarkers(temp_pbmc, 
                                  ident.1 = j, 
                                  ident.2 = "Con", 
                                  logfc.threshold = 0.25, 
                                  min.pct = 0.1)
      
      # 处理结果
      if(nrow(pbmc.markers) > 0) {
        pbmc.markers <- pbmc.markers %>% 
          rownames_to_column(var = "gene") %>%
          mutate(cluster = i, group = j)
        
        combined_table <- bind_rows(combined_table, pbmc.markers)
      }
    }, error = function(e) {
      message(paste("处理", i, "细胞类型时出错:", e$message))
    })
  }
  
  # 保存结果
  if(nrow(combined_table) > 0) {
    write.csv(combined_table, 
              file = paste0("/home/lin/gfp/火山图/", j, "_combined.csv"), 
              row.names = FALSE)
  }
}

file = "/home/lin/gfp/火山图/RT_combined.csv"
combined_markers = read.csv(file,header = T,stringsAsFactors = F)
path = paste("/home/lin/gfp/火山图/Volcano",".png", sep = "")


color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
png(path,units = "in",width = 30,height = 20,res = 600)
jjVolcano(diffData = combined_markers,
          log2FC.cutoff = 0.1, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#87CEFA','#EEA2AD'), #设置点的颜色
          tile.col = color.pals, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 20 #设置展示topN的基因
)
dev.off()

png(path, units = "in", width = 30, height = 20, res = 600)
jjVolcano(diffData = combined_markers,
          log2FC.cutoff = 0.1, 
          size = 3.5,                # 点的大小
          fontface = 'italic',       # 字体形式
          aesCol = c('#87CEFA','#EEA2AD'), # 点的颜色
          tile.col = color.pals,     # 簇的颜色
          topGeneN = 20,             # 展示top N的基因
          
          # 增加标签文字大小
          label.size = 5,            # 基因标签的大小
          celltypeSize = 5,          # 细胞类型标签的大小
          
          # 调整主题文字大小
          base_size = 18,            # 整体主题的基础字体大小
          
          # 调整细胞簇框的大小
          expand = c(-1.5, 1.5)      # 增加y轴扩展范围，使框更大
)
dev.off()
