# 加载必要的R包
library(Seurat)
library(dplyr)

# 清空环境
rm(list = ls())

# 设置工作目录
setwd("/home/lin/GGRM_hep/result/")

# 1. 读取数据
pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_ident_11m27d.rds")

# 检查一下当前的 Idents 是否是你想要的 seurat_clusters
# 如果不是，显式设定一下：
Idents(pbmc) <- "seurat_clusters"

# 2. 如果你的Seurat版本是V5，且之前做过Integration，这步 JoinLayers 很有必要
# 它可以确保你在做 FindAllMarkers 时用的是整合后的 counts/data，或者是合并后的 RNA层
pbmc <- JoinLayers(pbmc)

# 3. 计算所有簇的 Marker 基因
# 专家建议调整：
# min.pct = 0.25: 基因至少在25%的细胞中表达（默认0.1有点太低，容易找到噪声）
# logfc.threshold = 0.25: 差异倍数阈值
combined_markers <- FindAllMarkers(object = pbmc, 
                                   only.pos = TRUE,      # 只找上调基因
                                   min.pct = 0.25,       # 建议加上这个过滤，减少运算量且结果更准
                                   logfc.threshold = 0.25)

# 4. 提取每个簇的 Top 25 基因 (根据 avg_log2FC 排序)
top25_markers <- combined_markers %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC) # 推荐用 slice_max，比 top_n 更明确

# 5. 查看结果
print(head(top25_markers))

# 6. 保存结果 (强烈建议保存，方便后续做热图或GO分析)
write.csv(top25_markers, "Top25_Markers_per_Cluster.csv", row.names = FALSE)
write.csv(combined_markers, "All_Markers_List.csv", row.names = FALSE)

# ==========================================
# 热图看看 Top 10 效果如何
# ==========================================
# 取 Top 10 画图更清晰
top10 <- combined_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()
