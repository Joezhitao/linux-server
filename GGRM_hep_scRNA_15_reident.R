# ==============================================================================
# 0. 环境准备与设置
# ==============================================================================
# 加载进行单细胞分析所需的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork) 

# 清空当前环境
rm(list = ls())

# 设置工作目录
setwd("/home/lin/GGRM_hep/result/")

# ==============================================================================
# 1. 数据读取与原始分布回顾
# ==============================================================================
# 读取之前处理好的单细胞RDS文件
pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_ident_11m27d.rds")

# 确认当前的聚类ID是原始的数字编号
Idents(pbmc) <- "seurat_clusters"

# 先画一张原始的UMAP图，回顾一下之前的13个簇的分布情况
# 这步用于确认之前观察到的拓扑结构（1/3在上，2/4在左下等）
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# ==============================================================================
# 2. 细胞亚群重定义 (基于观察到的拓扑结构与生物学功能)
# ==============================================================================
# 根据在UMAP图上观察到的空间分布（Topology），结合Marker基因功能，
# 决定将这13个簇合并为5个具有明确生物学意义的“状态（States）”。

# 合并策略如下：
# 1. 右下方延伸区 (0, 6, 8, 11): 这是一个混合的大群，兼具祖细胞特征和代谢活性。
# 2. 上方区域 (1, 3): 这两群连在一起，表现为分化较好的上皮和分泌特征。
# 3. 左下方独立岛 (2, 4): 这群细胞离得最远，且富集EMT和干性基因，将其定义为耐药相关的核心群。
# 4. 中间连接区 (5, 10, 12): 这几群聚在一起，主要表现为应激反应，可能是状态过渡。
# 5. 底部独立区 (7, 9): 这是纯粹的增殖细胞群。

my_topology_ids <- c(
  # --- 区域 3: 右下方混合区 ---
  "0" = "Progenitor-Metabolic", 
  "6" = "Progenitor-Metabolic", 
  "8" = "Progenitor-Metabolic", 
  "11" = "Progenitor-Metabolic",
  
  # --- 区域 1: 上方区域 ---
  "1" = "Epithelial-Secretory",
  "3" = "Epithelial-Secretory",
  
  # --- 区域 4 (岛屿B): 左下方独立岛 (耐药核心) ---
  "2" = "Mesenchymal-Stem-like", 
  "4" = "Mesenchymal-Stem-like", 
  
  # --- 区域 2: 中间连接区 ---
  "5" = "Stress-Adaptive",
  "10" = "Stress-Adaptive", 
  "12" = "Stress-Adaptive",
  
  # --- 区域 4 (岛屿A): 底部独立区 ---
  "7" = "Cycling",
  "9" = "Cycling"
)

# 将上述定义好的名称应用到Seurat对象中
pbmc <- RenameIdents(pbmc, my_topology_ids)

# 将新的细胞状态注释保存到 metadata 中，名为 'cell_state'
pbmc$cell_state <- Idents(pbmc)

# ==============================================================================
# 3. (可选步骤) 数据清洗
# ==============================================================================
# 如果认为 Stress-Adaptive 这群细胞主要是死细胞干扰，可以选择剔除它们。
# 目前先保留，视后续图形效果决定是否运行下面这行代码：
# pbmc <- subset(pbmc, idents = "Stress-Adaptive", invert = TRUE)

# ==============================================================================
# 4. 最终图谱可视化
# ==============================================================================
# 定义了一套符合 SCI 发表标准的配色方案
# 特意用红色标记 Mesenchymal-Stem-like，因为这是关注的耐药重点
my_colors <- c(
  "Epithelial-Secretory" = "#3C5488CC",     # 深蓝
  "Progenitor-Metabolic" = "#4DBBD5CC",     # 浅蓝
  "Mesenchymal-Stem-like" = "#E64B35CC",    # 红色 (重点高亮)
  "Stress-Adaptive" = "#00A087CC",          # 绿色
  "Cycling" = "#F39B7FCC"                   # 橙色
)

# 绘制最终的细胞状态 UMAP 图
p_umap <- DimPlot(pbmc, reduction = "umap", 
                  cols = my_colors,
                  label = TRUE,       
                  label.size = 4,     
                  repel = TRUE,       
                  pt.size = 0.8) +    
  ggtitle("Cellular States Landscape") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 将图片保存为 PDF，以便后续插入文章
ggsave("02_Final_CellStates_UMAP.pdf", plot = p_umap, width = 8, height = 6)

# ==============================================================================
# 5. 组间差异验证 (耐药 vs 敏感)
# ==============================================================================
# 需要验证在耐药样本中，是否特定的细胞状态（如 Mesenchymal-Stem-like）比例增加了。

# 绘制堆叠柱状图来展示各样本的细胞组成
p_bar <- ggplot(pbmc@meta.data, aes(x = orig.ident, fill = cell_state)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = my_colors) +
  labs(y = "Proportion", x = "Sample", title = "Cell State Composition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 同时绘制拆分的 UMAP 图，直观对比不同组别的分布
p_split <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident", cols = my_colors, ncol = 4) + NoLegend()

# 将两张图拼在一起并保存
final_plot <- p_split / p_bar
ggsave("03_Composition_Analysis.pdf", plot = final_plot, width = 12, height = 10)

# 最后，保存处理好的 Seurat 对象
saveRDS(pbmc, "GGRM_hep_Final_Annotated.rds")

# 提示：分析完成
print("分析已完成，结果已保存。")
