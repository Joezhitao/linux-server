# 加载必要的包
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(Seurat)
library(tidyverse)
library(future)

# ============================================================================
# 0. 加载自定义的computeCommunProb函数并确保它被使用
# ============================================================================
# 加载自定义函数
source("/home/lin/R_linux_server/linux_server/computeCommunProb_1.R")

# 创建一个包装函数，确保调用我们的自定义函数
my_computeCommunProb <- computeCommunProb

# 替换CellChat包中的函数
assignInNamespace("computeCommunProb", my_computeCommunProb, ns = "CellChat")

# ============================================================================
# 1. 数据准备
# ============================================================================
rm(list = ls())
# 检查Seurat版本
print(paste("Seurat版本:", packageVersion("Seurat")))

# 加载合并的数据集
merged_seurat <- readRDS("/home/lin/c_group/HSCCD4T_cells.rds")

# 检查数据结构
print("数据维度:")
print(dim(merged_seurat))
print("细胞类型分布:")
print(table(merged_seurat$cell_type))

# 合并Seurat对象的层(适用于Seurat v5)
merged_seurat <- JoinLayers(merged_seurat)

# 提取必要的数据组件
counts <- GetAssayData(merged_seurat, slot = "counts", layer = "counts")
data <- GetAssayData(merged_seurat, slot = "data", layer = "data")  # 标准化数据
metadata <- merged_seurat@meta.data

# 创建新的Seurat对象(可选步骤，用于验证)
new_seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)
new_seurat[["RNA"]]$data <- data

# 验证数据结构
print("新Seurat对象维度:")
print(dim(new_seurat))
print("元数据行名与数据列名匹配情况:")
print(all(rownames(new_seurat@meta.data) == colnames(new_seurat)))

# ============================================================================
# 2. 创建CellChat对象
# ============================================================================

# 准备CellChat输入数据
data_input <- as.matrix(data)
meta_data <- tibble(cell_type = merged_seurat$cell_type) %>%
  as.data.frame()
rownames(meta_data) <- colnames(data_input)

# 创建CellChat对象
cellchat <- createCellChat(object = data_input, meta = meta_data, group.by = "cell_type")
print("CellChat对象创建完成，细胞类型包括:")
print(levels(cellchat@idents))

# ============================================================================
# 3. 设置配体-受体数据库
# ============================================================================

# 设置数据库(根据物种选择human或mouse)
CellChatDB <- CellChatDB.mouse  # 如果是人类样本，使用CellChatDB.human
print("数据库类别:")
showDatabaseCategory(CellChatDB)

# 使用除"Non-protein Signaling"外的所有CellChatDB
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use

# ============================================================================
# 4. 数据预处理与通讯网络推断
# ============================================================================

# 设置并行计算参数
options(future.globals.maxSize = 300 * 1024^3)  # 设置全局变量大小限制
plan("multisession", workers = 30)  # 设置并行处理

# 数据预处理
cellchat <- subsetData(cellchat)  # 提取信号基因表达数据
cellchat <- identifyOverExpressedGenes(cellchat)  # 识别过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)  # 识别过表达的相互作用

# 输出信息，表明将使用优化版的computeCommunProb函数
cat("使用优化版的computeCommunProb函数进行计算...\n")

# 计算通讯概率并推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, type = "triMean")  # 计算通讯概率
cellchat <- filterCommunication(cellchat, min.cells = 10)  # 过滤细胞数量过少的通讯
