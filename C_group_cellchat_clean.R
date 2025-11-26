# ============================================================================
# 肝脏细胞通讯分析 - 使用CellChat
# ============================================================================

# 加载必要的包
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(Seurat)
library(tidyverse)
library(future)

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

# 计算通讯概率并推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, type = "triMean")  # 计算通讯概率
cellchat <- filterCommunication(cellchat, min.cells = 10)  # 过滤细胞数量过少的通讯

# 在信号通路水平上推断细胞通讯
cellchat <- computeCommunProbPathway(cellchat)  # 计算通路水平通讯概率
cellchat <- aggregateNet(cellchat)  # 计算聚合的细胞通讯网络

# 保存CellChat对象
saveRDS(cellchat, file = "/home/lin/c_group/cellchat/cellchat_HSCCD4T.rds")

# ============================================================================
# 5. 可视化细胞通讯网络
# ============================================================================

# 获取所有显著通讯的通路
pathways.show.all <- cellchat@netP$pathways
print(paste("识别到的显著通路数量:", length(pathways.show.all)))
print("前5个显著通路:")
print(pathways.show.all[1:min(5, length(pathways.show.all))])

# 获取细胞类型信息
cell_types <- levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

# 5.1 基本网络可视化
pdf("/home/lin/c_group/cellchat/network_visualization.pdf", width = 10, height = 8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge= FALSE, title.name = "互作数量")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge= FALSE, title.name = "互作强度")
dev.off()

# 5.2 气泡图可视化
pdf("/home/lin/c_group/cellchat/bubble_plots.pdf", width = 12, height = 10)

# 创建细胞类型索引的命名列表，方便后续引用
cell_indices <- setNames(1:length(cell_types), cell_types)

# 为每个细胞类型对绘制气泡图
cell_pairs <- expand.grid(source = cell_types, target = cell_types) %>%
  filter(source != target) %>%
  as_tibble()

for (i in 1:nrow(cell_pairs)) {
  source_name <- cell_pairs$source[i]
  target_name <- cell_pairs$target[i]
  source_idx <- cell_indices[source_name]
  target_idx <- cell_indices[target_name]
  
  p <- netVisual_bubble(cellchat, sources.use = source_idx, targets.use = target_idx, 
                        remove.isolate = FALSE)
  print(p + ggtitle(paste0(source_name, " 到 ", target_name)))
}

# 总体互作气泡图
p <- netVisual_bubble(cellchat, remove.isolate = FALSE)

pdf("/home/lin/c_group/cellchat/bubble_plots.pdf", width = 12, height = 25)
print(p + ggtitle("所有细胞类型间互作"))
dev.off()

# 如果有显著通路，为第一个通路绘制气泡图
if(length(pathways.show.all) > 0) {
  pathway_name <- pathways.show.all[1]
  p <- netVisual_bubble(cellchat, sources.use = 1:length(cell_types), 
                        targets.use = 1:length(cell_types),
                        signaling = pathway_name, remove.isolate = FALSE)
  print(p + ggtitle(paste0(pathway_name, " 信号通路")))
}

dev.off()

# 5.3 和弦图可视化
pdf("/home/lin/c_group/cellchat/chord_diagrams.pdf", width = 10, height = 10)

# 为前3个通路绘制和弦图(如果存在)
for (i in 1:min(3, length(pathways.show.all))) {
  pathway_name <- pathways.show.all[i]
  netVisual_aggregate(cellchat, signaling = pathway_name, layout = "chord")
  title(main = paste0(pathway_name, " 信号通路"))
}

# 如果有显著通路，为第一个通路绘制基因级和弦图
if(length(pathways.show.all) > 0) {
  pathway_name <- pathways.show.all[1]
  pairLR.use <- extractEnrichedLR(cellchat, signaling = pathway_name, geneLR.return = FALSE)
  
  if(nrow(pairLR.use) > 0) {
    netVisual_chord_gene(cellchat, sources.use = 1:2, targets.use = 2:3, 
                         signaling = pathway_name, legend.pos.x = 8)
    title(main = paste0(pathway_name, " 配体-受体对"))
  }
}

dev.off()

# 5.4 热图可视化
pdf("/home/lin/c_group/cellchat/heatmap_pathways.pdf", width = 10, height = 8)
# 为前5个通路绘制热图(如果存在)
if(length(pathways.show.all) >= 1) {
  pathways_to_show <- pathways.show.all[1:min(5, length(pathways.show.all))]
  netVisual_heatmap(cellchat, signaling = pathways_to_show, color.heatmap = "Reds")
}
dev.off()

# 5.5 基因表达可视化
pdf("/home/lin/c_group/cellchat/gene_expression_plots.pdf", width = 10, height = 8)

# 如果有显著通路，检查第一个通路中的基因表达
if(length(pathways.show.all) > 0) {
  pathway_name <- pathways.show.all[1]
  plotGeneExpression(cellchat, signaling = pathway_name)
  
  # 显示特定配体-受体对的表达
  pairLR.use <- extractEnrichedLR(cellchat, signaling = pathway_name, geneLR.return = FALSE)
  if(nrow(pairLR.use) > 0) {
    plotGeneExpression(cellchat, features = c(pairLR.use[1,]))
  }
}

dev.off()

# ============================================================================
# 6. 系统级分析
# ============================================================================

pdf("/home/lin/c_group/cellchat/systems_analysis.pdf", width = 10, height = 8)

# 6.1 计算网络中心性得分
cellchat <- netAnalysis_computeCentrality(cellchat)
netAnalysis_signalingRole_network(cellchat, width = 8, height = 6, font.size = 10)

# 6.2 识别通讯模式
# 识别发送细胞的传出通讯模式
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

# 识别接收细胞的传入通讯模式
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming")
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

dev.off()

# 打印完成信息
print("CellChat分析完成！所有结果已保存至/home/lin/c_group/cellchat/目录")
