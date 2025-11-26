library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(Seurat)

rm(list = ls())

# Load the merged dataset
merged_seurat <- readRDS("/home/lin/c_group/merged_liver_cells.rds")

# Check the structure of the merged object
print(dim(merged_seurat))
print(table(merged_seurat$cell_type))

# Extract the necessary components - adjusted for Seurat v5
merged_seurat = JoinLayers(merged_seurat)
counts <- GetAssayData(merged_seurat, slot = "counts", layer = "counts")
data <- GetAssayData(merged_seurat, slot = "data", layer = "data")  # Normalized data
metadata <- merged_seurat@meta.data

# Create a new Seurat object with just the essential components
new_seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)

# For Seurat v5, use LayerData to add normalized data
new_seurat[["RNA"]]$data <- data

# Verify the structure
print(dim(new_seurat))
print(all(rownames(new_seurat@meta.data) == colnames(new_seurat)))

# Alternative approach: extract data directly for CellChat
data_input <- as.matrix(data)
meta_data <- data.frame(cell_type = merged_seurat$cell_type)
rownames(meta_data) <- colnames(data_input)

# Create CellChat object directly from matrices
cellchat <- createCellChat(object = data_input, meta = meta_data, group.by = "cell_type")

# Set the database
CellChatDB <- CellChatDB.mouse
# Show the structure of the database
showDatabaseCategory(CellChatDB)

# Use all CellChatDB except for "Non-protein Signaling"
CellChatDB.use <- subsetDB(CellChatDB)
# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing

# 增加全局变量大小限制到 2GB
options(future.globals.maxSize = 150 * 1024^3)  # 设置为 2GB

# 然后再运行并行处理
future::plan("multisession", workers = 4)

cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes
future::plan("multisession", workers = 4) # Do parallel processing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
# Filter out the cell-cell communication if there are only few cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Save the CellChat object
saveRDS(cellchat, file = "/home/lin/c_group/cellchat/cellchat_liver_cells.rds")

# Visualization of the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
pdf("/home/lin/c_group/cellchat/network_visualization.pdf", width = 10, height = 8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
# 绘制气泡图
pdf("/home/lin/c_group/cellchat/bubble_plots.pdf", width = 10, height = 8)

# (1) 显示从某些细胞组（由'sources.use'定义）到其他细胞组（由'targets.use'定义）的所有显著相互作用（L-R对）
# 假设Hepatocyte是第1个细胞类型，HSC是第2个，LESC是第3个
# 显示从Hepatocyte到HSC和LESC的相互作用
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2,3), remove.isolate = FALSE)

# 显示从HSC到Hepatocyte和LESC的相互作用
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3), remove.isolate = FALSE)

# 显示从LESC到Hepatocyte和HSC的相互作用
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2), remove.isolate = FALSE)

# 总图
p <- netVisual_bubble(cellchat, remove.isolate = FALSE)
pdf("/home/lin/c_group/cellchat/bubble_plots.pdf", width = 10, height = 18)
print(p)
dev.off()
# (2) 显示从某些细胞组到其他细胞组的所有显著信号通路
# 使用前面示例中的相同设置
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2,3), signaling = pathways.show.all, remove.isolate = FALSE)

# (3) 显示特定信号通路从某些细胞组到其他细胞组的所有L-R对
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all[1], geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(2,3), pairLR.use = pairLR.use, remove.isolate = FALSE)

dev.off()

# (B) 和弦图
# 和弦图用于可视化由多个配体-受体或信号通路介导的细胞-细胞通讯
pdf("/home/lin/c_group/cellchat/chord_diagrams.pdf", width = 10, height = 10)

# 使用netVisual_chord_gene显示多个配体-受体对
# 从某些细胞组到其他细胞组的所有显著配体-受体对
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(2,3), slot.name = "net", legend.pos.x = 8)

# 显示特定信号通路的所有显著配体-受体对
netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = c(2,3), signaling = pathways.show.all[1], legend.pos.x = 8)

# 显示用户提供的配体-受体对
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all[1], geneLR.return = FALSE)
netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = c(2,3), pairLR.use = pairLR.use, legend.pos.x = 8)

dev.off()

# 使用小提琴图/点图绘制信号基因表达分布
pdf("/home/lin/c_group/cellchat/gene_expression_plots.pdf", width = 10, height = 8)

# 检查信号通路中的基因表达分布
plotGeneExpression(cellchat, signaling = pathways.show.all[1])

# 显示特定配体-受体对的表达
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all[1], geneLR.return = FALSE)
plotGeneExpression(cellchat, features = c(pairLR.use[1,]))

dev.off()

