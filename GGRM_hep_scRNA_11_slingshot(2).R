# 加载必要的包
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(tradeSeq)
library(ggplot2)
library(dplyr)
library(ggrepel)

rm(list = ls())

#----------------------------------------------------------------
# 1. 数据准备与轨迹分析
#----------------------------------------------------------------
# 读取数据并设置输出路径
scRNAsub <- readRDS("/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")
output_dir <- "/home/lin/GGRM_hep/result/"
dir.create(output_dir, showWarnings = FALSE)

# 转换为SingleCellExperiment对象
counts_matrix <- as.matrix(scRNAsub[["RNA"]]$counts)
data_matrix <- as.matrix(scRNAsub[["RNA"]]$data)
umap_coords <- as.matrix(scRNAsub@reductions$umap@cell.embeddings)
clusters <- as.character(scRNAsub$seurat_clusters)

sce <- SingleCellExperiment(assays = list(counts = counts_matrix, logcounts = data_matrix))
colData(sce)$clusters <- clusters
reducedDims(sce) <- list(UMAP = umap_coords)

# 运行slingshot轨迹分析
sce <- slingshot(sce, clusterLabels = 'clusters', reducedDim = 'UMAP')
n_lineages <- length(slingLineages(sce))
cat("识别到", n_lineages, "条轨迹\n")

#----------------------------------------------------------------
# 2. 轨迹可视化
#----------------------------------------------------------------
# 创建获取曲线数据的函数
get_curve_data <- function(sce) {
  curves <- slingCurves(SlingshotDataSet(sce))
  curve_data <- lapply(1:length(curves), function(i) {
    tryCatch({
      points <- curves[[i]]$s[, 1:2]
      n_points <- nrow(points)
      directions <- rbind(
        points[2:n_points, ] - points[1:(n_points-1), ],
        points[n_points, ] - points[n_points-1, ]
      )
      magnitudes <- sqrt(rowSums(directions^2))
      directions <- directions / pmax(magnitudes, 1e-10)
      data.frame(
        x = points[, 1], y = points[, 2],
        dx = directions[, 1], dy = directions[, 2],
        segment = 1:n_points, curve_id = paste0("Lineage", i)
      )
    }, error = function(e) NULL)
  })
  do.call(rbind, Filter(Negate(is.null), curve_data))
}

curve_data <- get_curve_data(sce)

# 创建颜色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
lineage_colors <- brewer.pal(min(9, n_lineages), "Set1")

# 2.1 绘制所有轨迹（黑白版）
plot(umap_coords, col = rgb(0,0,0,0.5), pch=16, asp = 1,
     main = "All Lineages", xlab = "UMAP 1", ylab = "UMAP 2")
lines(SlingshotDataSet(sce), lwd=2, col=lineage_colors)
# 添加方向箭头
for(i in 1:n_lineages) {
  curve_i <- curve_data[curve_data$curve_id == paste0("Lineage", i),]
  if(nrow(curve_i) > 0) {
    arrow_indices <- unique(c(seq(10, nrow(curve_i), by = 10), nrow(curve_i)))
    for(idx in arrow_indices) {
      arrows(curve_i$x[idx] - curve_i$dx[idx]*0.5, 
             curve_i$y[idx] - curve_i$dy[idx]*0.5,
             curve_i$x[idx] + curve_i$dx[idx]*0.5, 
             curve_i$y[idx] + curve_i$dy[idx]*0.5,
             length = 0.1, col = lineage_colors[i], lwd = 2)
    }
  }
}
legend("right", legend = paste0("Lineage ", 1:n_lineages),
       col = lineage_colors, lwd = 2, cex = 0.8)

# 2.2 绘制带有伪时间颜色的细胞
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plotcol[is.na(plotcol)] <- "lightgrey"  # 保持灰色细胞
plot(umap_coords, col = plotcol, pch=16, asp = 1,
     main = "All Lineages with Pseudotime Colors",
     xlab = "UMAP 1", ylab = "UMAP 2")
lines(SlingshotDataSet(sce), lwd=2, col="black")
# 添加方向箭头
for(i in 1:n_lineages) {
  curve_i <- curve_data[curve_data$curve_id == paste0("Lineage", i),]
  if(nrow(curve_i) > 0) {
    arrow_indices <- unique(c(seq(10, nrow(curve_i), by = 10), nrow(curve_i)))
    for(idx in arrow_indices) {
      arrows(curve_i$x[idx] - curve_i$dx[idx]*0.5, 
             curve_i$y[idx] - curve_i$dy[idx]*0.5,
             curve_i$x[idx] + curve_i$dx[idx]*0.5, 
             curve_i$y[idx] + curve_i$dy[idx]*0.5,
             length = 0.1, col = "black", lwd = 2)
    }
  }
}
legend_colors <- colors[seq(1, 100, length.out = 5)]
legend_values <- round(seq(min(sce$slingPseudotime_1, na.rm=TRUE), 
                           max(sce$slingPseudotime_1, na.rm=TRUE), 
                           length.out = 5), 1)
legend("bottomright", legend = legend_values, 
       fill = legend_colors, title = "Pseudotime", cex = 0.8)

# 2.3 使用ggplot2绘制
plot_data <- data.frame(
  UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2],
  Pseudotime = sce$slingPseudotime_1,
  Cluster = clusters, InTrajectory = !is.na(sce$slingPseudotime_1)
)

p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ifelse(InTrajectory, Pseudotime, NA)), size = 1, alpha = 0.8) +
  scale_color_gradientn(colors = colors, name = "Pseudotime", na.value = "lightgrey") +
  geom_path(data = curve_data, aes(x = x, y = y, group = curve_id), color = "black", size = 1.2) +
  geom_segment(data = curve_data[seq(1, nrow(curve_data), by = 10),], 
               aes(x = x, y = y, xend = x + dx*0.5, yend = y + dy*0.5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  geom_text_repel(data = aggregate(cbind(UMAP_1, UMAP_2) ~ Cluster, plot_data, mean),
                  aes(label = Cluster), size = 5, fontface = "bold") +
  theme_minimal() + theme(panel.grid = element_blank()) +
  labs(title = "Cell Trajectory Analysis", x = "UMAP 1", y = "UMAP 2")
print(p)

#----------------------------------------------------------------
# 3. 为每条轨迹创建单独的可视化
#----------------------------------------------------------------
for(i in 1:n_lineages) {
  output_file_lineage <- paste0(output_dir, "slingshot_lineage_", i, ".png")
  png(output_file_lineage, width = 2400, height = 2000, res = 300)
  
  # 获取当前轨迹的伪时间和颜色
  pseudotime_col <- paste0("slingPseudotime_", i)
  current_pseudotime <- sce[[pseudotime_col]]
  plotcol <- rep("lightgrey", length(current_pseudotime))
  in_lineage <- !is.na(current_pseudotime)
  if(sum(in_lineage) > 0) {
    plotcol[in_lineage] <- colors[cut(current_pseudotime[in_lineage], breaks=100)]
  }
  
  # 绘制细胞和轨迹
  plot(umap_coords, col = plotcol, pch=16, asp = 1,
       main = paste("Lineage", i, "with Pseudotime"),
       xlab = "UMAP 1", ylab = "UMAP 2")
  lines(SlingshotDataSet(sce), lwd=1, col="lightgrey")
  lines(SlingshotDataSet(sce), lwd=3, col=lineage_colors[i], linInd=i)
  
  # 添加方向箭头
  curve_i <- curve_data[curve_data$curve_id == paste0("Lineage", i),]
  if(nrow(curve_i) > 0) {
    arrow_indices <- unique(c(seq(10, nrow(curve_i), by = 10), nrow(curve_i)))
    for(idx in arrow_indices) {
      arrows(curve_i$x[idx] - curve_i$dx[idx]*0.5, 
             curve_i$y[idx] - curve_i$dy[idx]*0.5,
             curve_i$x[idx] + curve_i$dx[idx]*0.5, 
             curve_i$y[idx] + curve_i$dy[idx]*0.5,
             length = 0.1, col = lineage_colors[i], lwd = 2)
    }
  }
  
  # 添加图例
  if(sum(in_lineage) > 0) {
    legend_colors <- colors[seq(1, 100, length.out = 5)]
    legend_values <- round(seq(min(current_pseudotime, na.rm=TRUE), 
                               max(current_pseudotime, na.rm=TRUE), 
                               length.out = 5), 1)
    legend("bottomright", legend = legend_values, 
           fill = legend_colors, title = "Pseudotime", cex = 0.8)
  }
  dev.off()
}

#----------------------------------------------------------------
# 4. 差异基因分析 - 针对特定轨迹
#----------------------------------------------------------------
# 选择要分析的轨迹
selected_lineage <- 1  # 可以改为1, 2, 或3

# 获取所选轨迹的伪时间和细胞
pseudotime_col <- paste0("slingPseudotime_", selected_lineage)
selected_pseudotime <- sce[[pseudotime_col]]
in_lineage <- !is.na(selected_pseudotime)
cat("轨迹", selected_lineage, "包含", sum(in_lineage), "个细胞\n")

# 过滤基因并准备数据
geneFilter <- apply(counts_matrix, 1, function(x) sum(x >= 3) >= ncol(counts_matrix) * 0.05)
filtered_counts <- counts_matrix[geneFilter, in_lineage]
lineage_pseudotime <- matrix(selected_pseudotime[in_lineage], ncol=1)
lineage_weights <- matrix(1, nrow=sum(in_lineage), ncol=1)

# 评估最佳knots数量
set.seed(123)
if(sum(geneFilter) > 500) {
  sample_genes <- sample(which(geneFilter), 500)
  eval_counts <- counts_matrix[sample_genes, in_lineage]
} else {
  eval_counts <- filtered_counts
}

cat("开始评估轨迹", selected_lineage, "的最佳knots数量...\n")
icMat <- evaluateK(counts = as.matrix(eval_counts), 
                   pseudotime = lineage_pseudotime,
                   cellWeights = lineage_weights,
                   k = 3:8, nGenes = 500, verbose = TRUE)

# 拟合GAM模型
cat("开始为轨迹", selected_lineage, "拟合GAM模型...\n")
sce_gam <- fitGAM(counts = as.matrix(filtered_counts),
                  pseudotime = lineage_pseudotime,
                  cellWeights = lineage_weights,
                  nknots = 6)

# 评估模型收敛情况
convergence_rate <- mean(rowData(sce_gam)$tradeSeq$converged, na.rm = TRUE)
cat("轨迹", selected_lineage, "的模型收敛率:", convergence_rate, "\n")

# 测试基因与伪时间的关联
ATres <- associationTest(sce_gam)
write.csv(ATres, file = paste0(output_dir, "lineage_", selected_lineage, "_association_test_results.csv"))

#----------------------------------------------------------------
# 5. 可视化差异基因
#----------------------------------------------------------------
# 热图可视化
top50 <- rownames(ATres[order(ATres$pvalue), ])[1:min(50, nrow(ATres))]
pst.ord <- order(selected_pseudotime[in_lineage])
heatdata <- log1p(as.matrix(counts_matrix[top50, in_lineage][, pst.ord]))
heatclus <- clusters[in_lineage][pst.ord]

output_file4 <- paste0(output_dir, "lineage_", selected_lineage, "_gene_expression_heatmap.png")
png(output_file4, width = 2500, height = 3000, res = 300)
heatmap(heatdata, Colv = NA,
        ColSideColors = brewer.pal(min(9, length(unique(heatclus))), "Set1")[as.numeric(factor(heatclus))],
        labCol = "", main = paste("Lineage", selected_lineage, "Gene Expression along Pseudotime"))
dev.off()

# 单个基因的表达曲线
top_genes <- rownames(ATres[order(ATres$pvalue), ])[1:5]
for(i in 1:length(top_genes)) {
  gene <- top_genes[i]
  
  # 准备数据
  gene_expr <- counts_matrix[gene, in_lineage]
  plot_data <- data.frame(
    pseudotime = selected_pseudotime[in_lineage],
    expression = gene_expr
  )
  
  # 绘制散点图和平滑曲线
  output_file_gene <- paste0(output_dir, "lineage_", selected_lineage, "_gene_", gene, "_dynamics.png")
  png(output_file_gene, width = 2000, height = 1500, res = 300)
  plot(plot_data$pseudotime, plot_data$expression,
       main = paste("Lineage", selected_lineage, "- Expression of", gene, "along Pseudotime"),
       xlab = "Pseudotime", ylab = "Expression",
       pch = 16, col = rgb(0,0,0,0.3))
  
  # 添加平滑曲线
  smooth_fit <- loess(expression ~ pseudotime, data = plot_data, span = 0.3)
  pseudotime_order <- order(plot_data$pseudotime)
  lines(plot_data$pseudotime[pseudotime_order], 
        predict(smooth_fit)[pseudotime_order], 
        col = "red", lwd = 3)
  dev.off()
}

cat("轨迹分析完成，所有图像已保存到:", output_dir, "\n")
