library(scTernary)
library(Seurat)
library(readxl)
rm(list = ls())
# 读取数据
seu_obj <- readRDS("/home/lin/c_group/CD8+T_ident.rds")
gene_list <- read_excel("/home/lin/c_group/san.xlsx")

# 设置样本顺序
sample_order <- c("Con", "Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R")

# 提取基因集
gene_sets <- list(
  Cytotoxicity = na.omit(gene_list$Cytotoxicity),
  Exhaustion = na.omit(gene_list$Exhaustion),
  Activation = na.omit(gene_list$`Activation:Effector function`)
)

# 创建注释基因数据框
anno_signature_genes <- data.frame(
  Symbol = unlist(gene_sets),
  gene_type = rep(names(gene_sets), lengths(gene_sets))
)

# 获取表达矩阵
if("counts" %in% names(seu_obj@assays$RNA@layers)) {
  expr_mat <- as.matrix(seu_obj@assays$RNA@layers$counts)
} else {
  expr_mat <- as.matrix(seu_obj@assays$RNA@data)
}
rownames(expr_mat) <- rownames(seu_obj@assays$RNA)
colnames(expr_mat) <- colnames(seu_obj@assays$RNA)

# 生成三元图数据
data_for_ternary <- generate_data_for_ternary(
  data_exp_mat = expr_mat,
  anno_signature_genes = anno_signature_genes,
  gene_name_col = "Symbol",
  gene_type_col = "gene_type",
  weight_by_gene_count = TRUE
)

# 转换为数据框并添加分组
ternary_df <- as.data.frame(data_for_ternary)
# 关键步骤：将group设置为因子并指定顺序
ternary_df$group <- factor(seu_obj$orig.ident, levels = sample_order)

# 绘制三元图
pdf("CD8T_ternary_plot.pdf", width = 8, height = 7)
vcdTernaryPlot(ternary_df,
               order_colnames = 1:3,
               point_size = 0.6,
               group = ternary_df$group,
               show_legend = TRUE)
dev.off()

# 使用 ggtern 包，它提供更多自定义选项
library(ggtern)

# 准备数据
plot_data <- data.frame(
  Cytotoxicity = ternary_df[,1],
  Exhaustion = ternary_df[,2],
  Activation = ternary_df[,3],
  group = ternary_df$group
)

# 绘制三元图
p <- ggtern(plot_data, aes(x = Cytotoxicity, y = Exhaustion, z = Activation, color = group)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = rainbow(length(levels(plot_data$group)))) +
  theme_bw() +
  theme(legend.position = "right") +
  labs(title = "CD8+ T Cell Functional States")
print(p)
dev.off()
ggsave("CD8T_ternary_ggtern.pdf", p, width = 10, height = 8)
