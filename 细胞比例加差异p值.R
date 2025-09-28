rm(list=ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(dunn.test)
library(rstatix)
library(FSA)  # for dunnTest function

# 设置工作目录
setwd("/home/lin/B_group/")
filepath <- "/home/lin/B_group/Hepatocytes.RDS"
pbmc <- readRDS(filepath)

# 定义颜色调色板
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

# 设置细胞类型和分组信息
pbmc$celltype <- pbmc$seurat_clusters
group_c = c("Bd30R","Bd15R","Bd8R","Bd1R","Ad1R")

# 创建颜色调色板
colors = brewer.pal(14, "Set3")
colors <- c(colors, "black", "white")

# 获取聚类ID并设置主题
cluster_id <- levels(pbmc@meta.data$seurat_clusters)
theme_set(theme_cowplot())
result <- setNames(color[1:length(cluster_id)], cluster_id)

# 排序细胞类型并子集数据
pbmc@meta.data$celltype <- ordered(pbmc@meta.data$celltype, 
                                   levels = cluster_id)
cell <- subset(pbmc, subset = celltype %in% cluster_id)
cell <- ScaleData(cell)

# 提取细胞类型和分组信息
cell_counts <- FetchData(cell, vars = c("celltype","group","orig.ident")) %>%
  mutate(group = factor(group, levels = group_c))

# 按细胞类型和组计数
cell_counts_tbl <- cell_counts %>%
  dplyr::count(celltype, group)

# 计算每个组中的比例
cell_proportions <- cell_counts %>%
  group_by(group) %>%
  mutate(total = n()) %>%
  group_by(group, celltype) %>%
  summarise(count = n(), 
            proportion = count/first(total),
            .groups = "drop")

# 保存数据
write.csv(cell_counts_tbl, file = "/home/lin/B_group/cell_counts_tbl.csv")
write.csv(cell_counts, file = "/home/lin/B_group/cell_counts.csv")
write.csv(cell_proportions, file = "/home/lin/B_group/cell_proportions.csv")

# 为统计分析创建"伪复制"样本
set.seed(123) # 设置随机种子以确保结果可重复
n_replicates <- 5 # 每组创建5个"伪复制"

# 为每个细胞分配一个伪复制ID
cell_counts$replicate <- sample(1:n_replicates, nrow(cell_counts), replace = TRUE)

# 计算每个伪复制中的细胞类型比例
pseudo_proportions <- cell_counts %>%
  group_by(group, replicate) %>%
  mutate(total = n()) %>%
  group_by(group, replicate, celltype) %>%
  summarise(count = n(), 
            proportion = count/first(total),
            .groups = "drop")

# 统计分析结果存储
stat_results <- data.frame(
  celltype = character(),
  p_value = numeric(),
  significance = character(),
  stringsAsFactors = FALSE
)

# 对每个细胞类型进行Kruskal-Wallis检验
for (ct in unique(cell_counts$celltype)) {
  # 提取当前细胞类型的数据
  ct_data <- pseudo_proportions %>% 
    filter(celltype == ct)
  
  # Kruskal-Wallis检验
  kw_result <- kruskal.test(proportion ~ group, data = ct_data)
  p_val <- kw_result$p.value
  
  # 添加到结果
  significance <- ""
  if (p_val < 0.001) significance <- "***"
  else if (p_val < 0.01) significance <- "**"
  else if (p_val < 0.05) significance <- "*"
  
  stat_results <- rbind(stat_results, data.frame(
    celltype = ct,
    p_value = p_val,
    significance = significance
  ))
  
  # 如果显著，进行Dunn's事后检验
  if (p_val < 0.05) {
    dunn_result <- dunnTest(proportion ~ group, data = ct_data, method = "bonferroni")
    # 保存Dunn's检验结果
    write.csv(dunn_result$res, 
              file = paste0("/Users/joe/Desktop/B组细胞数据/dunn_test_", 
                            gsub(" ", "_", ct), ".csv"),
              row.names = FALSE)
  }
}

# 创建更详细的统计摘要表
stat_summary <- stat_results %>%
  arrange(p_value) %>%
  mutate(
    p_value_formatted = sprintf("%.4f", p_value),
    significant = ifelse(p_value < 0.05, "Yes", "No")
  ) %>%
  select(celltype, p_value_formatted, significance, significant)

colnames(stat_summary) <- c("Cell Type", "P-value", "Significance", "Significantly Different")

# 保存统计摘要表
write.csv(stat_summary, file = "/home/lin/B_group/cell_type_statistical_summary.csv", 
          row.names = FALSE)

# 修改cluster_id以包含显著性标记
cluster_id_with_sig <- cluster_id
for (i in 1:length(cluster_id)) {
  ct <- cluster_id[i]
  sig_row <- stat_results[stat_results$celltype == ct, ]
  if (nrow(sig_row) > 0 && sig_row$significance != "") {
    cluster_id_with_sig[i] <- paste0(cluster_id[i], " ", sig_row$significance)
  }
}

# 创建带有更新标签的新颜色映射
result_with_sig <- setNames(color[1:length(cluster_id)], cluster_id_with_sig)

# 创建带有图例中显著性的图
p <- ggplot(data = cell_counts, aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = result, labels = cluster_id_with_sig) +
  coord_flip() +
  scale_y_reverse() +
  theme(legend.text = element_text(size = 20),
        text = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  labs(y = "Proportion", x = "Group", 
       fill = "Cell Type\n* p<0.05, ** p<0.01, *** p<0.001")

# 保存图
png("/Users/joe/Desktop/B组细胞数据/B组细胞比例_with_significance_legend.png", 
    units = "in", width = 20, height = 8, res = 600)
print(p)
dev.off()

# 创建统计结果的表格可视化
# 首先，按显著性排序
stat_table <- stat_summary %>%
  arrange(desc(Significance), `Cell Type`)

# 创建表格图
table_plot <- tableGrob(stat_table, rows = NULL, theme = ttheme_minimal(
  core = list(fg_params=list(fontsize=14),
              bg_params=list(fill=c("white", "grey95"))),
  colhead = list(fg_params=list(fontsize=16, fontface="bold"))))

# 将表格保存为图像
png("/Users/joe/Desktop/B组细胞数据/cell_type_statistical_table.png", 
    units = "in", width = 10, height = length(unique(cell_counts$celltype))*0.4 + 1, res = 300)
grid.draw(table_plot)
dev.off()

# 创建包含图和表的组合可视化
png("/Users/joe/Desktop/B组细胞数据/B组细胞比例_with_table.png", 
    units = "in", width = 24, height = 12, res = 300)
grid.arrange(p, table_plot, ncol=2, widths=c(2,1))
dev.off()