# ==============================================================================
# 1. 环境准备与数据加载
# ==============================================================================
rm(list = ls())

# 加载必要的包
library(Seurat)
library(tidyverse) # 包含了 ggplot2, dplyr, tidyr, stringr 等
library(patchwork)
library(xlsx)
library(schard)
library(RColorBrewer) # 显式加载颜色包

# 设置工作目录 (根据实际情况调整)
setwd("/home/lin/gfpdata")

# 读取h5ad文件并转换为Seurat对象
# 注意：确保文件路径正确
seu_obj <- schard::h5ad2seurat('/home/lin/gfpdata/CD8 T cells_processed.h5ad')

# 修正降维信息 (如果转换后名字是 Xumap_)
if("Xumap_" %in% names(seu_obj@reductions)){
  seu_obj@reductions$umap <- seu_obj@reductions$Xumap_
}

# 确保默认Assay是RNA或其他标准化后的Assay，防止提取到raw counts
DefaultAssay(seu_obj) <- "RNA" 

# ==============================================================================
# 2. 数据提取与清洗
# ==============================================================================
# 获取基因列表 (支持多个基因，例如 c("Pdcd1", "Cd8a"))
gene <- c("Pdcd1", "Layn", "Havcr2", "Lag3", "Cd244", "Ctla4", 
          "Lilrb1", "Tigit", "Tox", "Vsir", "Btla", "Entpd1", 
          "Cd160", "Lair1")

# 检查基因是否存在于对象中，防止报错
valid_genes <- gene[gene %in% rownames(seu_obj)]
if(length(valid_genes) == 0) stop("指定的基因在Seurat对象中不存在")

# 从Seurat对象中提取细胞注释以及基因表达量
vln.dat <- FetchData(seu_obj, vars = c(valid_genes, "orig.ident"))

# 定义因子顺序 (关键步骤：确保画图顺序符合逻辑)
# 这一步只需做一次即可
target_levels <- c("Con", "HT", "RT", "RT+HT")
vln.dat$orig.ident <- factor(vln.dat$orig.ident, levels = target_levels)

# 数据长宽转换 (使用 pivot_longer 代替 melt，语法更清晰)
vln.dat.melt <- vln.dat %>% 
  # 将基因列转换为长数据
  pivot_longer(
    cols = all_of(valid_genes), # 指定要转换的基因列
    names_to = "Gene",          # 新的列名用来存放基因名
    values_to = "value"         # 新的列名用来存放表达量
  ) %>%
  # 计算每组的均值用于填充颜色
  group_by(orig.ident, Gene) %>%
  mutate(fillcolor = mean(value)) %>%
  ungroup() # 解除分组，是个好习惯

# ==============================================================================
# 3. 绘图 (Stacked Violin Plot)
# ==============================================================================

# 定义颜色板 (YlOrRd 渐变色)
pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))(100)

# 绘图
p1 <- ggplot(vln.dat.melt, aes(x = orig.ident, y = value, fill = fillcolor)) +
  # 小提琴图：scale="width" 保证不同组宽度一致，linetype="blank" 去除边框
  geom_violin(scale = "width", trim = TRUE, linetype = "blank", alpha = 0.9) +
  
  # 颜色映射
  scale_fill_gradientn(colors = pal, name = "Avg Expression") +
  
  # 分面展示：每个基因一行
  facet_grid(Gene ~ ., scales = "free_y", switch = "y") +
  
  # SCI 风格主题设置
  theme_classic() +
  theme(
    # 去除网格线和背景
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    
    # x轴设置
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12),
    axis.title.x = element_blank(),
    
    # y轴设置 (通常这种堆积图不需要显示具体的y轴数值，或者只显示简化的)
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.title.y = element_blank(),
    
    # 分面标签位置 (基因名在左侧或右侧)
    strip.text.y = element_text(angle = 0, size = 12, face = "bold", hjust = 0),
    strip.placement = "outside", # 配合 switch="y" 让基因名在轴外侧
    
    # 图例设置
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

# 显示图片
print(p1)

# 保存图片 (可选)
# ggsave("Stacked_Violin_Pdcd1.pdf", p1, width = 6, height = 4)
