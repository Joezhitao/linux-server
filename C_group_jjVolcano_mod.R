# 清除环境和函数缓存
rm(list = ls())
# 重新加载必要的库
library(readxl)
library(ggplot2)
library(corrplot)
library(dplyr)
library(ggrepel)
library(geomtextpath)
library(jjAnno)

# 重新加载修改后的函数
source("/home/lin/R_linux_server/linux_server/jjVolcano_mod(魔改火山图).R")

# 读取数据
data <- read_excel("/home/lin/c_group/subpopulation_analysis/Hepatocytes/DiffGenes_LR_comparison/All_LR_comparisons.xlsx")

# 创建图形
p <- jjVolcano_mod(diffData = data,
                   tile.col = corrplot::COL2('RdYlBu', 15)[4:12],
                   size = 3.5,
                   topGeneN = 10,
                   fontface = 'italic',
                   polar = T,
                   label_box = TRUE,  # 启用带框的标签
                   nudge_y = 2,  # 增加标签向外的推力
                   segment.size = 0.3,  # 减小连线宽度
                   min.segment.length = 0) +  # 允许短线段
  ylim(-8,12)  # 增大y轴范围，给标签留出更多空间

# 显示图形
print(p)

# 保存图形
ggsave("volcano_plot.pdf", p, width = 10, height = 8)
