# 清空环境并加载必要的包
rm(list = ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(xlsx)

# 设置工作目录和读取数据
setwd("/home/lin/c_group/")

pbmc <- readRDS("/home/lin/c_group/CD8+T_ident.rds")
pbmc@meta.data$orig.ident
group <- levels(pbmc@meta.data$orig.ident)

# 读取marker基因数据框
g_df <- read.xlsx("/home/lin/c_group/gene.xlsx",
                  sheetIndex = 1,
                  header = T,
                  encoding = "UTF-8")
# 获取基因列表
gene <- g_df$gene
gene <- gene[gene %in% rownames(pbmc)]

#从Seurat对象中提取细胞注释以及基因表达量
vln.dat=FetchData(pbmc,c(gene,"orig.ident"))

# 定义因子顺序，防止画图时对细胞注释进行重排
vln.dat$orig.ident=factor(vln.dat$orig.ident,
                              levels = c("Ad1R", 
                                         "Cd1L", 
                                         "Cd1R", 
                                         "Cd15L", 
                                         "Cd15R", 
                                         "Cd30L", 
                                         "Cd30R"))
vln.dat=vln.dat[order(vln.dat$orig.ident),]
vln.dat$orig.ident=factor(vln.dat$orig.ident,levels = unique(vln.dat$orig.ident))
vln.dat.melt=vln.dat %>% 
  reshape2::melt(,gene) %>%
  rename("Gene"="variable") %>%
  group_by(orig.ident,Gene) %>%
  mutate(fillcolor=mean(value))

# 小提琴图的填充颜色
pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)

# 堆积小提琴图
p1 = ggplot(vln.dat.melt,aes(x=orig.ident,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(Gene~.)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left"
  )
p1
