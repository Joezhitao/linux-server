# ==============================================================================
# 1. 环境准备与数据加载
# ==============================================================================
rm(list = ls())
library(Seurat)
library(scMetabolism)
devtools::install_github("YosefLab/VISION@v2.1.0")
library(dplyr)
library(tidyr)
library(networkD3)
library(viridis)
library(ggplot2)

# 加载数据 (使用你提供的路径)
seu_obj <- readRDS("/home/lin/c_group/hep_ident_filtered.rds")

# 重新设置你的细胞身份 (根据你提供的代码)
# 注意：这里假设原来的levels顺序是0到10，共11个。
# 如果你的seurat_clusters本身顺序不同，请确保这里的对应关系是正确的。
new_names <- c("Stressed-Anabolic-Hep","Regenerative-Unit-Hep", "Core-Metabolic-Hep",
               "Regenerative-Unit-Hep", "Regenerative-Unit-Hep", "Core-Metabolic-Hep",
               "Core-Metabolic-Hep", "Stressed-Anabolic-Hep", "Stressed-Anabolic-Hep", 
               "Core-Metabolic-Hep", "Regenerative-Unit-Hep")
levels(seu_obj@meta.data$seurat_clusters) <- new_names
Idents(seu_obj) <- "seurat_clusters" # 确保Idents是你的新名字

# ==============================================================================
# 2. 执行代谢分析 (scMetabolism)
# ==============================================================================
# 定义代谢分析核心函数 (精简版，仅保留核心功能)
run_scMetabolism <- function(obj, metabolism.type = "KEGG", ncores = 4) {
  countexp <- as.matrix(GetAssayData(obj, layer = "counts"))
  
  # 获取代谢基因集文件路径
  gmtFile <- system.file("data", paste0(metabolism.type, "_metabolism_nc.gmt"), package = "scMetabolism")
  if(gmtFile == "") { # 如果找不到，尝试Reactome或报错提示
    gmtFile <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")
  }
  
  cat("正在计算代谢评分 (使用 AUCell 方法)...\n")
  library(AUCell)
  library(GSEABase)
  
  # 1. 构建排名
  cells_rankings <- AUCell_buildRankings(countexp, nCores = ncores, plotStats = F)
  # 2. 读取基因集
  geneSets <- getGmt(gmtFile)
  # 3. 计算AUC
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # 将结果存入Seurat对象
  signature_exp <- data.frame(getAUC(cells_AUC))
  obj[["METABOLISM"]] <- CreateAssayObject(data = as.matrix(signature_exp))
  
  return(obj)
}

# 运行分析
seu_obj <- run_scMetabolism(seu_obj, metabolism.type = "KEGG", ncores = 4)

# ==============================================================================
# 3. 准备桑葚图数据
# ==============================================================================
# 提取代谢评分矩阵
metabolism.matrix <- GetAssayData(seu_obj, assay = "METABOLISM", layer = "data")

# --- 策略A: 自动筛选差异最大的通路 (推荐) ---
# 计算每个亚群的平均代谢活性
avg_metabolism <- AverageExpression(seu_obj, assays = "METABOLISM")$METABOLISM

# 计算通路的变异程度 (标准差)，取变异最大的前 20 条通路
pathway_sd <- apply(avg_metabolism, 1, sd)
top_pathways <- names(sort(pathway_sd, decreasing = TRUE))[1:20]

# --- 策略B: 如果你想指定特定通路，请取消下面注释并修改 ---
# top_pathways <- c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", 
#                   "Fatty acid metabolism", "Oxidative phosphorylation")

cat("选定的绘图通路:\n")
print(top_pathways)

# 构建绘图数据框
# 我们要展示：通路 (Source) -> 细胞亚群 (Target)
plot_data <- as.data.frame(avg_metabolism[top_pathways, ])
plot_data$Pathway <- rownames(plot_data)

# 转换为长格式 (Long format)
links <- plot_data %>%
  pivot_longer(cols = -Pathway, names_to = "CellType", values_to = "Score") %>%
  filter(Score > quantile(Score, 0.25)) %>% # 过滤掉过低的分数，让图更干净
  mutate(Value = Score * 100) # 放大数值以便绘图

# 整理 Source 和 Target
# networkD3 需要 Source 和 Target 都是从 0 开始的索引
nodes <- data.frame(name = unique(c(links$Pathway, links$CellType)))

# 区分节点组别 (用于上色)
nodes$group <- ifelse(nodes$name %in% links$Pathway, "Pathway", "CellType")

# 将名称转换为索引 ID
links$IDsource <- match(links$Pathway, nodes$name) - 1
links$IDtarget <- match(links$CellType, nodes$name) - 1

# ==============================================================================
# 4. 绘制并保存桑葚图
# ==============================================================================
# 定义颜色方案 (Javascript 语法)
# 这里的颜色代码可以根据喜好修改
my_color <- 'd3.scaleOrdinal() .domain(["Pathway", "CellType"]) .range(["#69b3a2", "#404080"])'

sankey <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "Value",
  NodeID = "name",
  NodeGroup = "group",      # 节点按组别上色
  LinkGroup = "Pathway",    # 连线跟随通路颜色 (或者改为 "CellType" 跟随细胞)
  units = "Score",
  fontSize = 14,
  nodeWidth = 30,
  nodePadding = 10,         # 节点间距
  sinksRight = FALSE,       # 设置为 FALSE 防止标签过于靠边
  height = 600,
  width = 900
)

# 显示图表
sankey

# 保存为 HTML 文件
saveNetwork(sankey, file = "Hep_Metabolism_Sankey.html")

cat("完成！桑葚图已保存为 'Hep_Metabolism_Sankey.html'\n")
