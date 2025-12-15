# 加载必要的包
library(Seurat)
library(dplyr)

# 1. 读取数据 (假设你已经读入，如果没有请取消注释)
# CD4T <- readRDS("/home/lin/c_group/CD4+T_ident.rds")
# CD8T <- readRDS("/home/lin/c_group/CD8+T_ident.rds")

# =========================================================
# 第一步：从 CD4T 中提取 Treg 细胞
# =========================================================

# 确认 Treg 对应的标签名称
# 根据你的描述，CD4T@meta.data$seurat_clusters 中包含 "CD4+ Regulatory T cells"
treg_label <- "CD4+ Regulatory T cells"

# 提取 Treg 细胞子集
# 方法：使用 subset 函数，条件是 seurat_clusters 等于 treg_label
CD4_Treg <- subset(CD4T, subset = seurat_clusters == treg_label)

# 给提取出来的 Treg 增加一个临时标签，方便合并后识别
CD4_Treg$temp_label <- "CD4+ Treg"

# =========================================================
# 第二步：准备 CD8T 数据
# =========================================================

# 给 CD8T 数据也增加一个临时标签
CD8T$temp_label <- "CD8+ T"

# =========================================================
# 第三步：合并数据
# =========================================================

# 使用 merge 函数合并两个 Seurat 对象
# add.cell.ids 可选，用于防止细胞条形码(barcode)重复冲突
seu_obj_combined <- merge(x = CD4_Treg, y = CD8T, add.cell.ids = c("Treg", "CD8"))

# =========================================================
# 第四步：设置最终的 celltype 标签
# =========================================================

# 将我们在第二、三步创建的 temp_label 赋值给 celltype
seu_obj_combined$celltype <- seu_obj_combined$temp_label

# 删除临时标签列 (可选，为了保持 metadata 整洁)
seu_obj_combined@meta.data$temp_label <- NULL

# =========================================================
# 第五步：验证结果
# =========================================================

# 检查 celltype 是否只有那两类
print(table(seu_obj_combined$celltype))

# 检查对象结构
print(seu_obj_combined)

# 如果你需要做后续分析（如重新标准化、降维），通常需要重新跑标准流程：
# seu_obj_combined <- NormalizeData(seu_obj_combined)
# seu_obj_combined <- FindVariableFeatures(seu_obj_combined)
# seu_obj_combined <- ScaleData(seu_obj_combined)
# seu_obj_combined <- RunPCA(seu_obj_combined)
# ...

# 保存结果（可选）
saveRDS(seu_obj_combined, "/home/lin/c_group/Treg_CD8_combined.rds")
