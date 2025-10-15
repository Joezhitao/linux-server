rm(list = ls())
# ---------------------------
# 1. 加载包（自动安装缺失包）
# ---------------------------
packages <- c("Seurat", "data.table", "here")
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ---------------------------
# 2. 路径管理（更改为你的数据路径）
# ---------------------------
data_dir <- "/home/lin/CRC/GSE/GSE200997"
count_file <- file.path(data_dir, "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz")
anno_file  <- file.path(data_dir, "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz")

# ---------------------------
# 3. 读取表达矩阵（兼容多种格式）
# ---------------------------
# 支持csv和tsv，自动识别分隔符
read_count_matrix <- function(file) {
  # 判断分隔符
  first_line <- readLines(file, n = 1)
  sep <- ifelse(grepl(",", first_line), ",", "\t")
  dt <- data.table::fread(file, data.table = FALSE, sep = sep)
  rownames(dt) <- dt[, 1]
  dt <- dt[, -1]
  as.matrix(dt)
}

expr_mat <- read_count_matrix(count_file)

# ---------------------------
# 4. 创建Seurat对象
# ---------------------------
seurat_obj <- Seurat::CreateSeuratObject(counts = expr_mat)

# ---------------------------
# 5. 读取注释并自动匹配元数据
# ---------------------------
read_metadata <- function(file, seurat_obj) {
  # 判断分隔符
  first_line <- readLines(file, n = 1)
  sep <- ifelse(grepl(",", first_line), ",", "\t")
  meta <- data.table::fread(file, data.table = FALSE, sep = sep, header = TRUE)
  # 设细胞barcode为行名（假设第一列为barcode）
  rownames(meta) <- meta[, 1]
  meta <- meta[, -1, drop = FALSE]
  # 只保留与Seurat对象匹配的细胞
  meta <- meta[colnames(seurat_obj), , drop = FALSE]
  return(meta)
}

meta <- read_metadata(anno_file, seurat_obj)

# ---------------------------
# 6. 添加元数据到Seurat对象
# ---------------------------
seurat_obj <- Seurat::AddMetaData(seurat_obj, meta)

# ---------------------------
# 7. 检查结果
# ---------------------------
print(seurat_obj)
head(seurat_obj@meta.data)

levels(as.factor(seurat_obj@meta.data$orig.ident))
