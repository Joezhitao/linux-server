# ======================================================
# R语言: Seurat 转 h5ad (修正版 - 强制导出 Raw Counts)
# ======================================================

rm(list = ls())
gc()

library(sceasy)
library(reticulate)
library(Seurat)

# 1. 加载数据
setwd("/home/lin/c_group/")
pbmc <- readRDS("/home/lin/c_group/hep.rds")

# 2. 确保默认 Assay 是 RNA (不要用 Integrated 或 SCT)
DefaultAssay(pbmc) <- "RNA"

# 3. 降级 Seurat v5 格式 (关键步骤，保持你之前的操作)
# 这一步是为了让 reticulate 能读懂 Seurat 对象
if (class(pbmc[["RNA"]])[1] == "Assay5") {
  pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")
}

# 4. 转换格式 (核心修改在这里！！！)
# main_layer = "counts" : 告诉它把 Seurat 的 counts 槽放到 Python 的 .X 里
# transfer_layers : 确保 data (lognorm) 和 counts 都被转过去
cgroup = system.time({
  sceasy::convertFormat(
    pbmc, 
    from = "seurat", 
    to = "anndata",
    assay = "RNA",               # 指定 RNA assay
    main_layer = "counts",       # 【重要】将 Raw Counts 放入 .X
    transfer_layers = c("counts", "data"), # 同时保留 counts 和 data
    outFile = '/home/lin/c_group/New Folder/scanpy/hep.h5ad' # 覆盖原文件
  )
})

print("转换完成！耗时：")
print(cgroup)
