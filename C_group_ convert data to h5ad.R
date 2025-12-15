library(reticulate)
# 1. 设置超时时间为 300秒 (5分钟)，解决网络慢的问题
Sys.setenv(UV_HTTP_TIMEOUT = "300")
# 1. 检查当前 R 使用的是哪个 Python
py_config() 
# (注意看输出的路径，比如 /usr/bin/python3 或者 ~/.virtualenvs/r-reticulate/bin/python)

# 2. 尝试直接通过 R 安装 anndata 到该环境中
py_install("anndata", pip = TRUE)

library(sceasy)
library(Seurat)
pbmc <- readRDS("/home/lin/c_group/hep.rds")
# 合并RNA层数据
# 将需要保留在 H5ad 文件中的主要数据层（通常是原始 counts 矩阵）整合，以确保数据被正确导出
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])

# 降级Seurat v5格式以确保兼容性，降级是Seurat转出h5ad必须的
pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")

# 转换为h5ad格式
pbmc = system.time({
  sceasy::convertFormat(pbmc, 
                        from = "seurat", 
                        to = "anndata",
                        outFile = '/home/lin/c_group/New Folder/scanpy/hep.h5ad')
})
print(pbmc)