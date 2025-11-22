# 加载单细胞RNA-seq分析所需的R包
suppressPackageStartupMessages({
  library(Seurat)       # 单细胞RNA-seq分析的核心包
  library(harmony)      # 用于批次效应校正
  library(tidyverse)    # 数据处理和可视化工具集
  library(dplyr)        # 数据操作和转换
  library(patchwork)    # 组合多个ggplot图形
  library(ggplot2)      # 数据可视化
  library(cowplot)      # 排列ggplot图形
  library(DoubletFinder) # 检测单细胞数据中的双胞体
  library(DropletUtils) # 用于液滴测序数据的处理，包括空液滴检测
  library(scater)       # 单细胞数据质控和可视化
  library(scran)        # 单细胞RNA数据归一化和分析
})

# 清理R环境中的所有变量，确保分析从干净的环境开始
rm(list = ls())

###############################################################################
## 1、批量读取数据并进行空液滴检测
###############################################################################

# 获取目录列表 - 查找所有样本文件夹
dir_name <- list.files('/home/lin/c_group/C_scRNA/')
print("发现的样本目录:")
print(dir_name)

# 初始化两个重要数据结构：
# 1. scRNAlist - 用于存储所有样本的Seurat对象
# 2. loading_stats - 用于记录每个样本的处理统计信息
scRNAlist <- list()
loading_stats <- data.frame(
  Sample = character(),        # 样本名称
  Raw_Droplets = numeric(),    # 原始液滴数量
  Called_Cells = numeric(),    # 鉴定为细胞的液滴数量
  Cell_Rate = numeric(),       # 细胞率(%)
  Method_Used = character(),   # 使用的细胞鉴定方法
  stringsAsFactors = FALSE
)

print("开始批量读取数据并进行空液滴检测...")

# 循环处理每个样本目录
for(i in 1:length(dir_name)){
  cat(paste0("处理样本 ", i, ": ", dir_name[i], "\n"))
  
  # 读取10X Genomics格式的原始计数矩阵
  counts <- Read10X(data.dir = file.path("/home/lin/c_group/C_scRNA/", dir_name[i]))
  raw_droplets <- ncol(counts)  # 记录原始液滴数量
  method_used <- "Direct"       # 默认方法标记
  
  cat(paste0("  原始液滴数: ", raw_droplets, "\n"))
  
  # 只对大样本(>1000个液滴)尝试空液滴检测
  if(raw_droplets > 1000) {
    cat("  尝试空液滴检测...\n")
    
    tryCatch({
      # 计算每个液滴的总UMI计数
      total_counts <- Matrix::colSums(counts)
      
      # 评估数据是否适合空液滴检测：需要足够的低计数和高计数液滴
      low_count_droplets <- sum(total_counts <= 100 & total_counts > 0)
      high_count_droplets <- sum(total_counts > 1000)
      
      cat(paste0("    低计数液滴 (≤100): ", low_count_droplets, "\n"))
      cat(paste0("    高计数液滴 (>1000): ", high_count_droplets, "\n"))
      
      # 如果同时有足够的低计数和高计数液滴，使用emptyDrops算法
      if(low_count_droplets >= 100 && high_count_droplets >= 100) {
        # 运行emptyDrops算法，基于环境RNA背景分布识别真实细胞
        e.out <- emptyDrops(counts, 
                            lower = 100,      # 低于此UMI计数的液滴被视为环境RNA
                            niters = 3000,    # 统计检验的迭代次数
                            test.ambient = TRUE)  # 测试环境RNA分布
        
        # 选择显著的细胞(FDR <= 1%)
        is_cell <- e.out$FDR <= 0.01 & !is.na(e.out$FDR)
        
        # 如果检测到足够的细胞(>50)，使用emptyDrops结果
        if(sum(is_cell) > 50) {
          counts <- counts[, is_cell]  # 只保留被识别为细胞的列
          method_used <- "EmptyDrops"  # 记录使用的方法
          cat(paste0("    空液滴检测成功，保留 ", sum(is_cell), " 个细胞\n"))
        } else {
          # 空液滴检测结果太少，使用简单UMI阈值法作为备用
          cat("    空液滴检测结果过少，使用UMI阈值方法\n")
          counts <- counts[, total_counts > 1000]
          method_used <- "UMI_Threshold_1000"
        }
      } else {
        # 数据分布不适合空液滴检测，使用UMI阈值法
        cat("    数据不适合空液滴检测，使用UMI阈值方法\n")
        counts <- counts[, total_counts > 1000]
        method_used <- "UMI_Threshold_1000"
      }
    }, error = function(e) {
      # 捕获并处理空液滴检测过程中的任何错误
      cat(paste0("    空液滴检测失败: ", e$message, "\n"))
      cat("    使用UMI阈值方法作为备用\n")
      
      total_counts <- Matrix::colSums(counts)
      
      # 根据数据分布动态选择合适的UMI阈值
      threshold <- if(median(total_counts) > 2000) 1000 
      else if(median(total_counts) > 1000) 500 
      else 200
      
      # <<- 操作符用于在错误处理函数中修改外部变量
      counts <<- counts[, total_counts > threshold]
      method_used <<- paste0("UMI_Threshold_", threshold)
    })
  } else {
    # 对于小样本(<1000个液滴)，直接使用UMI阈值法
    cat("  样本较小，直接使用UMI阈值方法\n")
    total_counts <- Matrix::colSums(counts)
    
    # 为小样本选择较低的阈值
    threshold <- ifelse(median(total_counts) > 500, 500, 200)
    counts <- counts[, total_counts > threshold]
    method_used <- paste0("UMI_Threshold_", threshold, "_Small")
  }
  
  # 记录最终鉴定的细胞数量和比例
  called_cells_final <- ncol(counts)
  cell_rate <- round((called_cells_final / raw_droplets) * 100, 2)
  
  # 最终安全检查：如果细胞数过少，尝试使用更低的阈值重新处理
  if(called_cells_final < 50) {
    cat("  警告：细胞数过少，尝试降低阈值\n")
    # 重新读取数据并使用更低的阈值(100 UMIs)
    counts <- Read10X(data.dir = file.path("/home/lin/c_group/C_scRNA/", dir_name[i]))
    total_counts <- Matrix::colSums(counts)
    counts <- counts[, total_counts > 100]
    called_cells_final <- ncol(counts)
    cell_rate <- round((called_cells_final / raw_droplets) * 100, 2)
    method_used <- paste0(method_used, "_Lowered")
  }
  
  # 创建Seurat对象，设置基本质控标准
  if(called_cells_final > 0) {
    scRNAlist[[i]] <- CreateSeuratObject(
      counts, 
      project = dir_name[i],
      min.cells = 3,     # 只保留在至少3个细胞中表达的基因
      min.features = 200 # 只保留表达至少200个基因的细胞
    )
    
    cat(paste0("  成功创建Seurat对象\n"))
  } else {
    # 如果没有检测到有效细胞，记录警告并设为NULL
    warning(paste0("样本 ", i, " 没有检测到有效细胞"))
    scRNAlist[[i]] <- NULL
  }
  
  # 将本样本的处理结果添加到统计表中
  loading_stats <- rbind(loading_stats, data.frame(
    Sample = dir_name[i],
    Raw_Droplets = raw_droplets,
    Called_Cells = called_cells_final,
    Cell_Rate = cell_rate,
    Method_Used = method_used
  ))
  
  cat(paste0("  最终细胞数: ", called_cells_final, " (", cell_rate, "%) - ", method_used, "\n"))
  cat(paste(rep("-", 40), collapse = ""), "\n")  # 打印分隔线
}

# 从列表中移除空(NULL)样本
scRNAlist <- scRNAlist[!sapply(scRNAlist, is.null)]

# 显示所有样本的数据加载统计
print("数据加载统计:")
print(loading_stats)

# 显示不同细胞鉴定方法的使用频率
print("\n使用方法统计:")
print(table(loading_stats$Method_Used))


###############################################################################
## 2、设置物种特定参数
###############################################################################

# 物种配置 - 设置当前分析的目标物种为小鼠
SPECIES <- "rat"

# 物种基因注释配置 - 为不同物种定义基因识别模式和重要基因列表
species_config <- list(
  human = list(
    mt_pattern = "^MT-",       # 人类线粒体基因前缀模式
    ribo_pattern = "^RP[SL]",  # 人类核糖体基因前缀模式
    hb_genes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"), # 人类血红蛋白基因
    species_name = "Human"     # 物种显示名称
  ),
  mouse = list(
    mt_pattern = "^mt-",       # 小鼠线粒体基因前缀模式
    ribo_pattern = "^Rp[sl]",  # 小鼠核糖体基因前缀模式
    hb_genes = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh1", "Hbb-bh2", "Hbb-y", "Hbe1"), # 小鼠血红蛋白基因
    species_name = "Mouse"     # 物种显示名称
  ),
  rat = list(
    mt_pattern = "^Mt-",        # 大鼠线粒体基因前缀模式
    ribo_pattern = "^Rp[sl]",  # 大鼠核糖体基因前缀模式
    hb_genes = c("Hba1", "Hba2", "Hbb"), # 大鼠血红蛋白基因
    species_name = "Rat"       # 物种显示名称
  ),
  zebrafish = list(
    mt_pattern = "^mt",        # 斑马鱼线粒体基因前缀模式
    ribo_pattern = "^rp[sl]",  # 斑马鱼核糖体基因前缀模式
    hb_genes = c("hbaa1", "hbaa2", "hbae1", "hbae3", "hbba1", "hbba2", "hbbe1"), # 斑马鱼血红蛋白基因
    species_name = "Zebrafish" # 物种显示名称
  )
)

# 验证物种配置 - 确保选择的物种在配置中存在
current_config <- species_config[[SPECIES]]
if(is.null(current_config)) {
  stop(paste0("不支持的物种: '", SPECIES, "'! 支持的物种: ", 
              paste(names(species_config), collapse = ", ")))
}

# 显示当前配置 - 输出选定物种的配置信息
cat("=== 当前物种配置 ===\n")
cat("物种:", current_config$species_name, "\n")
cat("线粒体基因模式:", current_config$mt_pattern, "\n")
cat("核糖体基因模式:", current_config$ribo_pattern, "\n")
cat("血红蛋白基因数量:", length(current_config$hb_genes), "\n")
cat("血红蛋白基因:", paste(current_config$hb_genes[1:min(5, length(current_config$hb_genes))], collapse = ", "))
if(length(current_config$hb_genes) > 5) cat("...") # 如果基因过多，只显示前5个
cat("\n")
cat("========================\n")

###############################################################################
## 3、计算质控指标
###############################################################################

print("开始计算质控指标...")

# 质控统计表 - 初始化用于存储所有样本质控指标的数据框
qc_summary <- data.frame(
  Sample = 1:length(scRNAlist),                                 # 样本编号
  Sample_Name = sapply(scRNAlist, function(x) x@project.name),  # 样本名称
  Cells = sapply(scRNAlist, ncol),                              # 细胞数量
  MT_genes_found = NA,                                          # 发现的线粒体基因数
  HB_genes_found = NA,                                          # 发现的血红蛋白基因数
  Ribo_genes_found = NA,                                        # 发现的核糖体基因数
  Mean_MT_percent = NA,                                         # 平均线粒体基因表达比例
  Mean_HB_percent = NA,                                         # 平均血红蛋白基因表达比例
  Mean_Ribo_percent = NA,                                       # 平均核糖体基因表达比例
  Median_nUMI = NA,                                             # 中位数UMI计数
  Median_nGene = NA,                                            # 中位数基因数
  stringsAsFactors = FALSE
)

# 处理每个样本
for(i in 1:length(scRNAlist)) {
  cat(paste0("处理样本 ", i, "/", length(scRNAlist), " (", current_config$species_name, ")...\n"))
  
  seu <- scRNAlist[[i]]  # 获取当前样本的Seurat对象
  
  # 获取基因名称列表
  all_genes <- rownames(seu)
  
  # 根据物种配置识别特定类型的基因
  mt_genes <- grep(current_config$mt_pattern, all_genes, value = TRUE, ignore.case = TRUE)  # 线粒体基因
  ribo_genes <- grep(current_config$ribo_pattern, all_genes, value = TRUE, ignore.case = TRUE)  # 核糖体基因
  hb_genes_present <- intersect(current_config$hb_genes, all_genes)  # 血红蛋白基因
  
  # 如果直接匹配未找到血红蛋白基因，尝试不区分大小写匹配
  if(length(hb_genes_present) == 0) {
    hb_genes_present <- all_genes[toupper(all_genes) %in% toupper(current_config$hb_genes)]
  }
  
  # 显示找到的基因数量
  cat(paste0("  发现基因: MT(", length(mt_genes), "), HB(", length(hb_genes_present), 
             "), Ribo(", length(ribo_genes), ")\n"))
  
  # 显示找到的基因示例（前几个）
  if(length(mt_genes) > 0) {
    cat(paste0("    MT基因示例: ", paste(head(mt_genes, 3), collapse = ", "), "\n"))
  }
  if(length(hb_genes_present) > 0) {
    cat(paste0("    HB基因示例: ", paste(head(hb_genes_present, 3), collapse = ", "), "\n"))
  }
  
  # 计算基本质控指标：每个细胞的UMI总数和基因数
  seu$nUMI <- Matrix::colSums(GetAssayData(seu, layer = "counts"))
  seu$nGene <- Matrix::colSums(GetAssayData(seu, layer = "counts") > 0)
  
  # 计算线粒体基因表达比例
  if(length(mt_genes) > 0) {
    seu <- PercentageFeatureSet(seu, pattern = current_config$mt_pattern, col.name = "mt_percent")
  } else {
    warning(paste0("样本 ", i, " 中未发现线粒体基因，请检查基因命名格式"))
    seu$mt_percent <- 0
  }
  
  # 计算血红蛋白基因表达比例
  if(length(hb_genes_present) > 0) {
    seu <- PercentageFeatureSet(seu, features = hb_genes_present, col.name = "HB_percent")
  } else {
    cat(paste0("    样本 ", i, " 中未发现血红蛋白基因（正常情况）\n"))
    seu$HB_percent <- 0
  }
  
  # 计算核糖体基因表达比例
  if(length(ribo_genes) > 0) {
    seu <- PercentageFeatureSet(seu, pattern = current_config$ribo_pattern, col.name = "ribo_percent")
  } else {
    warning(paste0("样本 ", i, " 中未发现核糖体基因，请检查基因命名格式"))
    seu$ribo_percent <- 0
  }
  
  # 计算衍生指标：基因复杂度和log10转换的计数值
  seu$gene_complexity <- seu$nGene / seu$nUMI  # 基因复杂度（基因数/UMI数）
  seu$log10_nUMI <- log10(seu$nUMI)            # log10转换的UMI计数
  seu$log10_nGene <- log10(seu$nGene)          # log10转换的基因数
  
  # 更新Seurat对象到列表
  scRNAlist[[i]] <- seu
  
  # 更新质控统计表
  qc_summary[i, c("MT_genes_found", "HB_genes_found", "Ribo_genes_found",
                  "Mean_MT_percent", "Mean_HB_percent", "Mean_Ribo_percent",
                  "Median_nUMI", "Median_nGene")] <- 
    c(length(mt_genes), length(hb_genes_present), length(ribo_genes),
      round(mean(seu$mt_percent), 2), round(mean(seu$HB_percent), 2), 
      round(mean(seu$ribo_percent), 2),
      round(median(seu$nUMI), 0), round(median(seu$nGene), 0))
  
  # 输出关键质控指标摘要
  cat(paste0("  质控指标 - MT:", round(mean(seu$mt_percent), 2), "%, ",
             "HB:", round(mean(seu$HB_percent), 2), "%, ",
             "Ribo:", round(mean(seu$ribo_percent), 2), "%\n"))
  cat(paste0("  中位数 - UMI:", median(seu$nUMI), ", Gene:", median(seu$nGene), "\n"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
}

print("质控指标计算完成!")
print("\n=== 质控统计汇总 ===")
print(qc_summary)

# 保存质控统计结果到CSV文件，文件名包含物种和日期信息
write.csv(qc_summary, paste0("/home/lin/c_group/C_scRNA/QC_summary_", SPECIES, "_", Sys.Date(), ".csv"), row.names = FALSE)
cat(paste0("\n质控统计已保存至: /home/lin/c_group/C_scRNA/QC_summary_", SPECIES, "_", Sys.Date(), ".csv\n"))

###############################################################################
## 4、质控前可视化
###############################################################################

print("绘制质控前小提琴图...")

# 初始化列表存储每个样本的小提琴图
violin_before <- list()

# 为每个样本创建质控小提琴图
for(i in 1:length(scRNAlist)){
  # 检查样本是否存在且有细胞 - 避免处理空样本
  if(!is.null(scRNAlist[[i]]) && ncol(scRNAlist[[i]]) > 0) {
    
    # 获取样本名称用于图表标题
    sample_name <- scRNAlist[[i]]@project.name
    
    # 创建多面板小提琴图，展示5个关键质控指标
    violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                  features = c("nFeature_RNA", "nCount_RNA", 
                                               "mt_percent", "HB_percent", "ribo_percent"), 
                                  pt.size = 0.01,  # 设置点的大小较小以减少视觉干扰
                                  ncol = 5) +      # 5个指标排成一行
      plot_annotation(title = paste0("Sample ", i, " (", sample_name, ") - 质控前"),
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
    
    cat(paste0("完成样本 ", i, " 的质控图\n"))
  } else {
    # 跳过无效样本
    cat(paste0("跳过样本 ", i, " (无有效数据)\n"))
    violin_before[[i]] <- NULL
  }
}

# 移除列表中的空元素（对应无效样本）
violin_before <- violin_before[!sapply(violin_before, is.null)]

# 显示所有样本的质控图
print("=== 质控前小提琴图 ===")
for(i in 1:length(violin_before)) {
  print(violin_before[[i]])  # 在R Markdown中逐个显示图形
}

# 输出生成的图形数量
print(paste0("共生成 ", length(violin_before), " 个样本的质控图"))

###############################################################################
## 5、质控过滤
###############################################################################

print("使用简单固定阈值过滤...")

# 初始化数据框用于保存过滤前后的细胞数量统计
filter_stats <- data.frame(
  Sample = character(),      # 样本名称
  Before = numeric(),        # 过滤前细胞数
  After = numeric(),         # 过滤后细胞数
  Retention_Rate = numeric(),# 细胞保留率(%)
  stringsAsFactors = FALSE
)

# 对每个样本进行质控过滤
for(i in 1:length(scRNAlist)) {
  x <- scRNAlist[[i]]        # 获取当前样本
  cells_before <- ncol(x)    # 记录过滤前细胞数
  
  # 使用固定的宽松阈值进行细胞过滤
  # 这些阈值基于常见单细胞数据的经验值设置
  x_filtered <- subset(x, 
                       subset = nFeature_RNA > 200 &     # 至少表达200个基因
                         nFeature_RNA < 8000 &     # 基因数上限，避免双胞体
                         nCount_RNA > 500 &        # 最低UMI计数要求
                         nCount_RNA < 50000 &      # UMI计数上限，避免异常细胞
                         mt_percent < 20 &         # 线粒体比例上限，过高表示细胞死亡
                         HB_percent < 10)          # 血红蛋白比例上限，避免红细胞污染
  
  # 计算过滤后的细胞数和保留率
  cells_after <- ncol(x_filtered)
  retention_rate <- round((cells_after / cells_before) * 100, 2)
  
  # 更新过滤统计信息
  filter_stats <- rbind(filter_stats, data.frame(
    Sample = x@project.name,
    Before = cells_before,
    After = cells_after,
    Retention_Rate = retention_rate
  ))
  
  # 输出每个样本的过滤结果
  cat(paste0("样本 ", i, ": ", cells_before, " -> ", cells_after, 
             " (", retention_rate, "%)\n"))
  
  # 用过滤后的对象替换原对象
  scRNAlist[[i]] <- x_filtered
}

# 输出过滤完成信息和统计表
print("质控过滤完成!")
print(filter_stats)

###############################################################################
## 6、Doublet检测
###############################################################################

# 根据细胞数量估计doublet率 - 10X Genomics官方建议的双胞体率估计函数
# 细胞数越多，双胞体率越高，这与微流控芯片上的细胞捕获原理相关
numb <- function(x) {
  if (x <= 2000) return(0.008)       # 0.8%
  else if (x <= 4000) return(0.016)  # 1.6%
  else if (x <= 6000) return(0.024)  # 2.4%
  else if (x <= 8000) return(0.032)  # 3.2%
  else if (x <= 10000) return(0.04)  # 4.0%
  else if (x <= 12000) return(0.048) # 4.8%
  else if (x <= 14000) return(0.056) # 5.6%
  else if (x <= 16000) return(0.064) # 6.4%
  else if (x <= 18000) return(0.072) # 7.2%
  else return(0.08)                  # 8.0%
}

# 初始化Doublet检测统计表
doublet_summary <- data.frame(
  Sample = character(),            # 样本名称
  Cells_Before = numeric(),        # 检测前细胞数
  Cells_After = numeric(),         # 检测后细胞数
  Doublets_Removed = numeric(),    # 移除的双胞体数
  Doublet_Rate = numeric(),        # 实际双胞体率(%)
  Expected_Rate = numeric(),       # 预期双胞体率(%)
  Optimal_pK = numeric(),          # 最优pK参数值
  Homotypic_Prop = numeric(),      # 同源双胞体比例(%)
  Processing_Time = numeric(),     # 处理时间(分钟)
  stringsAsFactors = FALSE
)

print("开始doublet检测...")

# 对每个样本进行双胞体检测
for(i in 1:length(scRNAlist)){
  start_time <- Sys.time()  # 记录开始时间
  
  cat(paste0("样本 ", i, "/", length(scRNAlist), " doublet检测中...\n"))
  
  seu <- scRNAlist[[i]]
  cells_before <- ncol(seu)
  
  # 跳过细胞数过少的样本，这些样本双胞体率通常较低，且检测可能不准确
  if(cells_before < 100) {
    warning(paste0("样本 ", i, " 细胞数过少，跳过doublet检测"))
    next
  }
  
  # 数据预处理：标准化、特征选择、降维
  suppressWarnings({
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", 
                         scale.factor = 10000, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", 
                                nfeatures = 2000, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  })
  
  # 根据细胞数量计算预期双胞体率和数量
  doublet_rate <- numb(cells_before)
  nExp <- round(doublet_rate * cells_before)
  
  # 进行聚类，为同源双胞体估计做准备
  suppressWarnings({
    seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
    # 根据细胞数量动态调整聚类分辨率
    resolution <- ifelse(cells_before < 1000, 0.3, 
                         ifelse(cells_before < 3000, 0.5, 0.8))
    seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  })
  
  # 计算同源双胞体比例并调整预期双胞体数量
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)  # 估计同源双胞体比例
  nExp_adj <- round(nExp * (1 - homotypic.prop))  # 调整预期双胞体数量
  # 确保调整后的双胞体数量在合理范围内
  nExp_adj <- max(1, min(nExp_adj, round(cells_before * 0.12)))
  
  # 寻找最优pK值（影响检测灵敏度的关键参数）
  optimal_pk <- 0.09  # 默认值
  tryCatch({
    sweep.res.list <- paramSweep(seu, PCs = 1:30, sct = FALSE)  # 参数扫描
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)   # 汇总结果
    bcmvn <- find.pK(sweep.stats)  # 找到最佳pK值
    optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  }, error = function(e) {
    cat("  pK优化失败，使用默认值\n")
  })
  
  # 运行DoubletFinder算法
  seu <- doubletFinder(seu = seu, 
                       PCs = 1:30,     # 使用前30个主成分
                       pN = 0.25,      # 模拟双胞体比例
                       pK = optimal_pk,# 最优pK值
                       nExp = nExp_adj,# 调整后的预期双胞体数量
                       sct = FALSE)    # 不使用SCTransform
  
  # 获取DoubletFinder结果并添加到元数据
  df_results <- colnames(seu@meta.data)[grep("DF.classifications", colnames(seu@meta.data))]
  seu@meta.data$doublet_info <- seu@meta.data[[df_results]]
  
  # 过滤出单细胞（移除双胞体）
  seu_filtered <- subset(seu, subset = doublet_info == "Singlet")
  
  # 统计结果
  cells_after <- ncol(seu_filtered)
  doublets_removed <- cells_before - cells_after
  actual_doublet_rate <- round(doublets_removed / cells_before * 100, 2)
  processing_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # 更新Seurat对象列表
  scRNAlist[[i]] <- seu_filtered
  
  # 更新统计表
  doublet_summary <- rbind(doublet_summary, data.frame(
    Sample = seu@project.name,
    Cells_Before = cells_before,
    Cells_After = cells_after,
    Doublets_Removed = doublets_removed,
    Doublet_Rate = actual_doublet_rate,
    Expected_Rate = round(doublet_rate * 100, 1),
    Optimal_pK = round(optimal_pk, 3),
    Homotypic_Prop = round(homotypic.prop * 100, 1),
    Processing_Time = round(processing_time, 2)
  ))
  
  # 输出处理结果
  cat(paste0("  完成: ", cells_before, " -> ", cells_after, 
             " (", actual_doublet_rate, "%, ", 
             round(processing_time, 1), "min)\n"))
}

print("Doublet检测完成!")
print(doublet_summary)

###############################################################################
## 7、合并所有样本
###############################################################################

# 合并所有样本 - 将多个样本的数据整合到一个Seurat对象中
print("合并所有样本...")

# 判断是否有多个样本需要合并
if(length(scRNAlist) > 1) {
  # 使用Seurat的merge函数合并样本
  # x参数为第一个样本，y参数为剩余所有样本
  # add.cell.ids参数为每个样本添加唯一前缀，防止细胞ID冲突
  scRNAlist_merge <- merge(x = scRNAlist[[1]], 
                           y = scRNAlist[-1], 
                           add.cell.ids = sapply(scRNAlist, function(x) x@project.name))
  
  # 输出合并结果信息
  print(paste0("成功合并 ", length(scRNAlist), " 个样本，共 ", ncol(scRNAlist_merge), " 个细胞"))
} else {
  # 如果只有一个样本，直接使用该样本
  scRNAlist_merge <- scRNAlist[[1]]
  print("只有一个样本，无需合并")
}

# 统计合并后各样本的细胞数量分布
sample_counts <- table(scRNAlist_merge$orig.ident)
print("合并后各样本细胞数:")
print(sample_counts)

# 将样本标识转换为因子类型，便于后续分析和可视化
# 这确保了样本在图表和分析中的顺序是可控的
scRNAlist_merge@meta.data$orig.ident <- as.factor(scRNAlist_merge@meta.data$orig.ident)
print("样本组别:")
print(levels(scRNAlist_merge@meta.data$orig.ident))

# 保存RDS文件
saveRDS(scRNAlist_merge, file = "/home/lin/c_group/C_group_11m14d.rds")
