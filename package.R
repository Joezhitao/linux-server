# 安装所需包的函数
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(paste("安装包:", package))
    
    # 对于特殊的Bioconductor包使用BiocManager
    if (package %in% c("celldex", "SingleR")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package, update = FALSE)
    } 
    # 对于GitHub上的包使用devtools
    else if (package == "scrapper") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_github("YosefLab/scrapper")
    }
    # 对于tidydr包，可能是CRAN或GitHub上的
    else if (package == "tidydr") {
      tryCatch({
        install.packages("tidydr")
      }, error = function(e) {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          install.packages("devtools")
        }
        devtools::install_github("cran/tidydr")
      })
    }
    # 其他标准CRAN包
    else {
      install.packages(package)
    }
  } else {
    message(paste("包已安装:", package))
  }
}

# 需要安装的包列表
required_packages <- c(
  # 基础分析包
  "Seurat",
  
  # 细胞类型注释相关包
  "celldex",
  "SingleR",
  "scrapper",
  
  # 数据处理和可视化包
  "dplyr",
  "ggsci",
  "tidydr",
  
  # 其他依赖包
  "ggplot2",
  "cowplot",
  "harmony",
  "clustree",
  "kableExtra",
  "readr",
  "here"
)

# 安装BiocManager（如果需要）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 安装每个包
message("开始安装所需的R包...")
for (pkg in required_packages) {
  install_if_missing(pkg)
}

# 验证安装
missing_packages <- character()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  message("所有包安装成功！")
} else {
  warning("以下包安装失败，请手动安装: ", 
          paste(missing_packages, collapse = ", "))
}
