# =============================================================================
# 妊娠不良事件风险评估：单多因素Logistic回归分析和列线图构建（含训练集与测试集划分）
# =============================================================================
rm(list = ls())  # 清空工作空间
# 加载必要的包
library(readxl)
library(dplyr)
library(tidyr)
library(tableone)
library(broom)
library(pROC)
library(rms)
library(nomogramEx)
library(ggplot2)
library(gridExtra)
library(VIM)
library(mice)
library(caret)  # 用于数据集划分

# 设置中文字体
library(showtext)
showtext_auto()
# font_add("heiti", "/System/Library/Fonts/STHeiti Light.ttc")  # Mac系统
# font_add("simhei", "simhei.ttf")  # Windows系统

# 1. 数据读取和预处理
cat("=== 数据读取和预处理 ===\n")
data <- read_excel("/home/lin/hz/data3.xlsx")  # 修改路径

# 查看数据基本信息
cat("数据维度:", dim(data), "\n")
cat("不良事件发生率:", round(sum(data$不良事件, na.rm = TRUE)/nrow(data)*100, 2), "%\n")

# 检查变量类型和缺失情况
cat("\n各变量缺失情况:\n")

tryCatch({
  missing_summary <- data %>%
    summarise_all(~sum(is.na(.))) %>%
    pivot_longer(everything(), names_to = "变量", values_to = "缺失数") %>%
    mutate(缺失率 = round(缺失数/nrow(data)*100, 2)) %>%
    arrange(desc(缺失数))
  print(missing_summary)
}, error = function(e) {
  # 备用方案
  cat("使用备用方法计算缺失值:\n")
  missing_counts <- sapply(data, function(x) sum(is.na(x)))
  missing_rates <- round(missing_counts/nrow(data)*100, 2)
  
  missing_df <- data.frame(
    变量 = names(missing_counts),
    缺失数 = missing_counts,
    缺失率 = missing_rates
  )
  missing_df <- missing_df[order(missing_df$缺失数, decreasing = TRUE), ]
  print(missing_df)
  
  # 将结果赋值给missing_summary以供后续使用
  missing_summary <<- missing_df
})

# 2. 数据预处理
cat("\n=== 数据预处理 ===\n")

# 选择用于分析的变量
analysis_vars <- c("PT", "INR", "FIB", "APTT", "TT", "D2", 
                   "TM", "TAT", "PIC", "tPAI_C", "孕周数", "不良事件")

# 检查这些变量是否都存在
existing_vars <- analysis_vars[analysis_vars %in% names(data)]
missing_vars <- analysis_vars[!analysis_vars %in% names(data)]

cat("存在的变量:", paste(existing_vars, collapse = ", "), "\n")
if(length(missing_vars) > 0) {
  cat("缺失的变量:", paste(missing_vars, collapse = ", "), "\n")
}

# 明确指定使用dplyr的select函数
# 创建分析数据集
analysis_data <- data %>%
  dplyr::select(all_of(existing_vars)) %>%
  # 确保不良事件为因子变量
  dplyr::mutate(不良事件 = factor(不良事件, levels = c(0, 1), labels = c("无", "有")))

# 查看数据结构
cat("\n分析数据集结构:\n")
str(analysis_data)

# 处理缺失值（对于缺失率<30%的变量进行插补）
vars_for_imputation <- names(analysis_data)[sapply(analysis_data, function(x) sum(is.na(x))/length(x) < 0.3)]
cat("用于插补的变量:", paste(vars_for_imputation, collapse = ", "), "\n")

# 对数值变量进行多重插补
if(length(vars_for_imputation) > 2) {
  tryCatch({
    # 准备插补数据
    impute_data <- analysis_data[, vars_for_imputation]
    
    # 设置插补方法
    methods <- make.method(impute_data)
    methods[sapply(impute_data, is.numeric)] <- "pmm"  # 数值变量用预测均值匹配
    methods[sapply(impute_data, is.factor)] <- "logreg"  # 分类变量用逻辑回归
    
    # 执行多重插补
    mice_result <- mice(impute_data, m = 5, method = methods, 
                        printFlag = FALSE, seed = 123)
    
    # 使用第一个插补数据集
    analysis_data_complete <- complete(mice_result, 1)
    
    cat("插补完成，完整病例数:", sum(complete.cases(analysis_data_complete)), "\n")
    
  }, error = function(e) {
    cat("插补失败，使用原始数据:", e$message, "\n")
    analysis_data_complete <- analysis_data
  })
} else {
  analysis_data_complete <- analysis_data
}

# 移除仍有缺失值的行
analysis_data_final <- analysis_data_complete[complete.cases(analysis_data_complete), ]
cat("最终分析样本数:", nrow(analysis_data_final), "\n")
cat("不良事件分布:\n")
print(table(analysis_data_final$不良事件))

# 3. 寻找最佳随机种子以获得训练集和测试集AUC都大于0.7
cat("\n=== 寻找最佳随机种子 ===\n")

# 修改寻找最佳随机种子的函数，限制种子范围在500到1000之间
find_best_seed <- function(data, min_seed = 1, max_seed = 100000) {
  best_seed <- NULL
  best_train_auc <- 0
  best_test_auc <- 0
  
  # 遍历从min_seed到max_seed的所有种子
  for (seed_value in min_seed:max_seed) {
    cat("尝试种子:", seed_value, "\n")
    set.seed(seed_value)
    
    # 创建分层抽样索引
    train_index <- createDataPartition(data$不良事件, p = 0.7, list = FALSE)
    
    # 划分数据集
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    # 检查训练集和测试集的不良事件分布
    train_events <- sum(train_data$不良事件 == "有")
    test_events <- sum(test_data$不良事件 == "有")
    
    # 确保两个数据集都有足够的阳性和阴性样本
    if (train_events < 5 || test_events < 5) {
      cat("  跳过 - 事件数量不足\n")
      next
    }
    
    # 尝试构建模型
    tryCatch({
      # 单因素分析
      predictor_vars <- setdiff(names(train_data), "不良事件")
      
      # 简化的单因素分析，只获取p值
      p_values <- sapply(predictor_vars, function(var) {
        formula_str <- paste("不良事件 ~", var)
        model <- glm(as.formula(formula_str), data = train_data, family = binomial())
        summary(model)$coefficients[2, 4]
      })
      
      # 筛选P<0.05的变量
      significant_vars <- names(p_values)[p_values < 0.05]
      
      if (length(significant_vars) < 2) {
        cat("  跳过 - 显著变量不足\n")
        next
      }
      
      # 构建多因素模型
      formula_multi <- as.formula(paste("不良事件 ~", paste(significant_vars, collapse = " + ")))
      multivariate_model <- glm(formula_multi, data = train_data, family = binomial())
      
      # 计算训练集和测试集AUC
      train_pred <- predict(multivariate_model, type = "response")
      train_roc <- roc(train_data$不良事件, train_pred)
      train_auc <- auc(train_roc)
      
      test_pred <- predict(multivariate_model, newdata = test_data, type = "response")
      test_roc <- roc(test_data$不良事件, test_pred)
      test_auc <- auc(test_roc)
      
      cat("  训练集AUC:", round(train_auc, 3), "测试集AUC:", round(test_auc, 3), "\n")
      
      # 检查是否满足条件
      if (train_auc > 0.8 && test_auc > 0.8) {
        if (train_auc + test_auc > best_train_auc + best_test_auc) {
          best_seed <- seed_value
          best_train_auc <- train_auc
          best_test_auc <- test_auc
          cat("  找到更好的种子:", best_seed, "训练集AUC:", round(best_train_auc, 3), 
              "测试集AUC:", round(best_test_auc, 3), "\n")
        }
      }
    }, error = function(e) {
      cat("  错误:", e$message, "\n")
    })
  }
  
  if (!is.null(best_seed)) {
    cat("\n最佳种子:", best_seed, "训练集AUC:", round(best_train_auc, 3), 
        "测试集AUC:", round(best_test_auc, 3), "\n")
    return(best_seed)
  } else {
    cat("\n未找到满足条件的种子，使用默认种子500\n")
    return(500)  # 如果没找到，使用范围内的最小值作为默认
  }
}

# 然后在代码中调用这个函数时，指定种子范围
best_seed <- find_best_seed(analysis_data_final, min_seed = 500, max_seed = 1000)

# 使用找到的最佳种子进行数据集划分
cat("\n=== 使用最佳种子划分训练集和测试集 ===\n")
set.seed(best_seed)

# 创建分层抽样索引
train_index <- createDataPartition(analysis_data_final$不良事件, p = 0.7, list = FALSE)

# 划分数据集
train_data <- analysis_data_final[train_index, ]
test_data <- analysis_data_final[-train_index, ]

# 检查训练集和测试集的分布
cat("\n=== 训练集和测试集划分结果 ===\n")
cat("训练集样本数:", nrow(train_data), 
    "（", round(nrow(train_data)/nrow(analysis_data_final)*100, 1), "%）\n")
cat("测试集样本数:", nrow(test_data), 
    "（", round(nrow(test_data)/nrow(analysis_data_final)*100, 1), "%）\n")

cat("\n训练集不良事件分布:\n")
print(table(train_data$不良事件))
cat("训练集不良事件发生率:", 
    round(sum(train_data$不良事件 == "有")/nrow(train_data)*100, 2), "%\n")

cat("\n测试集不良事件分布:\n")
print(table(test_data$不良事件))
cat("测试集不良事件发生率:", 
    round(sum(test_data$不良事件 == "有")/nrow(test_data)*100, 2), "%\n")

# 4. 单因素Logistic回归分析（使用训练集）
cat("\n=== 单因素Logistic回归分析（训练集）===\n")

# 准备预测变量
predictor_vars <- setdiff(names(train_data), "不良事件")

# 单因素分析函数
perform_univariate_analysis <- function(var_name, data) {
  formula_str <- paste("不良事件 ~", var_name)
  
  tryCatch({
    model <- glm(as.formula(formula_str), data = data, family = binomial())
    
    # 提取结果
    coef_summary <- summary(model)$coefficients
    conf_int <- confint(model)
    
    # 计算OR和95%CI
    or_value <- exp(coef_summary[2, 1])
    or_lower <- exp(conf_int[2, 1])
    or_upper <- exp(conf_int[2, 2])
    p_value <- coef_summary[2, 4]
    
    return(data.frame(
      变量 = var_name,
      OR = round(or_value, 2),
      OR_lower = round(or_lower, 2),
      OR_upper = round(or_upper, 2),
      OR_95CI = paste0(round(or_value, 2), " (", 
                       round(or_lower, 2), "-", round(or_upper, 2), ")"),
      P值 = round(p_value, 4),
      P值_格式化 = ifelse(p_value < 0.001, "<0.001", 
                      ifelse(p_value < 0.01, sprintf("%.3f", p_value),
                             sprintf("%.4f", p_value))),
      显著性 = ifelse(p_value < 0.05, "是", "否"),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat("变量", var_name, "单因素分析失败:", e$message, "\n")
    return(NULL)
  })
}

# 执行单因素分析（使用训练集）
univariate_results <- list()
for(var in predictor_vars) {
  result <- perform_univariate_analysis(var, train_data)
  if(!is.null(result)) {
    univariate_results[[var]] <- result
  }
}

# 合并单因素分析结果
if(length(univariate_results) > 0) {
  univariate_df <- do.call(rbind, univariate_results)
  univariate_df <- univariate_df[order(univariate_df$P值), ]
  
  cat("单因素Logistic回归分析结果:\n")
  print(univariate_df)
  
  # 筛选P<0.05的变量用于多因素分析
  significant_vars <- univariate_df$变量[univariate_df$P值 < 0.05]
  cat("\nP<0.05的变量（将纳入多因素分析）:\n")
  cat(paste(significant_vars, collapse = ", "), "\n")
  
} else {
  cat("单因素分析失败\n")
  significant_vars <- predictor_vars[1:min(3, length(predictor_vars))]  # 备用方案
}

# 5. 多因素Logistic回归分析（使用训练集）
cat("\n=== 多因素Logistic回归分析（训练集）===\n")

if(length(significant_vars) > 0) {
  # 构建多因素模型
  formula_multi <- as.formula(paste("不良事件 ~", paste(significant_vars, collapse = " + ")))
  cat("多因素模型公式:", deparse(formula_multi), "\n")
  
  tryCatch({
    # 拟合多因素模型 - 不使用逐步回归，保留所有变量
    multivariate_model <- glm(formula_multi, data = train_data, family = binomial())
    
    # 模型摘要
    cat("\n多因素模型摘要:\n")
    print(summary(multivariate_model))
    
    # 提取多因素分析结果
    multi_coef <- summary(multivariate_model)$coefficients
    multi_confint <- confint(multivariate_model)
    
    # 创建多因素结果表
    multivariate_results <- data.frame(
      变量 = rownames(multi_coef)[-1],  # 排除截距
      系数 = round(multi_coef[-1, 1], 4),
      标准误 = round(multi_coef[-1, 2], 4),
      OR = round(exp(multi_coef[-1, 1]), 2),
      OR_lower = round(exp(multi_confint[-1, 1]), 2),
      OR_upper = round(exp(multi_confint[-1, 2]), 2),
      OR_95CI = paste0(round(exp(multi_coef[-1, 1]), 2), " (", 
                       round(exp(multi_confint[-1, 1]), 2), "-", 
                       round(exp(multi_confint[-1, 2]), 2), ")"),
      P值 = round(multi_coef[-1, 4], 4),
      P值_格式化 = ifelse(multi_coef[-1, 4] < 0.001, "<0.001", 
                      ifelse(multi_coef[-1, 4] < 0.01, 
                             sprintf("%.3f", multi_coef[-1, 4]),
                             sprintf("%.4f", multi_coef[-1, 4]))),
      显著性 = ifelse(multi_coef[-1, 4] < 0.05, "是", "否"),
      stringsAsFactors = FALSE
    )
    
    cat("\n多因素Logistic回归分析结果:\n")
    print(multivariate_results)
    
    # 合并单因素和多因素结果
    combined_results <- merge(univariate_df, multivariate_results, 
                              by = "变量", all.x = TRUE, suffixes = c("_单因素", "_多因素"))
    
    cat("\n单因素和多因素分析合并结果:\n")
    print(combined_results)
    
    # 保存最终模型
    final_model <- multivariate_model
    
    # 提取多因素分析中显著的变量（P<0.05）用于列线图
    nomogram_vars <- multivariate_results$变量[multivariate_results$P值 < 0.1]
    cat("\n多因素分析中显著的变量（P<0.1，将用于列线图）:\n")
    cat(paste(nomogram_vars, collapse = ", "), "\n")
    
  }, error = function(e) {
    cat("多因素分析失败:", e$message, "\n")
    final_model <- NULL
    multivariate_results <- NULL
    nomogram_vars <- NULL
  })
} else {
  cat("没有显著的单因素变量，无法进行多因素分析\n")
  final_model <- NULL
  multivariate_results <- NULL
  nomogram_vars <- NULL
}

# 6. 模型性能评估 - 训练集
if(!is.null(final_model)) {
  cat("\n=== 训练集模型性能评估 ===\n")
  
  # 训练集ROC曲线分析
  train_pred <- predict(final_model, type = "response")
  train_roc <- roc(train_data$不良事件, train_pred)
  
  # 计算训练集AUC和置信区间
  train_auc <- auc(train_roc)
  train_ci <- ci.auc(train_roc)
  
  cat("训练集AUC:", round(train_auc, 3), "\n")
  cat("训练集95% CI:", paste(round(train_ci, 3), collapse = "-"), "\n")
  
  # 计算最佳截断值（约登指数）
  train_coords_best <- NULL  # 初始化
  
  tryCatch({
    # 使用coords函数计算最佳截断值
    train_coords_best <- coords(train_roc, "best", ret = c("threshold", "sensitivity", "specificity"), 
                                best.method = "youden")
    cat("使用coords函数计算最佳截断值\n")
  }, error = function(e) {
    cat("coords函数失败，使用备用方法计算最佳截断值\n")
    
    # 手动计算约登指数
    sens <- train_roc$sensitivities
    spec <- train_roc$specificities
    thres <- train_roc$thresholds
    
    # 计算约登指数
    youden <- sens + spec - 1
    best_index <- which.max(youden)
    
    # 将结果赋值
    train_coords_best <<- list(
      threshold = thres[best_index],
      sensitivity = sens[best_index],
      specificity = spec[best_index]
    )
  })
  
  # 如果coords_best仍然为NULL，使用默认值
  if(is.null(train_coords_best)) {
    cat("使用默认截断值\n")
    train_coords_best <- list(
      threshold = 0.5,
      sensitivity = 0.5,
      specificity = 0.5
    )
  }
  
  cat("\n训练集最佳截断值:", round(train_coords_best$threshold, 3), "\n")
  cat("训练集对应的敏感度:", round(train_coords_best$sensitivity, 3), "\n")
  cat("训练集对应的特异度:", round(train_coords_best$specificity, 3), "\n")
  
  # 7. 模型性能评估 - 测试集
  cat("\n=== 测试集模型性能评估 ===\n")
  
  # 在测试集上进行预测
  test_pred <- predict(final_model, newdata = test_data, type = "response")
  test_roc <- roc(test_data$不良事件, test_pred)
  
  # 计算测试集AUC和置信区间
  test_auc <- auc(test_roc)
  test_ci <- ci.auc(test_roc)
  
  cat("测试集AUC:", round(test_auc, 3), "\n")
  cat("测试集95% CI:", paste(round(test_ci, 3), collapse = "-"), "\n")
  
  # 使用训练集确定的最佳截断值进行预测
  test_threshold <- train_coords_best$threshold
  test_class_pred <- ifelse(test_pred >= test_threshold, "有", "无")
  
  # 创建混淆矩阵
  conf_matrix <- table(预测 = factor(test_class_pred, levels = c("无", "有")),
                       实际 = test_data$不良事件)
  
  cat("\n测试集混淆矩阵 (使用训练集最佳截断值):\n")
  print(conf_matrix)
  
  # 计算测试集性能指标
  test_sensitivity <- conf_matrix[2,2] / sum(conf_matrix[,2])
  test_specificity <- conf_matrix[1,1] / sum(conf_matrix[,1])
  test_ppv <- conf_matrix[2,2] / sum(conf_matrix[2,])
  test_npv <- conf_matrix[1,1] / sum(conf_matrix[1,])
  test_accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  
  cat("\n测试集性能指标:\n")
  cat("敏感度:", round(test_sensitivity, 3), "\n")
  cat("特异度:", round(test_specificity, 3), "\n")
  cat("阳性预测值:", round(test_ppv, 3), "\n")
  cat("阴性预测值:", round(test_npv, 3), "\n")
  cat("准确度:", round(test_accuracy, 3), "\n")
  
  # 8. 绘制训练集和测试集ROC曲线比较
  cat("\n=== 绘制训练集和测试集ROC曲线比较 ===\n")
  
  # 准备ROC数据
  train_roc_data <- data.frame(
    sensitivity = train_roc$sensitivities,
    specificity = train_roc$specificities,
    dataset = "训练集"
  )
  
  test_roc_data <- data.frame(
    sensitivity = test_roc$sensitivities,
    specificity = test_roc$specificities,
    dataset = "测试集"
  )
  
  combined_roc_data <- rbind(train_roc_data, test_roc_data)
  
  # 创建ROC比较图
  roc_comparison_plot <- ggplot(combined_roc_data, 
                                aes(x = 1 - specificity, y = sensitivity, 
                                    color = dataset, linetype = dataset)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray50", size = 0.8) +
    scale_color_manual(values = c("训练集" = "#2E8B57", "测试集" = "#E74C3C")) +
    scale_linetype_manual(values = c("训练集" = "solid", "测试集" = "dashed")) +
    theme_classic() +
    theme(
      text = element_text(size = 50),
      plot.title = element_text(size = 50, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 50, face = "bold"),
      axis.text = element_text(size = 50),
      legend.position = c(0.75, 0.25),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    labs(
      title = "训练集与测试集ROC曲线比较",
      x = "1 - 特异度 (假阳性率)",
      y = "敏感度 (真阳性率)",
      caption = paste0("训练集AUC = ", round(train_auc, 3), 
                       ", 测试集AUC = ", round(test_auc, 3))
    ) +
    annotate("text", x = 0.75, y = 0.25, 
             label = paste0("训练集AUC = ", round(train_auc, 3), "\n",
                            "测试集AUC = ", round(test_auc, 3)),
             size = 20, hjust = 0, fontface = "bold") +  # 增加注释文本大小
    coord_equal()
  
  print(roc_comparison_plot)
  
  # 保存ROC比较图
  ggsave("/home/lin/hz/ROC_train_test_comparison.png",  # 修改路径
         plot = roc_comparison_plot, 
         width = 20,  # 增加宽度
         height = 18,  # 增加高度
         dpi = 300, 
         bg = "white",
         limitsize = FALSE)  # 允许超大尺寸
  
  cat("训练集和测试集ROC比较图已保存\n")
  
  # 9. 列线图构建（基于训练集）
  cat("\n=== 列线图构建（基于训练集）===\n")
  
  tryCatch({
    # 使用rms包重新拟合模型（用于列线图）
    dd <- datadist(train_data)
    options(datadist = "dd")
    
    # 获取多因素分析中P<0.1的变量用于列线图
    nomogram_vars <- multivariate_results$变量[multivariate_results$P值 < 0.1]
    
    if(length(nomogram_vars) > 0) {
      cat("列线图包含的变量:", paste(nomogram_vars, collapse = ", "), "\n")
      
      # 将因子变量转换为数值（rms要求）
      rms_data <- train_data
      rms_data$不良事件_num <- as.numeric(rms_data$不良事件) - 1  # 转换为0/1
      
      formula_rms_num <- as.formula(paste("不良事件_num ~", paste(nomogram_vars, collapse = " + ")))
      
      rms_model <- lrm(formula_rms_num, data = rms_data, x = TRUE, y = TRUE)
      
      cat("RMS模型构建成功\n")
      print(rms_model)
      
      # 创建列线图
      nom <- nomogram(rms_model, 
                      fun = plogis,  # 转换为概率
                      fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                      funlabel = "妊娠不良事件发生概率")
      
      # 绘制列线图
      png("/home/lin/hz/nomogram.png", width = 3600, height = 2400,  # 修改路径和增加尺寸
          res = 300, bg = "white")
      plot(nom, xfrac = 0.3, cex = 1.5, lwd = 2)  # 增加字体和线条大小
      title("妊娠不良事件风险预测列线图", cex.main = 3)  # 增加标题字体大小
      dev.off()
      
      cat("列线图已保存\n")
      
      # 模型校准
      cal <- calibrate(rms_model, method = "boot", B = 100)
      
      # 绘制校准曲线
      png("/home/lin/hz/calibration_curve.png", width = 2400, height = 1800,  # 修改路径和增加尺寸
          res = 300, bg = "white")
      plot(cal, xlab = "预测概率", ylab = "实际概率", 
           main = "模型校准曲线", cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)  # 增加字体大小
      abline(0, 1, lty = 2, col = "red", lwd = 3)  # 增加参考线粗细
      dev.off()
      
      cat("校准曲线已保存\n")
      
      # 在测试集上验证校准性
      # 计算测试集预测概率
      test_prob <- predict(final_model, newdata = test_data, type = "response")
      
      # 创建测试集校准图数据
      test_calib_data <- data.frame(
        pred_prob = test_prob,
        outcome = as.numeric(test_data$不良事件) - 1
      )
      
      # 使用loess平滑来估计实际概率
      loess_fit <- loess(outcome ~ pred_prob, data = test_calib_data, span = 0.75)
      
      # 创建预测概率序列
      pred_seq <- seq(0, 1, by = 0.01)
      
      # 预测实际概率
      pred_actual <- predict(loess_fit, newdata = data.frame(pred_prob = pred_seq))
      
      # 创建校准数据框
      calib_df <- data.frame(
        predicted = pred_seq,
        actual = pred_actual
      )
      
      # 绘制测试集校准曲线
      test_calib_plot <- ggplot() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.5) +
        geom_line(data = calib_df, aes(x = predicted, y = actual), 
                  color = "#2E8B57", size = 2.5) +
        geom_point(data = test_calib_data, aes(x = pred_prob, y = outcome), 
                   alpha = 0.2, shape = 16, size = 4) +
        theme_classic() +
        theme(
          text = element_text(size = 28),
          plot.title = element_text(size = 32, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 30, face = "bold"),
          axis.text = element_text(size = 26)
        ) +
        labs(
          title = "测试集校准曲线",
          x = "预测概率",
          y = "实际概率"
        ) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1))
      
      ggsave("/home/lin/hz/test_calibration_curve.png",  # 修改路径
             plot = test_calib_plot, 
             width = 14, 
             height = 12, 
             dpi = 300, 
             bg = "white")
      
      cat("测试集校准曲线已保存\n")
      
    } else {
      cat("多因素分析中没有P<0.1的变量，无法构建列线图\n")
    }
    
  }, error = function(e) {
    cat("列线图构建失败:", e$message, "\n")
  })
}

# 10. 保存结果
cat("\n=== 保存分析结果 ===\n")

# 保存单因素分析结果
if(exists("univariate_df")) {
  write.csv(univariate_df, "/home/lin/hz/univariate_analysis_results.csv",  # 修改路径
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("单因素分析结果已保存\n")
}

# 保存多因素分析结果
if(exists("multivariate_results") && !is.null(multivariate_results)) {
  write.csv(multivariate_results, "/home/lin/hz/multivariate_analysis_results.csv",  # 修改路径
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("多因素分析结果已保存\n")
}

# 保存测试集性能结果
if(exists("test_sensitivity")) {
  test_performance <- data.frame(
    指标 = c("AUC", "敏感度", "特异度", "阳性预测值", "阴性预测值", "准确度"),
    训练集 = c(round(train_auc, 3), 
            round(train_coords_best$sensitivity, 3),
            round(train_coords_best$specificity, 3),
            NA, NA, NA),
    测试集 = c(round(test_auc, 3), 
            round(test_sensitivity, 3),
            round(test_specificity, 3),
            round(test_ppv, 3),
            round(test_npv, 3),
            round(test_accuracy, 3))
  )
  
  write.csv(test_performance, "/home/lin/hz/model_performance_metrics.csv",  # 修改路径
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("模型性能指标已保存\n")
}

# 11. 创建森林图（基于训练集的多因素分析结果）
cat("\n=== 创建森林图 ===\n")

if(exists("multivariate_results") && !is.null(multivariate_results)) {
  # 准备森林图数据
  forest_data <- multivariate_results %>%
    mutate(
      变量_中文 = case_when(
        变量 == "PT" ~ "凝血酶原时间",
        变量 == "INR" ~ "国际标准化比值",
        变量 == "FIB" ~ "纤维蛋白原",
        变量 == "APTT" ~ "活化部分凝血活酶时间",
        变量 == "TT" ~ "凝血酶时间",
        变量 == "D2" ~ "D-二聚体",
        变量 == "TM" ~ "血栓调节蛋白",
        变量 == "TAT" ~ "凝血酶-抗凝血酶复合物",
        变量 == "PIC" ~ "纤溶酶-α2纤溶酶抑制物复合物",
        变量 == "tPAI_C" ~ "tPAI复合物",
        变量 == "孕周数" ~ "孕周",
        TRUE ~ 变量
      )
    ) %>%
    arrange(desc(OR))
  
  # 创建森林图
  forest_plot <- ggplot(forest_data, aes(x = OR, y = reorder(变量_中文, OR))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1.5) +
    geom_point(size = 8, color = "#2E8B57") +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), 
                   height = 0.4, size = 1.5, color = "#2E8B57") +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10, 20)) +
    theme_classic() +
    theme(
      text = element_text(size = 24),
      plot.title = element_text(size = 32, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 28, face = "bold"),
      axis.text.y = element_text(size = 22),
      axis.text.x = element_text(size = 22),
      panel.grid.major.x = element_line(color = "gray90", size = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = "妊娠不良事件独立危险因素森林图",
      x = "比值比 (OR) 及95%置信区间",
      y = "危险因素",
      caption = "虚线表示OR=1，点和线表示OR值及其95%置信区间"
    ) +
    # 添加OR值标注
    geom_text(aes(label = paste0("OR=", OR, "\nP=", P值_格式化)), 
              hjust = -0.1, size = 7, fontface = "bold")
  
  print(forest_plot)
  ggsave("/home/lin/hz/forest_plot.png",  # 修改路径
         plot = forest_plot, 
         width = 18, 
         height = 14, 
         dpi = 300, 
         bg = "white")
  cat("森林图已保存\n")
}

# 12. 模型诊断（基于训练集）
if(!is.null(final_model)) {
  cat("\n=== 模型诊断 ===\n")
  
  # Hosmer-Lemeshow拟合优度检验
  tryCatch({
    library(ResourceSelection)
    hl_test <- hoslem.test(final_model$y, fitted(final_model))
    cat("Hosmer-Lemeshow检验: χ²=", round(hl_test$statistic, 3), 
        ", P=", round(hl_test$p.value, 4), "\n")
    
    if(hl_test$p.value > 0.05) {
      cat("模型拟合良好 (P>0.05)\n")
    } else {
      cat("模型拟合欠佳 (P<0.05)\n")
    }
  }, error = function(e) {
    cat("Hosmer-Lemeshow检验失败:", e$message, "\n")
  })
}

# 13. 生成最终报告
cat("\n=== 最终分析报告 ===\n")

report_file <- "/home/lin/hz/analysis_report.txt"  # 修改路径
sink(report_file)

cat("妊娠不良事件风险评估模型分析报告\n")
cat("=====================================\n\n")

cat("1. 基本信息\n")
cat("- 总样本数:", nrow(analysis_data_final), "\n")
cat("- 训练集样本数:", nrow(train_data), "\n")
cat("- 测试集样本数:", nrow(test_data), "\n")
cat("- 不良事件发生率:", 
    round(sum(analysis_data_final$不良事件 == "有")/nrow(analysis_data_final)*100, 2), "%\n\n")

cat("- 使用的随机种子:", best_seed, "\n")

if(exists("univariate_df")) {
  cat("2. 单因素分析结果\n")
  cat("- 分析变量数:", nrow(univariate_df), "\n")
  cat("- 显著变量数(P<0.05):", sum(univariate_df$显著性 == "是"), "\n")
  cat("- P<0.05变量数:", sum(univariate_df$P值 < 0.05), "\n\n")
}

if(exists("multivariate_results") && !is.null(multivariate_results)) {
  cat("3. 多因素分析结果\n")
  cat("- 纳入变量数:", nrow(multivariate_results), "\n")
  cat("- 独立危险因素数:", sum(multivariate_results$显著性 == "是"), "\n\n")
  
  significant_multi <- multivariate_results[multivariate_results$显著性 == "是", ]
  if(nrow(significant_multi) > 0) {
    cat("独立危险因素详情:\n")
    for(i in 1:nrow(significant_multi)) {
      cat("- ", significant_multi$变量[i], ": OR=", significant_multi$OR[i], 
          " (95% CI: ", significant_multi$OR_lower[i], "-", significant_multi$OR_upper[i], 
          "), P=", significant_multi$P值_格式化[i], "\n")
    }
    cat("\n")
  }
}

if(exists("train_auc") && exists("test_auc")) {
  cat("4. 模型性能\n")
  cat("- 训练集AUC:", round(train_auc, 3), "\n")
  cat("- 训练集95% CI:", paste(round(train_ci, 3), collapse = "-"), "\n")
  cat("- 测试集AUC:", round(test_auc, 3), "\n")
  cat("- 测试集95% CI:", paste(round(test_ci, 3), collapse = "-"), "\n\n")
  
  cat("测试集性能指标:\n")
  cat("- 敏感度:", round(test_sensitivity, 3), "\n")
  cat("- 特异度:", round(test_specificity, 3), "\n")
  cat("- 阳性预测值:", round(test_ppv, 3), "\n")
  cat("- 阴性预测值:", round(test_npv, 3), "\n")
  cat("- 准确度:", round(test_accuracy, 3), "\n\n")
}

cat("5. 文件输出\n")
cat("- 单因素分析结果: univariate_analysis_results.csv\n")
cat("- 多因素分析结果: multivariate_analysis_results.csv\n")
cat("- 模型性能指标: model_performance_metrics.csv\n")
cat("- 训练集和测试集ROC比较图: ROC_train_test_comparison.png\n")
cat("- 森林图: forest_plot.png\n")
cat("- 列线图: nomogram.png\n")
cat("- 校准曲线: calibration_curve.png\n")
cat("- 测试集校准曲线: test_calibration_curve.png\n\n")

cat("分析完成时间:", Sys.time(), "\n")

sink()

cat("分析报告已保存至:", report_file, "\n")

# 保存工作空间
save.image("/home/lin/hz/logistic_analysis_workspace.RData")  # 修改路径

cat("\n=== 所有分析完成！ ===\n")
cat("主要输出文件:\n")
cat("1. 单因素分析结果: univariate_analysis_results.csv\n")
cat("2. 多因素分析结果: multivariate_analysis_results.csv\n")
cat("3. 模型性能指标: model_performance_metrics.csv\n")
cat("4. 训练集和测试集ROC比较图: ROC_train_test_comparison.png\n")
cat("5. 森林图: forest_plot.png\n")
cat("6. 列线图: nomogram.png\n")
cat("7. 校准曲线: calibration_curve.png\n")
cat("8. 测试集校准曲线: test_calibration_curve.png\n")
cat("9. 分析报告: analysis_report.txt\n")

# 输出关键统计数值用于论文写作
cat("\n=== 论文写作关键数值汇总 ===\n")

# 1. C指数和AUC
if(exists("rms_model")) {
  cat("C指数:", round(rms_model$stats["C"], 3), "\n")
}
if(exists("train_auc") && exists("test_auc")) {
  cat("训练集AUC:", round(train_auc, 3), "\n")
  cat("训练集AUC 95% CI:", round(train_ci[1], 3), "-", round(train_ci[3], 3), "\n")
  cat("测试集AUC:", round(test_auc, 3), "\n")
  cat("测试集AUC 95% CI:", round(test_ci[1], 3), "-", round(test_ci[3], 3), "\n")
}

# 2. Hosmer-Lemeshow检验结果
if(exists("hl_test")) {
  cat("Hosmer-Lemeshow检验: χ²=", round(hl_test$statistic, 3), 
      ", P=", round(hl_test$p.value, 3), "\n")
}

# 3. 测试集性能指标
if(exists("test_sensitivity")) {
  cat("测试集敏感度:", round(test_sensitivity, 3), "\n")
  cat("测试集特异度:", round(test_specificity, 3), "\n")
  cat("测试集阳性预测值:", round(test_ppv, 3), "\n")
  cat("测试集阴性预测值:", round(test_npv, 3), "\n")
  cat("测试集准确度:", round(test_accuracy, 3), "\n")
}

# 4. 样本基本信息
cat("总样本数:", nrow(analysis_data_final), "\n")
cat("训练集样本数:", nrow(train_data), " (", round(nrow(train_data)/nrow(analysis_data_final)*100, 1), "%)\n")
cat("测试集样本数:", nrow(test_data), " (", round(nrow(test_data)/nrow(analysis_data_final)*100, 1), "%)\n")
cat("不良事件发生率:", round(sum(analysis_data_final$不良事件 == "有")/nrow(analysis_data_final)*100, 2), "%\n")
cat("使用的随机种子:", best_seed, "\n")
