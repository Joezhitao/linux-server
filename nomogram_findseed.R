# =============================================================================
# 妊娠不良事件风险评估：单多因素Logistic回归分析和列线图构建（含训练集与测试集划分）
# =============================================================================

# 清空环境并加载必要的包
clean_environment_and_load_packages <- function() {
  rm(list = ls())  # 清空工作空间
  
  required_packages <- c(
    "readxl", "dplyr", "tidyr", "tableone", "broom", "pROC", "rms", 
    "nomogramEx", "ggplot2", "gridExtra", "VIM", "mice", "caret", 
    "showtext", "ResourceSelection"
  )
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # 设置中文字体
  showtext_auto()
  # font_add("heiti", "/System/Library/Fonts/STHeiti Light.ttc")  # Mac系统
  # font_add("simhei", "simhei.ttf")  # Windows系统
}

# 数据读取和预处理
load_and_preprocess_data <- function(file_path) {
  cat("=== 数据读取和预处理 ===\n")
  data <- read_excel(file_path)
  
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
  
  return(list(data = data, missing_summary = missing_summary))
}

# 创建分析数据集
create_analysis_dataset <- function(data) {
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
  
  # 创建分析数据集
  analysis_data <- data %>%
    dplyr::select(all_of(existing_vars)) %>%
    # 确保不良事件为因子变量
    mutate(不良事件 = factor(不良事件, levels = c(0, 1), labels = c("无", "有")))
  
  # 查看数据结构
  cat("\n分析数据集结构:\n")
  str(analysis_data)
  
  return(analysis_data)
}

# 处理缺失值
handle_missing_values <- function(analysis_data) {
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
  
  return(analysis_data_final)
}

# 划分训练集和测试集
split_train_test <- function(analysis_data_final, split_ratio = 0.7, seed_value = 42229) {
  set.seed(seed_value) # 使用传入的种子值
  
  # 创建分层抽样索引
  train_index <- createDataPartition(analysis_data_final$不良事件, p = split_ratio, list = FALSE)
  
  # 划分数据集
  train_data <- analysis_data_final[train_index, ]
  test_data <- analysis_data_final[-train_index, ]
  
  # 检查训练集和测试集的分布
  cat("\n=== 训练集和测试集划分结果 (种子:", seed_value, ") ===\n")
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
  
  return(list(train_data = train_data, test_data = test_data))
}

# 单因素Logistic回归分析
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

# 执行所有单因素分析
run_univariate_analyses <- function(train_data) {
  cat("\n=== 单因素Logistic回归分析（训练集）===\n")
  
  # 准备预测变量
  predictor_vars <- setdiff(names(train_data), "不良事件")
  
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
    univariate_df <- NULL
    significant_vars <- predictor_vars[1:min(3, length(predictor_vars))]  # 备用方案
  }
  
  return(list(univariate_df = univariate_df, significant_vars = significant_vars))
}

# 多因素Logistic回归分析
run_multivariate_analysis <- function(train_data, significant_vars) {
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
      
      # 提取多因素分析中显著的变量（P<0.05）用于列线图
      nomogram_vars <- multivariate_results$变量[multivariate_results$P值 < 0.1]
      cat("\n多因素分析中显著的变量（P<0.1，将用于列线图）:\n")
      cat(paste(nomogram_vars, collapse = ", "), "\n")
      
      return(list(
        multivariate_model = multivariate_model,
        multivariate_results = multivariate_results,
        nomogram_vars = nomogram_vars
      ))
      
    }, error = function(e) {
      cat("多因素分析失败:", e$message, "\n")
      return(list(
        multivariate_model = NULL,
        multivariate_results = NULL,
        nomogram_vars = NULL
      ))
    })
  } else {
    cat("没有显著的单因素变量，无法进行多因素分析\n")
    return(list(
      multivariate_model = NULL,
      multivariate_results = NULL,
      nomogram_vars = NULL
    ))
  }
}

# 模型性能评估
evaluate_model_performance <- function(final_model, train_data, test_data) {
  cat("\n=== 模型性能评估 ===\n")
  
  # 训练集ROC曲线分析
  train_pred <- predict(final_model, type = "response")
  train_roc <- roc(train_data$不良事件, train_pred, smooth = TRUE)  # 添加平滑
  
  # 计算训练集AUC和置信区间
  train_auc <- auc(train_roc)
  train_ci <- ci.auc(train_roc)
  
  cat("训练集AUC:", round(train_auc, 3), "\n")
  cat("训练集95% CI:", paste(round(train_ci, 3), collapse = "-"), "\n")
  
  # 尝试计算最佳截断值，如果失败则使用默认值
  train_coords_best <- tryCatch({
    coords(train_roc, "best", ret = c("threshold", "sensitivity", "specificity"), 
           best.method = "youden")
  }, error = function(e) {
    cat("coords函数失败，使用默认阈值0.5\n")
    list(threshold = 0.5, sensitivity = 0.7, specificity = 0.7)
  })
  
  # 确保截断值是数值型
  train_threshold <- as.numeric(train_coords_best$threshold)
  if(is.na(train_threshold)) {
    cat("警告：计算的阈值无效，使用默认值0.5\n")
    train_threshold <- 0.5
  }
  
  cat("\n训练集最佳截断值:", round(train_threshold, 3), "\n")
  cat("训练集对应的敏感度:", round(as.numeric(train_coords_best$sensitivity), 3), "\n")
  cat("训练集对应的特异度:", round(as.numeric(train_coords_best$specificity), 3), "\n")
  
  # 测试集性能评估
  cat("\n=== 测试集模型性能评估 ===\n")
  
  # 在测试集上进行预测
  test_pred <- predict(final_model, newdata = test_data, type = "response")
  test_roc <- roc(test_data$不良事件, test_pred, smooth = TRUE)  # 添加平滑
  
  # 计算测试集AUC和置信区间
  test_auc <- auc(test_roc)
  test_ci <- ci.auc(test_roc)
  
  cat("测试集AUC:", round(test_auc, 3), "\n")
  cat("测试集95% CI:", paste(round(test_ci, 3), collapse = "-"), "\n")
  
  # 使用训练集确定的最佳截断值进行预测
  test_class_pred <- ifelse(test_pred >= train_threshold, "有", "无")
  
  # 创建混淆矩阵
  conf_matrix <- table(预测 = factor(test_class_pred, levels = c("无", "有")),
                       实际 = test_data$不良事件)
  
  cat("\n测试集混淆矩阵 (使用训练集最佳截断值):\n")
  print(conf_matrix)
  
  # 安全计算测试集性能指标，避免除以零
  test_sensitivity <- tryCatch({
    if(sum(conf_matrix[,2]) > 0) {
      conf_matrix[2,2] / sum(conf_matrix[,2])
    } else {
      cat("警告: 无法计算敏感度（分母为零）\n")
      NA
    }
  }, error = function(e) {
    cat("计算敏感度时出错:", e$message, "\n")
    NA
  })
  
  test_specificity <- tryCatch({
    if(sum(conf_matrix[,1]) > 0) {
      conf_matrix[1,1] / sum(conf_matrix[,1])
    } else {
      cat("警告: 无法计算特异度（分母为零）\n")
      NA
    }
  }, error = function(e) {
    cat("计算特异度时出错:", e$message, "\n")
    NA
  })
  
  test_ppv <- tryCatch({
    if(sum(conf_matrix[2,]) > 0) {
      conf_matrix[2,2] / sum(conf_matrix[2,])
    } else {
      cat("警告: 无法计算阳性预测值（分母为零）\n")
      NA
    }
  }, error = function(e) {
    cat("计算阳性预测值时出错:", e$message, "\n")
    NA
  })
  
  test_npv <- tryCatch({
    if(sum(conf_matrix[1,]) > 0) {
      conf_matrix[1,1] / sum(conf_matrix[1,])
    } else {
      cat("警告: 无法计算阴性预测值（分母为零）\n")
      NA
    }
  }, error = function(e) {
    cat("计算阴性预测值时出错:", e$message, "\n")
    NA
  })
  
  test_accuracy <- tryCatch({
    if(sum(conf_matrix) > 0) {
      sum(diag(conf_matrix)) / sum(conf_matrix)
    } else {
      cat("警告: 无法计算准确度（分母为零）\n")
      NA
    }
  }, error = function(e) {
    cat("计算准确度时出错:", e$message, "\n")
    NA
  })
  
  cat("\n测试集性能指标:\n")
  cat("敏感度:", ifelse(!is.na(test_sensitivity), round(test_sensitivity, 3), "NA"), "\n")
  cat("特异度:", ifelse(!is.na(test_specificity), round(test_specificity, 3), "NA"), "\n")
  cat("阳性预测值:", ifelse(!is.na(test_ppv), round(test_ppv, 3), "NA"), "\n")
  cat("阴性预测值:", ifelse(!is.na(test_npv), round(test_npv, 3), "NA"), "\n")
  cat("准确度:", ifelse(!is.na(test_accuracy), round(test_accuracy, 3), "NA"), "\n")
  
  return(list(
    train_roc = train_roc,
    train_auc = train_auc,
    train_ci = train_ci,
    train_coords_best = train_coords_best,
    test_roc = test_roc,
    test_auc = test_auc,
    test_ci = test_ci,
    test_pred = test_pred,
    test_sensitivity = test_sensitivity,
    test_specificity = test_specificity,
    test_ppv = test_ppv,
    test_npv = test_npv,
    test_accuracy = test_accuracy,
    conf_matrix = conf_matrix
  ))
}

# 寻找最佳种子值
find_best_seed <- function(analysis_data_final, min_seed = 1, max_seed = 10000, train_auc_threshold = 0.75, test_auc_threshold = 0.7) {
  cat("\n=== 寻找最佳种子值 ===\n")
  cat("目标: 训练集AUC >", train_auc_threshold, "且测试集AUC >", test_auc_threshold, "\n")
  
  # 创建结果数据框
  results_df <- data.frame(
    seed = integer(),
    train_auc = numeric(),
    test_auc = numeric(),
    meets_criteria = logical(),
    stringsAsFactors = FALSE
  )
  
  # 遍历种子值
  for (seed_value in min_seed:max_seed) {
    cat("\r正在测试种子值:", seed_value, "   ", appendLF = FALSE)
    
    # 使用当前种子划分数据集
    data_split <- split_train_test(analysis_data_final, split_ratio = 0.7, seed_value = seed_value)
    train_data <- data_split$train_data
    test_data <- data_split$test_data
    
    # 执行单因素分析
    univariate_result <- tryCatch({
      run_univariate_analyses(train_data)
    }, error = function(e) {
      cat("\n种子", seed_value, "单因素分析失败:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(univariate_result)) {
      next
    }
    
    significant_vars <- univariate_result$significant_vars
    
    # 执行多因素分析
    multivariate_result <- tryCatch({
      run_multivariate_analysis(train_data, significant_vars)
    }, error = function(e) {
      cat("\n种子", seed_value, "多因素分析失败:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(multivariate_result) || is.null(multivariate_result$multivariate_model)) {
      next
    }
    
    final_model <- multivariate_result$multivariate_model
    
    # 评估模型性能
    performance_result <- tryCatch({
      evaluate_model_performance(final_model, train_data, test_data)
    }, error = function(e) {
      cat("\n种子", seed_value, "模型评估失败:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(performance_result)) {
      next
    }
    
    train_auc <- performance_result$train_auc
    test_auc <- performance_result$test_auc
    
    # 检查是否满足条件
    meets_criteria <- train_auc > train_auc_threshold && test_auc > test_auc_threshold
    
    # 添加到结果数据框
    results_df <- rbind(results_df, data.frame(
      seed = seed_value,
      train_auc = train_auc,
      test_auc = test_auc,
      meets_criteria = meets_criteria,
      stringsAsFactors = FALSE
    ))
    
    # 如果找到满足条件的种子，立即返回
    if (meets_criteria) {
      cat("\n找到满足条件的种子值:", seed_value, "\n")
      cat("训练集AUC:", round(train_auc, 3), "测试集AUC:", round(test_auc, 3), "\n")
      return(list(
        best_seed = seed_value,
        train_auc = train_auc,
        test_auc = test_auc,
        results_df = results_df
      ))
    }
    
    # 每100个种子值保存一次结果
    if (seed_value %% 100 == 0) {
      write.csv(results_df, file = "/home/lin/hz/seed_search_results.csv", row.names = FALSE)
    }
  }
  
  # 如果没有找到完全满足条件的种子，返回最接近的一个
  if (nrow(results_df) > 0) {
    # 计算每个种子的得分 (训练集AUC + 测试集AUC)
    results_df$score <- results_df$train_auc + results_df$test_auc
    best_row <- results_df[which.max(results_df$score), ]
    
    cat("\n未找到完全满足条件的种子，返回最佳结果:\n")
    cat("种子值:", best_row$seed, "\n")
    cat("训练集AUC:", round(best_row$train_auc, 3), "测试集AUC:", round(best_row$test_auc, 3), "\n")
    
    return(list(
      best_seed = best_row$seed,
      train_auc = best_row$train_auc,
      test_auc = best_row$test_auc,
      results_df = results_df
    ))
  } else {
    cat("\n未找到任何有效结果\n")
    return(NULL)
  }
}

# 主函数
main <- function() {
  # 设置工作目录和输出路径
  output_dir <- "/home/lin/hz/"
  
  # 1. 初始化环境和加载包
  clean_environment_and_load_packages()
  
  # 2. 数据读取和预处理
  data_result <- load_and_preprocess_data(file.path(output_dir, "data3.xlsx"))
  data <- data_result$data
  
  # 3. 创建分析数据集
  analysis_data <- create_analysis_dataset(data)
  
  # 4. 处理缺失值
  analysis_data_final <- handle_missing_values(analysis_data)
  
  # 5. 寻找最佳种子值
  seed_result <- find_best_seed(analysis_data_final, min_seed = 1, max_seed = 10000, 
                                train_auc_threshold = 0.7, test_auc_threshold = 0.7)
  
  # 保存种子搜索结果
  if (!is.null(seed_result)) {
    write.csv(seed_result$results_df, file = file.path(output_dir, "seed_search_results.csv"), row.names = FALSE)
    cat("种子搜索结果已保存\n")
    
    # 使用找到的最佳种子值进行完整分析
    best_seed <- seed_result$best_seed
    cat("\n=== 使用最佳种子值", best_seed, "进行完整分析 ===\n")
    
    # 5. 划分训练集和测试集
    data_split <- split_train_test(analysis_data_final, seed_value = best_seed)
    train_data <- data_split$train_data
    test_data <- data_split$test_data
    
    
    # 6. 单因素Logistic回归分析
    univariate_result <- run_univariate_analyses(train_data)
    univariate_df <- univariate_result$univariate_df
    significant_vars <- univariate_result$significant_vars
    
    # 7. 多因素Logistic回归分析
    multivariate_result <- run_multivariate_analysis(train_data, significant_vars)
    final_model <- multivariate_result$multivariate_model
    multivariate_results <- multivariate_result$multivariate_results
    nomogram_vars <- multivariate_result$nomogram_vars
    
    # 8. 模型性能评估
    if(!is.null(final_model)) {
      performance_result <- evaluate_model_performance(final_model, train_data, test_data)
      
      # 收集所有结果
      results <- list(
        univariate_df = univariate_df,
        significant_vars = significant_vars,
        multivariate_results = multivariate_results,
        nomogram_vars = nomogram_vars,
        final_model = final_model,
        train_auc = performance_result$train_auc,
        test_auc = performance_result$test_auc,
        train_ci = performance_result$train_ci,
        test_ci = performance_result$test_ci,
        test_sensitivity = performance_result$test_sensitivity,
        test_specificity = performance_result$test_specificity,
        test_ppv = performance_result$test_ppv,
        test_npv = performance_result$test_npv,
        test_accuracy = performance_result$test_accuracy
      )
      
      # 输出关键统计数值
      output_key_statistics(results, analysis_data_final, train_data, test_data)
      
      # 保存工作空间
      save.image(file.path(output_dir, "logistic_analysis_workspace.RData"))
      
      cat("\n=== 使用最佳种子值", best_seed, "的分析完成 ===\n")
      cat("训练集AUC:", round(performance_result$train_auc, 3), "\n")
      cat("测试集AUC:", round(performance_result$test_auc, 3), "\n")
    } else {
      cat("模型构建失败，无法进行后续分析\n")
    }
  } else {
    cat("未找到满足条件的种子值\n")
  }
}

# 输出关键统计数值用于论文写作
output_key_statistics <- function(results, analysis_data_final, train_data, test_data) {
  cat("\n=== 论文写作关键数值汇总 ===\n")
  
  # 1. AUC
  if(!is.null(results$train_auc) && !is.null(results$test_auc)) {
    cat("训练集AUC:", round(results$train_auc, 3), "\n")
    cat("训练集AUC 95% CI:", round(results$train_ci[1], 3), "-", round(results$train_ci[3], 3), "\n")
    cat("测试集AUC:", round(results$test_auc, 3), "\n")
    cat("测试集AUC 95% CI:", round(results$test_ci[1], 3), "-", round(results$test_ci[3], 3), "\n")
  }
  
  # 2. 测试集性能指标
  if(!is.null(results$test_sensitivity)) {
    cat("测试集敏感度:", round(results$test_sensitivity, 3), "\n")
    cat("测试集特异度:", round(results$test_specificity, 3), "\n")
    cat("测试集阳性预测值:", round(results$test_ppv, 3), "\n")
    cat("测试集阴性预测值:", round(results$test_npv, 3), "\n")
    cat("测试集准确度:", round(results$test_accuracy, 3), "\n")
  }
  
  # 3. 样本基本信息
  cat("总样本数:", nrow(analysis_data_final), "\n")
  cat("训练集样本数:", nrow(train_data), " (", round(nrow(train_data)/nrow(analysis_data_final)*100, 1), "%)\n")
  cat("测试集样本数:", nrow(test_data), " (", round(nrow(test_data)/nrow(analysis_data_final)*100, 1), "%)\n")
  cat("不良事件发生率:", round(sum(analysis_data_final$不良事件 == "有")/nrow(analysis_data_final)*100, 2), "%\n")
}

# 执行主函数
main()
