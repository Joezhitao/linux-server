# 将以下内容保存到 /home/lin/R_linux_server/linux_server/jjVolcano_mod(魔改火山图).R

jjVolcano_mod <- function (diffData = NULL, myMarkers = NULL, order.by = c("avg_log2FC"), 
                           log2FC.cutoff = 0.25, pvalue.cutoff = 0.05, adjustP.cutoff = 0.01, 
                           topGeneN = 5, col.type = "updown", back.col = "grey93", 
                           pSize = 0.75, aesCol = c("#0099CC", "#CC3333"), legend.position = c(0.7, 
                                                                                               0.9), base_size = 14, tile.col = jjAnno::useMyCol("paired", 
                                                                                                                                                 n = 9), cluster.order = NULL, polar = FALSE, expand = c(-1, 
                                                                                                                                                                                                         1), flip = FALSE, celltypeSize = 3, label_box = TRUE, y_limits = NULL, reverse_y = FALSE, ...) 
{
  diff.marker <- diffData %>% dplyr::filter(abs(avg_log2FC) >= 
                                              log2FC.cutoff & p_val < pvalue.cutoff)
  
  # 如果需要反转y轴，反转avg_log2FC的符号
  if(reverse_y) {
    diff.marker$avg_log2FC <- -diff.marker$avg_log2FC
  }
  
  diff.marker <- diff.marker %>% dplyr::mutate(type = ifelse(avg_log2FC >= 
                                                               log2FC.cutoff, "sigUp", "sigDown")) %>% dplyr::mutate(type2 = ifelse(p_val_adj < 
                                                                                                                                      adjustP.cutoff, paste("adjust Pvalue < ", adjustP.cutoff, 
                                                                                                                                                            sep = ""), paste("adjust Pvalue >= ", adjustP.cutoff, 
                                                                                                                                                                             sep = "")))
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
  }
  back.data <- purrr::map_df(unique(diff.marker$cluster), 
                             function(x) {
                               tmp <- diff.marker %>% dplyr::filter(cluster == 
                                                                      x)
                               new.tmp <- data.frame(cluster = x, min = min(tmp$avg_log2FC) - 
                                                       0.2, max = max(tmp$avg_log2FC) + 0.2)
                               return(new.tmp)
                             })
  top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
  top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = topGeneN, 
                                                        order_by = get(order.by))
  top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = topGeneN, 
                                                        order_by = get(order.by))
  top.marker <- rbind(top.marker.max, top.marker.min)
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>% dplyr::filter(gene %in% 
                                                  myMarkers)
  }
  else {
    top.marker <- top.marker
  }
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = cluster, 
                                                  y = avg_log2FC)) + ggplot2::geom_col(data = back.data, 
                                                                                       ggplot2::aes(x = cluster, y = min), fill = back.col) + 
    ggplot2::geom_col(data = back.data, ggplot2::aes(x = cluster, 
                                                     y = max), fill = back.col)
  if (col.type == "updown") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type), 
                                    size = pSize) + ggplot2::scale_color_manual(values = c(sigDown = aesCol[1], 
                                                                                           sigUp = aesCol[2]))
  }
  else if (col.type == "adjustP") {
    p2 <- p1 + ggplot2::geom_jitter(ggplot2::aes(color = type2), 
                                    size = pSize) + ggplot2::scale_color_manual(values = c(aesCol[2], 
                                                                                           aesCol[1]))
  }
  
  # 更新legend.position的使用方式
  if (is.numeric(legend.position) && length(legend.position) == 2) {
    # 对于ggplot2 3.5.0及以上版本
    legend_theme <- ggplot2::theme(
      legend.position.inside = legend.position,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
  } else {
    # 对于较旧版本或非数值位置（如"right", "none"等）
    legend_theme <- ggplot2::theme(
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
  }
  
  p3 <- p2 + ggplot2::scale_y_continuous(n.breaks = 6) + 
    ggplot2::theme_classic(base_size = base_size) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank()) + 
    legend_theme + 
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") + 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  
  # 首先添加tile和文本标签
  p4 <- p3 + ggplot2::geom_tile(ggplot2::aes(x = cluster, 
                                             y = 0, fill = cluster), color = "black", height = log2FC.cutoff * 
                                  2, alpha = 0.3, show.legend = FALSE) + ggplot2::scale_fill_manual(values = tile.col) + 
    ggrepel::geom_text_repel(data = top.marker, ggplot2::aes(x = cluster, 
                                                             y = avg_log2FC, label = gene), max.overlaps = 50, 
                             ...)
  
  # 设置y轴范围（如果提供）
  if (!is.null(y_limits)) {
    p4 <- p4 + ggplot2::ylim(y_limits)
  }
  
  if (polar == TRUE) {
    # 基本极坐标图，使用expansion而不是重新设置scale_y_continuous
    p5_base <- p4 + 
      ggplot2::scale_y_continuous(n.breaks = 6, 
                                  expand = ggplot2::expansion(mult = expand)) + 
      ggplot2::theme_void(base_size = base_size) + 
      legend_theme + 
      ggplot2::coord_polar(clip = "off", theta = "x")
    
    # 在最上层添加带框的cluster标签
    if (label_box) {
      # 创建一个数据框用于标签
      label_data <- data.frame(
        cluster = unique(diff.marker$cluster),
        y = 0  # 在y=0处放置标签
      )
      
      # 添加带框的标签在最上层
      p5 <- p5_base + 
        ggplot2::geom_label(data = label_data, 
                            ggplot2::aes(x = cluster, y = 0, label = cluster),
                            size = celltypeSize,
                            fontface = "bold",
                            fill = "white",
                            alpha = 0.8,
                            label.size = 0.5)  # 添加边框
    } else {
      # 使用原始的文本路径
      p5 <- p5_base + geomtextpath::geom_textpath(ggplot2::aes(x = cluster, 
                                                               y = 0, label = cluster))
    }
  }
  else {
    if (flip == TRUE) {
      if (label_box) {
        p5 <- p4 + 
          ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                           y = 0, label = cluster), 
                              size = celltypeSize,
                              fontface = "bold",
                              fill = "white",
                              alpha = 0.8,
                              label.size = 0.5) + 
          ggplot2::theme(axis.line.y = ggplot2::element_blank(), 
                         axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
          ggplot2::coord_flip()
      } else {
        p5 <- p4 + 
          ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                           y = 0, label = cluster), size = celltypeSize) + 
          ggplot2::theme(axis.line.y = ggplot2::element_blank(), 
                         axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + 
          ggplot2::coord_flip()
      }
    }
    else {
      if (label_box) {
        p5 <- p4 + 
          ggplot2::geom_label(ggplot2::aes(x = cluster, 
                                           y = 0, label = cluster), 
                              size = celltypeSize,
                              fontface = "bold",
                              fill = "white",
                              alpha = 0.8,
                              label.size = 0.5) + 
          ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                         axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
      } else {
        p5 <- p4 + 
          ggplot2::geom_text(ggplot2::aes(x = cluster, 
                                          y = 0, label = cluster), size = celltypeSize) + 
          ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                         axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
      }
    }
  }
  return(p5)
}
