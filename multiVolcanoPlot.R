color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


multiVolcanoPlot = function(dat, color.arr=NULL, onlyAnnotateUp=T,
                            log2Foldchang=0.58, adjp=0.05, top_marker=15, 
                            max_overlaps=10, width=0.9){
  library(dplyr)
  library(ggrepel)
  # set default color list
  if(is.null(color.arr)){
    len = length(unique(dat$cluster))
    color.arr=scales::hue_pal()(len)
  }
  
  dat.plot <- dat %>% mutate(
    "significance"=case_when(p_val_adj < adjp & avg_log2FC >= log2Foldchang  ~ 'Up',
                             p_val_adj < adjp & avg_log2FC <= -log2Foldchang  ~ 'Down',
                             TRUE ~ 'None'))
  tbl = table(dat.plot$significance)
  print( tbl )
  background.dat <- data.frame(
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC>0) %>%
      summarise("y.localup"=max(avg_log2FC)),
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC<=0) %>%
      summarise("y.localdown"=min(avg_log2FC)),
    x.local=seq(1:length(unique(dat.plot$cluster)))
  ) %>% select(-cluster.1)
  
  x.number <- background.dat %>% select(cluster, x.local)
  dat.plot <- dat.plot%>% left_join(x.number,by = "cluster")
  
  #selecting top-up and top-down proteins
  dat.marked.up <- dat.plot %>% filter(significance=="Up") %>%
    group_by(cluster) %>% arrange(-avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked.down <- dat.plot %>% filter(significance=="Down") %>%
    group_by(cluster) %>% arrange(avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked <- dat.marked.up %>% bind_rows(dat.marked.down)
  
  # 增加框的高度
  tile_height = log2Foldchang * 3.0
  
  # 调整信息行的位置，将细胞簇标签下移
  dat.infor <- background.dat %>%
    mutate("y.infor"=rep(-0.5, length(cluster)))  # 将标签下移到框的底部
  
  ##plotting:
  vol.plot <- ggplot()+
    # background
    geom_col(background.dat,mapping=aes(x.local, y.localup),
             fill="grey80", alpha=0.2, width=0.95, just = 0.5)+
    geom_col(background.dat,mapping=aes(x.local,y.localdown),
             fill="grey80", alpha=0.2, width=0.95, just = 0.5)+
    # point plot
    geom_jitter(dat.plot, mapping=aes(x.local, avg_log2FC,
                                      color=significance),
                size=1.2, width = 0.4, alpha= 1)+
    scale_color_manual(name="significance", 
                       breaks = c('Up', 'None', 'Down'),
                       values = c("#d56e5e","#cccccc", "#5390b5")) +
    # 增大细胞簇框的高度和宽度
    geom_tile(dat.infor, mapping=aes(x.local, y.infor),
              height = tile_height,
              fill = color.arr[1:length(unique(dat.plot$cluster))],
              alpha = 0.5,
              width = width * 1.1) +
    labs(x=NULL,y="log2 Fold change")+
    # 增大细胞类型标签字体并添加背景框，使其更突出
    geom_label(dat.infor, mapping=aes(x.local, y=-1.5, label=cluster),  # 将标签位置下移
               size=8, fontface="bold", 
               fill = "white", color = "black",
               label.size = 0.5,  # 添加边框
               label.padding = unit(0.5, "lines")) +  # 增加内边距
    # 基因标签
    ggrepel::geom_label_repel(data=if(onlyAnnotateUp) dat.marked.up else dat.marked,
                              mapping=aes(x=x.local, y=avg_log2FC, label=gene),
                              size=5.5,
                              force = 4,  # 增加排斥力
                              max.overlaps = 25,
                              label.size = 0.2,
                              fill = "white",
                              alpha = 0.85,
                              seed = 123,
                              min.segment.length = 0,
                              point.padding = 0.5,
                              box.padding = 0.7,  # 增加标签间距
                              segment.size = 0.5,
                              segment.linetype = 1)+
    # 注释文本
    annotate("text", x=1.5, y=max(background.dat$y.localup)+1,
             label=paste0("|log2FC|>=", log2Foldchang, " & FDR<", adjp),
             size=6, fontface="bold")+
    # 主题设置
    theme_classic(base_size = 18)+
    
    theme(
      axis.title = element_text(size = 22, color = "black", face="bold"),
      axis.text = element_text(size = 20, color = "black"),
      axis.line.y = element_line(color = "black", size = 1.2),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_text(size=20, face="bold"),
      legend.text = element_text(size=18),
      legend.spacing.x = unit(0.2,'cm'),
      legend.key.width = unit(1.0,'cm'),
      legend.key.height = unit(1.0,'cm'),
      legend.background = element_blank(),
      legend.box = "horizontal",
      legend.position = c(0.9, 0.9),
      legend.justification = c(1, 1)
    )+
    guides(
      color=guide_legend(override.aes = list(size=10), title="Change")
    )
  vol.plot
}


group <- c("RT")


for (i in group) {
  # 更改文件路径和扩展名为CSV
  file = paste("/home/lin/gfp/火山图/", i, "_combined.csv", sep = "")
  
  # 使用read.csv代替read.xlsx
  DEG <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # 其余代码保持不变
  path = paste("/home/lin/gfp/火山图/Volcano", i, ".png", sep = "")
  plot <- multiVolcanoPlot(DEG, top_marker = 10, color.arr = color.pals, onlyAnnotateUp = F)
  png(path, units = "in", width = 26, height = 13, res = 800)
  print(plot)
  dev.off()
}


