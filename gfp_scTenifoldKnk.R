library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)  # 用于防止标签重叠

setwd("/home/lin/c_group/")
pbmc <- readRDS('/home/lin/c_group/HSCs_ident.rds')

# 提取表达矩阵
countMatrix <- GetAssayData(pbmc, layer = "counts")

# 检查目标基因是否在行名里
"Tgfb1" %in% rownames(countMatrix)  # 应该是TRUE


#选择敲除的基因进行计算,敲除基因需要从之前高变基因前1000内的进行。
knk_gene <- "Tgfb1"

result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = knk_gene, #需要敲除的基因
                        qc = TRUE,#是否进行QC
                        qc_mtThreshold = 0.1,#mt阈值
                        qc_minLSize = 2000,#文库阈值(细胞测到的基因总数)
                        nc_nNet = 10, #子网络数量
                        nc_nCells = 500, #每个网络中随机抽取的细胞数
                        nc_nComp = 3 #PCA 的主成分数量
)

top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 20)
ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes",
       x="Gene", y="FC") +
  theme_minimal()

df <- result$diffRegulation
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.01)
ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value)",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()
