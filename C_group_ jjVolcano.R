library(scRNAtoolVis)
library(readxl)

# test
data('pbmc.markers')

data <- read_excel("/home/lin/c_group/subpopulation_analysis/Hepatocytes/DiffGenes_LR_comparison/All_LR_comparisons.xlsx")
# expand limits
jjVolcano(diffData = data,
          tile.col = corrplot::COL2('RdYlBu', 15)[4:12],
          size  = 3.5,
          topGeneN = 10,
          fontface = 'italic',
          polar = T) +
  ylim(-8,10)
