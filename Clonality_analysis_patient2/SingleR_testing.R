library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

table(hpca.se$label.main)
table(hpca.se$label.fine)
table(hpca.se$label.ont)
table(hpca.se$labels)


rownames(hpca.se[, hpca.se$label.fine == "T_cell:Treg:Naive"])
colData(hpca.se)
rownames(hpca.se)
rowRanges(hpca.se)
