library(dplyr)
library(openxlsx)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(data.table)

#Read in raw data and run doupletFinder
skin.fetus <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/skin.fetus.rds")


#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
suppressMessages(require(DoubletFinder))

skin.fetus = FindVariableFeatures(skin.fetus, verbose = F)
skin.fetus = ScaleData(skin.fetus, vars.to.regress = c("nFeature_RNA", "percent.mito"),
                      verbose = F)
skin.fetus = RunPCA(skin.fetus, verbose = F, npcs = 20)
skin.fetus = RunUMAP(skin.fetus, dims = 1:10, verbose = F)

# Re-visualize the clusters separated by orig.ident
DimPlot(object = skin.fetus, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)


# Can run parameter optimization with paramSweep (here one could possibly adjust params more and play around with them...)
sweep.res <- paramSweep_v3(skin.fetus)
sweep.stats <-summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cellscells.
fetus.nExp <- round(ncol(skin.fetus) * 0.04)  # expect 4% doublets
skin.fetus <- doubletFinder_v3(skin.fetus, pN = 0.25, pK = 0.09, nExp = fetus.nExp, PCs = 1:10)

# name of the DF prediction can change
skin.fetus.DF.name = colnames(skin.fetus@meta.data)[grepl("DF.classification", colnames(skin.fetus@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(skin.fetus, group.by = "seurat_clusters") + NoAxes(),
                   DimPlot(skin.fetus, group.by = skin.fetus.DF.name) + NoAxes())
