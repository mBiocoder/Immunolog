#=================================================================================================================#
# Exploratory analysis of fetal human skin
# Author: Mahima Arunkumar
# Date: 28.01.2022 
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Data for fetal human skin from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156972
# Rscript purpose:
# - Identify the alpha and beta T-cells in the skin 
# - Compare this to healthy donor (EX0004)
# - DEG in order to identify a signature for "naive" T cells from the author's data
# (of the paper above)
# - Take this signature to identify "naive" T cells in our healthy donor (EX0004) and allo-HSCT donor
# - Check for developmental pathway using velocity and diffusion map tool
# - Clonality analysis to assess the expansion of the cells
# - Test for possible ligands
# - UMAP and violin plots to check for enrichment of residency signatures and the markers CD69, CD103 etc. 
# in the fetal T cells and compare to adult
# Goals and Background:
# - We would like to see whether naive T cells can persist/be resident in the adult skin
# - We want to know if they might even be remnants from fetal live (persisting precursor cells)
# - We want to see if there is a developmental pathway as there are also some memory T cells in the fetal skin
# - Only cells with the same TCR belong to the same family that might have differentiated. Expectation:
# Naive T cells should not be clonal, which is why this analysis is relevant to see what could have stimulated 
# them in utero
# - We want to see if certain residency signatures are enriched
#=================================================================================================================#

library(dplyr)
library(openxlsx)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)

####################### Load the PBMC dataset ############################################

pbmc.data <- Read10X_h5("../data/Paper/GSE156972_raw_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "fetal", min.cells = 3, min.features = 200)
pbmc

##################### Define experiment specific variables ###############################

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('sc_fetal_skin')
an.descs <- c('sc_fetal')

# Create new folder to store figures from this analysis
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/qc'))
#dir.create(paste0('figures/', analysisid, '_', an.desc, '/pca'))
dir.figs <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs)
dir.figs.qc <- paste0('figures/', analysisid, '_', an.desc, '/qc/', analysisid, '_', an.descs)
#dir.figs.pca <- paste0('figures/', analysisid, '_', an.desc, '/pca/', analysisid, '_', an.descs)


###################### QC and selecting cells for further analysis #######################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

############################### Normalizing the data #######################################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)

saveRDS(pbmc, file = "./output/fetal_result.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Find clusters
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster1.markers <- FindMarkers(pbmc, ident.2 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster2.markers <- FindMarkers(pbmc, ident.3 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster3.markers <- FindMarkers(pbmc, ident.4 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster4.markers <- FindMarkers(pbmc, ident.5 = 4, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster5.markers <- FindMarkers(pbmc, ident.6 = 5, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster6.markers <- FindMarkers(pbmc, ident.7 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster7.markers <- FindMarkers(pbmc, ident.8 = 7, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster8.markers <- FindMarkers(pbmc, ident.9 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster9.markers <- FindMarkers(pbmc, ident.10 = 9, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster10.markers <- FindMarkers(pbmc, ident.11 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster11.markers <- FindMarkers(pbmc, ident.12 = 11, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster12.markers <- FindMarkers(pbmc, ident.13 = 12, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster13.markers <- FindMarkers(pbmc, ident.14 = 13, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster14.markers <- FindMarkers(pbmc, ident.15 = 14, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster15.markers <- FindMarkers(pbmc, ident.16 = 15, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


#trm marker genes -> CD103 not found in our data
VlnPlot(pbmc, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
FeaturePlot(pbmc, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
RidgePlot(pbmc, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))

#separate plots for tcf7 gene
VlnPlot(pbmc, features = c("TCF7"))
FeaturePlot(pbmc, features = c("TCF7"))
RidgePlot(pbmc, features = c("TCF7"))

#check if dataset contains only T cells or all cells
VlnPlot(pbmc, features = c("GNLY", "GZMB", "MS4A1"))


#Heatmap of top 10
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)


################################## Adult data set EX004 #######################################

################################## Load the PBMC dataset ############################################

adult.data <- Read10X("../data/EX0004/raw_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
adult <- CreateSeuratObject(counts = adult.data, project = "adult", min.cells = 3, min.features = 200)
adult

# Define Experiment_ID according to ngs_sample_list
experimentid <- c("EX0004")

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx("../../ngs_sample_list.xlsx", sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = "NA", fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], summarise, sum = NA)

# Sample data
metadata <- data.frame(row.names = nsl$Sample_name, 
                       reshape2::colsplit(nsl$Sample_name, "_", c("exp", "celltype", "source")))


#Reference sample name
tr <- c("skin")


##################### Define experiment specific variables ###############################

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])


# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('sc_adult_skin')
an.descs <- c('sc_adult')

# Create new folder to store figures from this analysis
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/qc'))
#dir.create(paste0('figures/', analysisid, '_', an.desc, '/pca'))
dir.figs <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs)
dir.figs.qc <- paste0('figures/', analysisid, '_', an.desc, '/qc/', analysisid, '_', an.descs)
#dir.figs.pca <- paste0('figures/', analysisid, '_', an.desc, '/pca/', analysisid, '_', an.descs)

###################### QC and selecting cells for further analysis #######################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
adult[["percent.mt"]] <- PercentageFeatureSet(adult, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(adult, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adult, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


adult <- subset(adult, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

############################### Normalizing the data #######################################

adult <- NormalizeData(adult, normalization.method = "LogNormalize", scale.factor = 10000)
adult <- NormalizeData(adult)

#Identification of highly variable features (feature selection)
adult <- FindVariableFeatures(adult, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adult), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(adult)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data

all.genes <- rownames(adult)
adult <- ScaleData(adult, features = all.genes)


#Perform linear dimensional reduction
adult <- RunPCA(adult, features = VariableFeatures(object = adult))

# Examine and visualize PCA results a few different ways
print(adult[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(adult, dims = 1:2, reduction = "pca")
DimPlot(adult, reduction = "pca")
DimHeatmap(adult, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(adult, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
adult <- JackStraw(adult, num.replicate = 100)
adult <- ScoreJackStraw(adult, dims = 1:20)

JackStrawPlot(adult, dims = 1:15)
ElbowPlot(adult)

#Cluster the cells
adult <- FindNeighbors(adult, dims = 1:10)
adult <- FindClusters(adult, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(adult), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
adult <- RunUMAP(adult, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(adult, reduction = "umap", label = TRUE)

saveRDS(adult, file = "./output/adult_results.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(adult, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
adult.markers <- FindAllMarkers(adult, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
adult.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Find clusters
cluster0.markers <- FindMarkers(adult, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster1.markers <- FindMarkers(adult, ident.2 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster2.markers <- FindMarkers(adult, ident.3 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster3.markers <- FindMarkers(adult, ident.4 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster4.markers <- FindMarkers(adult, ident.5 = 4, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster5.markers <- FindMarkers(adult, ident.6 = 5, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster6.markers <- FindMarkers(adult, ident.7 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster7.markers <- FindMarkers(adult, ident.8 = 7, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster8.markers <- FindMarkers(adult, ident.9 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster9.markers <- FindMarkers(adult, ident.10 = 9, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


#trm marker genes -> CD103 not found in our data
VlnPlot(adult, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
FeaturePlot(adult, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
RidgePlot(adult, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))

#tcf7 gene separately
VlnPlot(adult, features = c("TCF7"))
FeaturePlot(adult, features = c("TCF7"))
RidgePlot(adult, features = c("TCF7"))


#Heatmap of top 10
adult.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(adult, features = top10$gene) + NoLegend()



################################### Integrated analysis adult and fetus ########################

#load rds file to object
adult <- readRDS("./output/adult_qc.rds")
skin <- readRDS("./output/fetal_qc.rds")

ifnb.list <- list(pbmc, adult)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)


#perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

#Integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")


#Identify conserved cell type markers
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers)

#Plots marker genes
VlnPlot(immune.combined, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
FeaturePlot(immune.combined, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
RidgePlot(immune.combined, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))

#tcf7 gene separately
VlnPlot(immune.combined, features = c("TCF7"))
FeaturePlot(immune.combined, features = c("TCF7"))
RidgePlot(immune.combined, features = c("TCF7"))



#Identify differential expressed genes across conditions
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)


FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


############################################# Like Gustavo ####################################
#Fetal
train_single <- trainSingleR(GetAssayData(pbmc, slot = 'data'),pbmc@active.ident, genes = 'de',  sd.thresh = 1,
                             de.method = c("classic", "wilcox", "t"),
                             de.n = NULL,
                             de.args = list(),
                             aggr.ref = FALSE,
                             aggr.args = list(),
                             recompute = TRUE,
                             restrict = NULL,
                             assay.type = "logcounts",
                             check.missing = TRUE  )

sr.bd <- classifySingleR(
  GetAssayData(pbmc, slot = 'data'),
  train_single,
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  sd.thresh = NULL,
  prune = TRUE,
  assay.type = "logcounts",
  check.missing = TRUE)

saveRDS(sr.bd, file = "fetal_singleR_results.rds")
# Load Seurat object
sk <- readRDS("./output/fetal_qc.rds")

sr.sk <- readRDS('fetal_singleR_results.rds')


# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
sr.sk$meta.data$orig.ident <- sk@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.sk$meta.data$xy <- sk@reductions$umap@cell.embeddings[,1] # the tSNE coordinates (Embeddings(sk, reduction = "tsne"))
sr.sk$meta.data$xy.um <- sk@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "umap"))
sr.sk$meta.data$clusters <- sk@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(sk))


# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
#sr.bd$seurat = seurat.object # (optional)
sr.bd$meta.data$orig.ident <- bd@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.bd$meta.data$xy <- bd@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = "tsne"))
sr.bd$meta.data$xy.um <- bd@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = "umap"))
sr.bd$meta.data$clusters <- bd@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(bd))
sr.bd$meta.data$celltype <- sr.bd$singler[[2]]$SingleR.single$labels
sr.bd$singler[[2]]$SingleR.single$labels <- plyr::mapvalues(sr.bd$singler[[2]]$SingleR.single$labels, 
                                                            from = levels(factor(sr.bd$singler[[2]]$SingleR.single$labels)), 
                                                            to = c("CD4 Tn", "CD4 Tcm", "CD4 Tem", 
                                                                   "CD8 Tn", "CD8 Tcm", "CD8 Tem", 
                                                                   "mB cell", "other", "mB cell", "other", "nB cell", "NK cell", "Treg"))


##################################### New try: adult ###############################################

pd <- list()


# Load Seurat object from healthy adult
pd[['tp.sk']] <- readRDS('./GA_AN0146_sc_skin_blood_10x_seurat_skin.rds')
pd$tp.sk@meta.data[c(6)] <- NULL
pd$tp.sk@meta.data[['sample']] <- c('Adult skin')
pd$tp.sk@meta.data$cluster <- pd$tp.sk@meta.data$RNA_snn_res.1.2
pd$tp.sk@meta.data[['RNA_snn_res.1.2']] <- NULL

# Load singleR object from healthy adult
sr <- list()

sr[['tp.sk']] <- readRDS('./GA_AN0146_sc_skin_blood_10x_singler_skin.rds')

table(sr$tp.sk[['singler']][[2]][['SingleR.single']][['labels']])

# Transfer singleR cell type annotation to Seurat objects
pd$tp.sk@meta.data[['celltype.sr']] <- sr$tp.sk[['singler']][[2]][['SingleR.single']][['labels']]


# Annotate cell types to general cell types
pd$tp.sk$subset <- plyr::mapvalues(pd$tp.sk$celltype.sr, from = levels(factor(pd$tp.sk$celltype.sr)), 
                                   to = c('CD4 T', 'CD4 Tcm', 'CD4 Tem', 'CD8 T', 'CD8 Tcm', 'CD8 Tem', 'CLP', 
                                          'DC', 'Macrophages M1','Monocytes', 'NK cells', 'Plasma cells', 'Tregs'))
pd$tp.sk$subset <- factor(pd$tp.sk$subset, levels = factor(pd$tp.sk$subset)[c(2,3,1,5,6,4,7)])
pd$tp.sk$celltype <- pd$tp.sk$subset


# Filter labels to remove non-T cells
sr$tp.sk$singler[[1]]$SingleR.single$labels[,1] <- ifelse(pd$tp.sk@meta.data$celltype %in% c('DC','Macrophages M1', 'Monocytes', 'NK cells', 'Plasma cells', 'CLP'), 'other', sr$tp.sk$singler[[1]]$SingleR.single$labels[,1])

# Transfer singleR cell type annotation to Seurat objects
pd$tp.sk@meta.data[['Compartment']] <- sr$tp.sk[['singler']][[1]][['SingleR.single']][['labels']]
pd$tp.sk@meta.data$Compartment <- factor(pd$tp.sk@meta.data$Compartment, levels = sort(as.character(unique(pd$tp.sk@meta.data$Compartment)))[c(4,5,2,3,1)])
table(pd$tp.sk@meta.data$Compartment)



#################### Visualization: adult ########################
theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = "black"),
                      axis.text.y = element_text(size = 12.8, color = "black"), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = "plain", hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = "black"), 
                      axis.ticks.y = element_line(size = 0.4, colour = "black"),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour="black"),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = "black", size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = "plain", margin = margin(0,0,10,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

# Rearrange cluster numbers to match subsets in all four samples
for(x in names(pd)) pd[[x]]$cluster_old <- pd[[x]]$celltype.sr
pd$tp.sk$celltype.sr <- plyr::mapvalues(pd$tp.sk$celltype.sr, from = levels(factor(pd$tp.sk$celltype.sr)), to = c(9,1,4,2,10,3,5,6,7,12,8,11,13))
pd$tp.sk$celltype.sr <- factor(pd$tp.sk$celltype.sr, levels = sort(as.numeric(levels(factor(pd$tp.sk$celltype.sr)))))

#-----------------------------------------------------------------------------------------------------------------#
# Vizualization in UMAP embeddings
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data to visualize UMAP embeddings
df.ts <- Map(x=names(pd), function(x) {
  data.frame(cells = colnames(pd[[x]]),
             tx = pd[[x]]@reductions$tsne@cell.embeddings[,1],
             ty = pd[[x]]@reductions$tsne@cell.embeddings[,2],
             ux = pd[[x]]@reductions$umap@cell.embeddings[,1],
             uy = pd[[x]]@reductions$umap@cell.embeddings[,2],
             sample = pd[[x]]@meta.data$sample,
             #cluster_old = pd[[x]]@meta.data$cluster_old,
             #cluster_new = factor(pd[[x]]$cluster_new, levels = sort(as.numeric(levels(factor(pd[[x]]$cluster_new))))),
             cluster = pd[[x]]@meta.data$cluster,
             compartment = pd[[x]]@meta.data$subset,
             subset = pd[[x]]@meta.data$celltype,
             row.names = 1:length(colnames(pd[[x]]))) })


# Plot tSNE/UMAP for genotypes over subsets
pal <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')
blues <- RColorBrewer::brewer.pal(n = 9, name = 'Blues')
reds <- RColorBrewer::brewer.pal(n = 9, name = 'Reds')
coln <- list(subset = list(hd.bd =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40', 'grey80'), 
                           hd.sk =  c(pal[c(2,6,8,4)], 'darkorchid3'),
                           #tp.bd =  c(pal[c(2,6,8,4)], 'darkorchid3', 'grey40'),
                           tp.bd =  c(pal[c(1,2,5,6,8,4)], 'darkorchid3', 'grey40'), 
                           tp.sk =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40')), 
             cluster = list(#hd.bd = c('cyan3', pal[c(1,2)],'royalblue4',pal[c(5)], pal[c(6,4,8)], 'darkorchid3', '#CD96CD', 'grey40', 'grey60', 'grey80'),
               hd.bd = c('cyan3', pal[c(1,2)],'royalblue4',pal[c(5)], pal[c(6,4,8)], '#CD96CD', 'darkorchid3', 'grey40', 'grey60', 'grey80'), 
               hd.sk = c(blues[c(3,4,5,7,9)],'cyan3',pal[c(5,6)],reds[c(8)], pal[c(4,7,8)],'darkorchid3'),
               tp.bd = c(pal[c(2,1)],'royalblue4',pal[c(5)],reds[c(8)],pal[c(6)],pal[c(4,7,8)],'darkorchid3','grey40','grey60','grey80'), 
               tp.sk = c(blues[c(3,5,7,9)],pal[c(5,6)],reds[c(8)],pal[c(4)],pal[c(7,8)],'#CD96CD','darkorchid3','grey40')) )

ts.sg <- Map(x=names(pd), function(x) {
  n <- length(levels(pd[[x]]@meta.data$subset))
  ggplot(df.ts[[x]] %>% arrange(subset), aes(x = ux, y = uy, color = subset, alpha = subset, size = subset)) + 
    geom_point(shape = 16, stroke = 0.5) + 
    scale_color_manual("Subset", values = c(coln$subset[[x]], 'black'), guide = guide_legend(override.aes = list(size = 4, alpha = c(rep(ifelse(x=='tp.sk',0.8,0.6),n),1)))) + 
    scale_alpha_manual("Subset", values = c(rep(ifelse(x=='tp.sk',0.8,0.3),n),1)) +
    scale_size_manual("Subset", values = c(rep(0.8,n),ifelse(x=='tp.sk',0.8,1.2))) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + #facet_wrap(~variable, ncol = 4) + 
    theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(pd[[x]]$sample) })

ts.sg.tsne <- Map(x=names(pd), function(x) {
  n <- length(levels(pd[[x]]@meta.data$subset))-1
  ggplot(df.ts[[x]] %>% arrange(subset), aes(x = tx, y = ty, color = subset, alpha = subset, size = subset)) + 
    geom_point(shape = 16, stroke = 0.5) + 
    scale_color_manual("Subset", values = c(coln$subset[[x]], 'black'), guide = guide_legend(override.aes = list(size = 4, alpha = c(rep(ifelse(x=='tp.sk',0.8,0.6),n),1)))) + 
    scale_alpha_manual("Subset", values = c(rep(ifelse(x=='tp.sk',0.8,0.3),n),1)) +
    scale_size_manual("Subset", values = c(rep(0.8,n),ifelse(x=='tp.sk',0.8,1.2))) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + #facet_wrap(~variable, ncol = 4) + 
    theme_custom + xlab('tsne 1') + ylab('tsne 2') + ggtitle(pd[[x]]$sample) })


# Define gene list
genelist <- sort(c('TCF7'))
all(genelist %in% rownames(tp.int)); genelist[!(genelist %in% rownames(tp.int))]


################################## Blood - Adult data set EX004 #######################################

####################### Load the blood dataset ############################################

blood.adult.data <- Read10X("../data/EX0004/Blood/raw_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
blood.adult <- CreateSeuratObject(counts = blood.adult.data, project = "blood", min.cells = 3, min.features = 200)
blood.adult

# Define Experiment_ID according to ngs_sample_list
experimentid <- c("EX0004")

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx("../../ngs_sample_list.xlsx", sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = "NA", fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], summarise, sum = NA)

# Sample data
metadata <- data.frame(row.names = nsl$Sample_name, 
                       reshape2::colsplit(nsl$Sample_name, "_", c("exp", "celltype", "source")))


#Reference sample name
tr <- c("blood")


##################### Define experiment specific variables ###############################

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])


# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('bd_adult_blood')
an.descs <- c('bd_adult')

# Create new folder to store figures from this analysis
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/qc'))
#dir.create(paste0('figures/', analysisid, '_', an.desc, '/pca'))
dir.figs <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs)
dir.figs.qc <- paste0('figures/', analysisid, '_', an.desc, '/qc/', analysisid, '_', an.descs)
#dir.figs.pca <- paste0('figures/', analysisid, '_', an.desc, '/pca/', analysisid, '_', an.descs)

###################### QC and selecting cells for further analysis #######################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
blood.adult[["percent.mt"]] <- PercentageFeatureSet(blood.adult, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(blood.adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(blood.adult, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(blood.adult, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


blood.adult <- subset(blood.adult, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

############################### Normalizing the data #######################################

blood.adult <- NormalizeData(blood.adult, normalization.method = "LogNormalize", scale.factor = 10000)
blood.adult <- NormalizeData(blood.adult)

#Identification of highly variable features (feature selection)
blood.adult <- FindVariableFeatures(blood.adult, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(blood.adult), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(blood.adult)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data

all.genes <- rownames(blood.adult)
blood.adult <- ScaleData(blood.adult, features = all.genes)


#Perform linear dimensional reduction
blood.adult <- RunPCA(blood.adult, features = VariableFeatures(object = blood.adult))

# Examine and visualize PCA results a few different ways
print(blood.adult[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(blood.adult, dims = 1:2, reduction = "pca")
DimPlot(blood.adult, reduction = "pca")
DimHeatmap(blood.adult, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(blood.adult, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
blood.adult <- JackStraw(blood.adult, num.replicate = 100)
blood.adult <- ScoreJackStraw(blood.adult, dims = 1:20)

JackStrawPlot(adult, dims = 1:15)
ElbowPlot(adult)

#Cluster the cells
blood.adult <- FindNeighbors(blood.adult, dims = 1:10)
blood.adult <- FindClusters(blood.adult, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(blood.adult), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
blood.adult <- RunUMAP(blood.adult, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(blood.adult, reduction = "umap", label = TRUE)

saveRDS(blood.adult, file = "./output/blood_adult_results.rds")
