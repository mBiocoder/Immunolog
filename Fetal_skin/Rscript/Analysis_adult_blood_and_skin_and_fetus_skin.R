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

################################## Making rds files for downstream analysis using Seurat pipeline ################

############## fetal skin #################
# Load UMI counts from 10x data
skin.fetus.data <- Read10X_h5("../data/Paper/GSE156972_raw_gene_bc_matrices_h5.h5")

# Initialize the Seurat object with the raw (non-normalized data).
skin.fetus <- CreateSeuratObject(counts = skin.fetus.data, min.cells = 0, min.features = 0, project = "skin-fetus")

# Calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`
mito.features <- grep(pattern = "^MT-", x = rownames(skin.fetus), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = skin.fetus, slot = 'counts')[mito.features, ]) / 
  Matrix::colSums(GetAssayData(object = skin.fetus, slot = 'counts'))
skin.fetus[['percent.mito']] <- percent.mito
#VlnPlot(skin.fetus, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Visualize feature-feature relationships
FeatureScatter(skin.fetus, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black", group.by = "orig.ident")
FeatureScatter(skin.fetus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black", group.by = "orig.ident")

# Filter out cells that have unique feature counts over 2000 or less than 200
skin.fetus <- subset(skin.fetus, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 0.2)
#set.seed(1); bd <- subset(bd, cells = sample(colnames(bd), size = 1000)) #length(colnames(bd)); head(colnames(bd),3)
vp.skin.fetus <- VlnPlot(skin.fetus, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Employ a global-scaling normalization method that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10000 by default), and log-transforms the result
skin.fetus <- NormalizeData(skin.fetus, normalization.method = "LogNormalize", scale.factor = 1e4)

# Detection of variable features across the single cells
# Calculate the average expression and dispersion for each feature and z-score for dispersion within each bin
#bd <- FindVariableFeatures(bd, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), 
#                           dispersion.cutoff = c(0.5, Inf))
skin.fetus <- FindVariableFeatures(skin.fetus, selection.method = "vst", nfeatures = 8000, 
                           mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(skin.fetus))
vf.skin.fetus <- VariableFeaturePlot(skin.fetus, cols = c("black", "red"), pt.size = 1, log = NULL, assay = NULL)

# Scaling the data and removing unwanted sources of variation (stored in the scale.data slot)
skin.fetus <- ScaleData(skin.fetus, features = rownames(skin.fetus), vars.to.regress = c("nCount_RNA", "percent.mito"))

# Perform linear dimensional reduction
skin.fetus <- RunPCA(skin.fetus, features = VariableFeatures(skin.fetus), verbose = F)

# Examine and visualize PCA results a few different ways (VizDimReduction, DimPlot, and DimHeatmap)
print(skin.fetus[['pca']], dims = 1:5, nfeatures = 5, projected = F)
VizDimLoadings(skin.fetus, dims = c(1,2))
DimPlot(skin.fetus, reduction = "pca")

# Scores each feature in the dataset (including features not included in the PCA) based on their correlation with the PCs
skin.fetus <- ProjectDim(skin.fetus, do.center = T, verbose = F)
print(skin.fetus[['pca']], dims = 1:5, nfeatures = 5, projected = T)
VizDimLoadings(skin.fetus, dims = c(1,2), projected = T)
DimPlot(skin.fetus, projected = T, reduction = "pca")

# Explore the primary sources of heterogeneity to decide which PCs to include for further downstream analyses
DimHeatmap(skin.fetus, dims = 1, cells = NULL, balanced = T)
#DimHeatmap(bd, dims = 1:20, cells = 500, balanced = T)
hm.pc.skin.fetus <- DimHeatmap(skin.fetus, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle("Heatmap of PC markers from 500 cells (blood)")

# Determine statistically significant principal components
#bd <- JackStraw(bd, dims = 20, prop.freq = 0.01, num.replicate = 100, maxit = 1000)
#bd <- ScoreJackStraw(bd, dims = 1:20)
#JackStrawPlot(bd, dims = 1:20)
ep.bd <- ElbowPlot(skin.fetus, ndims = 25)

# Cluster cells (the ‘granularity’ (number of clusters) of the downstream clustering is set to 0.6-1.2 for 3K cells)
skin.fetus <- FindNeighbors(skin.fetus, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
skin.fetus <- FindClusters(skin.fetus, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Run Non-linear dimensional reduction (tSNE and UMAP)
skin.fetus <- RunTSNE(skin.fetus, reduction = "pca", dims = 1:20, seed.use = 10, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30)
skin.fetus <- RunUMAP(skin.fetus, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
DimPlot(skin.fetus, reduction = "tsne", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 8, pt.size = 0.8, cols = NULL)
DimPlot(skin.fetus, reduction = "umap", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 6, pt.size = 0.8, cols = NULL) 
DimPlot(skin.fetus, reduction = "pca") 

skin.fetus[["clust_1"]] <- Idents(object = skin.fetus)
ts.sp.skin.fetus <- DimPlot(skin.fetus, reduction = "umap", dims = c(1,2), split.by = "clust_1", 
                    label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_y_continuous(breaks = c(-25,0,25)) 
vp.skin.fetusc <- VlnPlot(skin.fetus, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# Finding differentially expressed features (cluster biomarkers)
skin.fetus.markers <- FindAllMarkers(skin.fetus, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                             test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
skin.fetus.top <- skin.fetus.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gn <- c("IL7R", "CD4", "CD8A", "CD8B", "GNLY", "MS4A1", "FCGR3A", "NKG7", "GZMB")

# Distinguish naive and memory
FeaturePlot(object = skin.fetus, features = c("S100A4", "CCR7"), cols = c("lightgrey", "blue"), coord.fixed = T)

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
VlnPlot(skin.fetus, features = skin.fetus.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
VlnPlot(skin.fetus, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(skin.fetus, features = gn, dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "umap", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.skin.fetus <- DoHeatmap(skin.fetus, features = skin.fetus.top$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers")

#Save rds file
saveRDS(skin.fetus, file = "./output/Analysis_adult_blood_skin_and_fetus_skin/skin.fetus.rds")

############## adult skin #################
# Load UMI counts from 10x data
skin.adult.data <- Read10X("../data/EX0004/raw_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
skin.adult <- CreateSeuratObject(counts = skin.adult.data, min.cells = 0, min.features = 0, project = "skin-adult")

# Calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`
mito.features <- grep(pattern = "^MT-", x = rownames(skin.adult), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = skin.adult, slot = 'counts')[mito.features, ]) / 
  Matrix::colSums(GetAssayData(object = skin.adult, slot = 'counts'))
skin.adult[['percent.mito']] <- percent.mito

# Visualize feature-feature relationships
FeatureScatter(skin.adult, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black", group.by = "orig.ident")
FeatureScatter(skin.adult, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black", group.by = "orig.ident")

# Filter out cells that have unique feature counts over 2000 or less than 200
skin.adult <- subset(skin.adult, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 0.2)
#set.seed(1); bd <- subset(bd, cells = sample(colnames(bd), size = 1000)) #length(colnames(bd)); head(colnames(bd),3)
vp.skin.adult <- VlnPlot(skin.adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Employ a global-scaling normalization method that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10000 by default), and log-transforms the result
skin.adult <- NormalizeData(skin.adult, normalization.method = "LogNormalize", scale.factor = 1e4)

# Detection of variable features across the single cells
# Calculate the average expression and dispersion for each feature and z-score for dispersion within each bin
#bd <- FindVariableFeatures(bd, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), 
#                           dispersion.cutoff = c(0.5, Inf))
skin.adult <- FindVariableFeatures(skin.adult, selection.method = "vst", nfeatures = 8000, 
                                   mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(skin.adult))
vf.skin.adult <- VariableFeaturePlot(skin.adult, cols = c("black", "red"), pt.size = 1, log = NULL, assay = NULL)

# Scaling the data and removing unwanted sources of variation (stored in the scale.data slot)
skin.adult <- ScaleData(skin.adult, features = rownames(skin.adult), vars.to.regress = c("nCount_RNA", "percent.mito"))

# Perform linear dimensional reduction
skin.adult <- RunPCA(skin.adult, features = VariableFeatures(skin.adult), verbose = F)

# Examine and visualize PCA results a few different ways (VizDimReduction, DimPlot, and DimHeatmap)
print(skin.adult[['pca']], dims = 1:5, nfeatures = 5, projected = F)
VizDimLoadings(skin.adult, dims = c(1,2))
DimPlot(skin.adult, reduction = "pca")

# Scores each feature in the dataset (including features not included in the PCA) based on their correlation with the PCs
skin.adult <- ProjectDim(skin.adult, do.center = T, verbose = F)
print(skin.adult[['pca']], dims = 1:5, nfeatures = 5, projected = T)
VizDimLoadings(skin.adult, dims = c(1,2), projected = T)
DimPlot(skin.adult, projected = T, reduction = "pca")

# Explore the primary sources of heterogeneity to decide which PCs to include for further downstream analyses
DimHeatmap(skin.adult, dims = 1, cells = NULL, balanced = T)
#DimHeatmap(bd, dims = 1:20, cells = 500, balanced = T)
hm.pc.skin.adult <- DimHeatmap(skin.adult, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle("Heatmap of PC markers from 500 cells (blood)")

# Determine statistically significant principal components
#bd <- JackStraw(bd, dims = 20, prop.freq = 0.01, num.replicate = 100, maxit = 1000)
#bd <- ScoreJackStraw(bd, dims = 1:20)
#JackStrawPlot(bd, dims = 1:20)
ep.bd <- ElbowPlot(skin.adult, ndims = 25)

# Cluster cells (the ‘granularity’ (number of clusters) of the downstream clustering is set to 0.6-1.2 for 3K cells)
skin.adult <- FindNeighbors(skin.adult, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
skin.adult <- FindClusters(skin.adult, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Run Non-linear dimensional reduction (tSNE and UMAP)
skin.adult <- RunTSNE(skin.adult, reduction = "pca", dims = 1:20, seed.use = 10, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30)
skin.adult <- RunUMAP(skin.adult, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
DimPlot(skin.adult, reduction = "tsne", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 8, pt.size = 0.8, cols = NULL)
DimPlot(skin.adult, reduction = "umap", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 6, pt.size = 0.8, cols = NULL) 
DimPlot(skin.adult, reduction = "pca") 

skin.adult[["clust_1"]] <- Idents(object = skin.adult)
ts.sp.skin.adult <- DimPlot(skin.adult, reduction = "umap", dims = c(1,2), split.by = "clust_1", 
                            label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_y_continuous(breaks = c(-25,0,25)) 
vp.skin.adultc <- VlnPlot(skin.adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# Finding differentially expressed features (cluster biomarkers)
skin.adult.markers <- FindAllMarkers(skin.adult, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                                     test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
skin.adult.top <- skin.adult.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gn <- c("IL7R", "CD4", "CD8A", "CD8B", "GNLY", "MS4A1", "FCGR3A", "NKG7", "GZMB")

# Distinguish naive and memory
FeaturePlot(object = skin.adult, features = c("S100A4", "CCR7"), cols = c("lightgrey", "blue"), coord.fixed = T)

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
VlnPlot(skin.adult, features = skin.adult.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
VlnPlot(skin.adult, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(skin.adult, features = gn, dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "umap", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.skin.adult <- DoHeatmap(skin.adult, features = skin.adult.top$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                              size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers")

#Save rds file
saveRDS(skin.adult, file = "./output/Analysis_adult_blood_skin_and_fetus_skin/skin.adult.rds")

############## adult blood #################
# Load UMI counts from 10x data
blood.adult.data <- Read10X("../data/EX0004/Blood/raw_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
blood.adult <- CreateSeuratObject(counts = blood.adult.data, min.cells = 0, min.features = 0, project = "blood-adult")

# Calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`
mito.features <- grep(pattern = "^MT-", x = rownames(blood.adult), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = blood.adult, slot = 'counts')[mito.features, ]) / 
  Matrix::colSums(GetAssayData(object = blood.adult, slot = 'counts'))
blood.adult[['percent.mito']] <- percent.mito

# Visualize feature-feature relationships
FeatureScatter(blood.adult, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black", group.by = "orig.ident")
FeatureScatter(blood.adult, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black", group.by = "orig.ident")

# Filter out cells that have unique feature counts over 2000 or less than 200
blood.adult <- subset(blood.adult, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 0.2)
#set.seed(1); bd <- subset(bd, cells = sample(colnames(bd), size = 1000)) #length(colnames(bd)); head(colnames(bd),3)
vp.blood.adult <- VlnPlot(blood.adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Employ a global-scaling normalization method that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10000 by default), and log-transforms the result
blood.adult <- NormalizeData(blood.adult, normalization.method = "LogNormalize", scale.factor = 1e4)

# Detection of variable features across the single cells
# Calculate the average expression and dispersion for each feature and z-score for dispersion within each bin
#bd <- FindVariableFeatures(bd, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), 
#                           dispersion.cutoff = c(0.5, Inf))
blood.adult <- FindVariableFeatures(blood.adult, selection.method = "vst", nfeatures = 8000, 
                                   mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(blood.adult))
vf.blood.adult <- VariableFeaturePlot(blood.adult, cols = c("black", "red"), pt.size = 1, log = NULL, assay = NULL)

# Scaling the data and removing unwanted sources of variation (stored in the scale.data slot)
blood.adult <- ScaleData(blood.adult, features = rownames(blood.adult), vars.to.regress = c("nCount_RNA", "percent.mito"))

# Perform linear dimensional reduction
blood.adult <- RunPCA(blood.adult, features = VariableFeatures(blood.adult), verbose = F)

# Examine and visualize PCA results a few different ways (VizDimReduction, DimPlot, and DimHeatmap)
print(blood.adult[['pca']], dims = 1:5, nfeatures = 5, projected = F)
VizDimLoadings(blood.adult, dims = c(1,2))
DimPlot(blood.adult, reduction = "pca")

# Scores each feature in the dataset (including features not included in the PCA) based on their correlation with the PCs
blood.adult <- ProjectDim(blood.adult, do.center = T, verbose = F)
print(blood.adult[['pca']], dims = 1:5, nfeatures = 5, projected = T)
VizDimLoadings(blood.adult, dims = c(1,2), projected = T)
DimPlot(blood.adult, projected = T, reduction = "pca")

# Explore the primary sources of heterogeneity to decide which PCs to include for further downstream analyses
DimHeatmap(blood.adult, dims = 1, cells = NULL, balanced = T)
#DimHeatmap(bd, dims = 1:20, cells = 500, balanced = T)
hm.pc.blood.adult <- DimHeatmap(blood.adult, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle("Heatmap of PC markers from 500 cells (blood)")

# Determine statistically significant principal components
#bd <- JackStraw(bd, dims = 20, prop.freq = 0.01, num.replicate = 100, maxit = 1000)
#bd <- ScoreJackStraw(bd, dims = 1:20)
#JackStrawPlot(bd, dims = 1:20)
ep.bd <- ElbowPlot(blood.adult, ndims = 25)

# Cluster cells (the ‘granularity’ (number of clusters) of the downstream clustering is set to 0.6-1.2 for 3K cells)
blood.adult <- FindNeighbors(blood.adult, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
blood.adult <- FindClusters(blood.adult, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Run Non-linear dimensional reduction (tSNE and UMAP)
blood.adult <- RunTSNE(blood.adult, reduction = "pca", dims = 1:20, seed.use = 10, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30)
blood.adult <- RunUMAP(blood.adult, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
DimPlot(blood.adult, reduction = "tsne", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 8, pt.size = 0.8, cols = NULL)
DimPlot(blood.adult, reduction = "umap", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 6, pt.size = 0.8, cols = NULL) 
DimPlot(blood.adult, reduction = "pca") 

blood.adult[["clust_1"]] <- Idents(object = blood.adult)
ts.sp.blood.adult <- DimPlot(blood.adult, reduction = "umap", dims = c(1,2), split.by = "clust_1", 
                            label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_y_continuous(breaks = c(-25,0,25)) 
vp.blood.adultc <- VlnPlot(blood.adult, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# Finding differentially expressed features (cluster biomarkers)
blood.adult.markers <- FindAllMarkers(blood.adult, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                                     test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
blood.adult.top <- blood.adult.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gn <- c("IL7R", "CD4", "CD8A", "CD8B", "GNLY", "MS4A1", "FCGR3A", "NKG7", "GZMB")

# Distinguish naive and memory
FeaturePlot(object = blood.adult, features = c("S100A4", "CCR7"), cols = c("lightgrey", "blue"), coord.fixed = T)

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
VlnPlot(blood.adult, features = blood.adult.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
VlnPlot(blood.adult, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(blood.adult, features = gn, dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "umap", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.blood.adult <- DoHeatmap(blood.adult, features = blood.adult.top$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                              size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers")

#Save rds file
saveRDS(blood.adult, file = "./output/Analysis_adult_blood_skin_and_fetus_skin/blood.adult.rds")

################################### Integrated analysis blood adult, skin adult and fetus ########################

#load rds file to object
skin_adult <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/skin.adult.rds") #skin adult
blood_adult <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/blood.adult.rds") # blood adult
skin_fetus <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/skin.fetus.rds") #fetal skin


ifnb.list <- list(blood_adult, skin_adult, skin_fetus)

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
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
immune.combined <- FindClusters(immune.combined, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

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

# Marker genes -same markers as Gustavo used

gn <- c('CD3E','CD4','CD8A', 'CD8B', 'CCR7', 'FCER1A', 'FCGR3A', 'NKG7', 'GZMB', 'GZMA', 
        'GNLY', 'MS4A1', 'CD14', 'LYZ', 'CD34', 'IL7R', 'PDCD1', 'B3GAT1', 'PECAM1', 'ICOS', 
        'IL2RA', 'ITGAE', 'CD69', 'CXCR6', 'ITGA1', 'PTPRC', 'CX3CR1', 'SELL', 'CD19', 'CXCR5', 
        'CCR6', 'FAS', 'KLRB1', 'CD27', 'CD28', 'TRDC', 'TRGC1', 'TRGC2', 'ITGB1', 'SELPLG', 
        'CCR2', 'CST3', 'MS4A7')#, '', '', '', '', '')
names(gn) <- c('T cell', 'CD4 T cell', 'CD8 T cell', 'CD8 T cell', 'Tn', 'DC', 'Monocyte', 'NK cell', 'NK cell', 'NK cell', 
               'NK cell', 'B cell', 'Monocyte', 'Monocyte', 'HSC', 'T cell', 'Tm', 'Tem', 'Tn', 'Treg', 
               'Treg', 'Trm', 'Trm', 'Trm', 'Trm', 'Monocyte', 'trafficking', 'Tn', 'B cell', 'B cell', 
               'B cell', 'Tm', 'RORgT cell', 'Tn and Tcm', 'Tn and Tcm', 'gdT cell', 'gdT cell', 'gdT cell', 'trafficking', 'trafficking', 
               'Monocyte', 'DC', 'Monocyte')

gn

#Save rds file
saveRDS(immune.combined,"./output/Analysis_adult_blood_skin_and_fetus_skin/immune_combined_all_three.rds")
#Read in rds file
immune.combined <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/immune_combined_all_three.rds")

#Plot PCA plot 
DimPlot(immune.combined, split.by = "orig.ident", reduction = "pca")
DimPlot(immune.combined, reduction = "pca")

#Plot UMAP for 3 origins separately
DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")

#Plot dot plot for residency markers
markers.to.plot <- c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1")
DotPlot(immune.combined, features = markers.to.plot, cols = c("red", "blue", "grey"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()

#Plot dot plot for cell type annotation; Genes from Gustavo's dotplot from him publication
markers.to.plot <- c("CD3E", "TRAC", "TRBC1", "TRBC2", "CD4", "IL7R", "CD8A", "CD8B", "CTLA4", "FOXP3", "ICOS", "TRDC", "TRGC1", "TRGC2", "KLRB1", "GNLY", "GZMB", "MS4A1")
DotPlot(immune.combined, features = markers.to.plot, cols = c("red", "blue", "grey"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()


# Plot results in order to assign the clusters the specific cell types
#B cells, T cells, Dendritic cells, NK cells, gdT, NKT and monocytes we will visualize first

#T cells alone
FeaturePlot(immune.combined, features = c("CD3E", "IL7R"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD3E", "IL7R"))

#Tn
FeaturePlot(immune.combined, features = c("TCF7", "CCR7", "PECAM1", "SELL"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("TCF7", "CCR7", "PECAM1", "SELL"))

#CD4 T
FeaturePlot(immune.combined, features = c("CD4"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD4"))

#CD8 T
FeaturePlot(immune.combined, features = c("CD8A", "CD8B"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD8A", "CD8B"))

#DC
FeaturePlot(immune.combined, features = c("FCER1A", "CST3"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("FCER1A", "CST3"))

#Monocyte
FeaturePlot(immune.combined, features = c("CD14", "LYZ", "PTPRC", "CCR2", "MS4A7"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD14", "LYZ", "PTPRC", "CCR2", "MS4A7"))

#NK
FeaturePlot(immune.combined, features = c("NKG7", "GZMB", "GZMA", "GNLY"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("NKG7", "GZMB", "GZMA", "GNLY"))

#NKT
FeaturePlot(immune.combined, features = c("KLRB1", "GZMK"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("KLRB1", "GZMK"))

#B cells
FeaturePlot(immune.combined, features = c("MS4A1", "CD19", "CXCR5"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("MS4A1", "CD19", "CXCR5"))

#HSC
FeaturePlot(immune.combined, features = c("CD34"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD34"))

#Tm
FeaturePlot(immune.combined, features = c("PDCD1", "FAS"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("PDCD1", "FAS"))

#Tem
FeaturePlot(immune.combined, features = c("B3GAT1"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("B3GAT1"))

#Treg
FeaturePlot(immune.combined, features = c("ICOS", "IL2RA"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("ICOS", "IL2RA"))

#Trm
FeaturePlot(immune.combined, features = c("ITGAE", "CD69", "CXCR6", "ITGA1"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("ITGAE", "CD69", "CXCR6", "ITGA1"))

#trafficing
FeaturePlot(immune.combined, features = c("CX3CR1", "ITGB1", "SELPLG"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CX3CR1", "ITGB1", "SELPLG"))

#RORgT
FeaturePlot(immune.combined, features = c("KLRB1"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("KLRB1"))

#abT
FeaturePlot(immune.combined, features = c("TRAC", "TRBC1", "TRBC2"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("TRAC", "TRBC1", "TRBC2"))

#Tn and Tcm
FeaturePlot(immune.combined, features = c("CD27", "CD28"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("CD27", "CD28"))

#gdT
FeaturePlot(immune.combined, features = c("TRDC", "TRGC1", "TRGC2"), label = TRUE, min.cutoff = 'q10')
VlnPlot(immune.combined, features = c("TRDC", "TRGC1", "TRGC2"))


#Assign cell type identities
immune.combined <- RenameIdents(immune.combined, `0` = "CD4 T", `1` = "CD4 T", `2` = "CD4 T",
                                `3` = "CD4 T", `4` = "CD8 T", `5` = "NKT", `6` = "NKT", `7` = "CD8 T", `8` = "CD4 T", `9` = "CD4 T",
                                `10` = "NK", `11` = "NKT", `12` = "Treg", `13` = "B", `14` = "other", `15` = "gdT",`16` = "CD4 T", 
                                `17` = "other", `18` = "other", `19` = "CD4 T", `20` = "other", `21` = "other", `22` = "other", 
                                `23` = "other", `24` = "other")

#Plot separated UMAP
DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")

#Plot UMAP
DimPlot(immune.combined, label = TRUE)

# How many cells are in each cluster
table(Idents(immune.combined))

# How many cells are in each condition?
table(immune.combined$orig.ident)

# What proportion of cells are in each cluster?
prop.table(table(Idents(immune.combined)))

# How does cluster membership vary by conditon?
table1 <- prop.table(table(Idents(immune.combined), immune.combined$orig.ident), margin = 2)
table1 <- melt(table1, id.vars = 0, measure.vars = c("adult", "blood", "fetal"),
               variable.name = "origin", value.name = "value")

#change table1 colnames
colnames(table1) <- c("subsets", "origin", "Percentage_of_total_cells")
table1

# Stacked barplot + percent
ggplot(table1, aes(fill=subsets, y=Percentage_of_total_cells, x=origin)) + 
  geom_bar(position="fill", stat="identity")  + scale_y_continuous(labels = scales::percent_format())
  
 

#Retain only Cd4 T, Cd8 T and Tregs as ab+T cells
#Remove all non abT+-cells
ab_Tcell_subset <- subset(immune.combined,
                          idents = c("B", "other", "gdT", "NK", "NKT"), invert = TRUE)


# Re-visualize the clusters
DimPlot(object = ab_Tcell_subset, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

# Re-visualize the clusters separated by orig.ident
DimPlot(object = ab_Tcell_subset, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE, split.by = "orig.ident")

#Save rds file
saveRDS(ab_Tcell_subset,"./output/Analysis_adult_blood_skin_and_fetus_skin/ab_Tcell_subset.rds") 
#Read in rds file
ab_Tcell_subset <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/ab_Tcell_subset.rds")

# How many cells are in each cluster
table(Idents(ab_Tcell_subset))

# How many cells are in each condition?
table(ab_Tcell_subset$orig.ident)

# What proportion of cells are in each cluster?
prop.table(table(Idents(ab_Tcell_subset)))

# How does cluster membership vary by condition?
prop.table(table(Idents(ab_Tcell_subset), ab_Tcell_subset$orig.ident), margin = 2)


#Plot IL1A and IL1B in UMAP
DefaultAssay(ab_Tcell_subset) <- "RNA"
FeaturePlot(ab_Tcell_subset, features = c("IL1A"), split.by = "orig.ident", min.cutoff = 'q10')

################ Identify differential expressed genes across conditions ###########

library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

#Subset ab_Tcell_subset to contain only blood adult data
#blood_adult_data <- subset(x = ab_Tcell_subset, subset = orig.ident == "blood-adult")

# Subset ab_Tcell_subset to contain only skin adult and skin fetus data
skin_fetus_and_adult_data <-  subset(x = ab_Tcell_subset, subset = orig.ident == "blood-adult", invert = TRUE)

avg.t.cells <- as.data.frame(log1p(AverageExpression(skin_fetus_and_adult_data, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(skin_fetus_and_adult_data)

#DEG skin-fetus vs. skin-adult, only search for positive markers and pick top 50(or 10)
Idents(skin_fetus_and_adult_data) <- skin_fetus_and_adult_data$orig.ident
deg.gene.pos <- FindMarkers(skin_fetus_and_adult_data, ident.1="skin-fetus" , ident.2 ="skin-adult", only.pos = TRUE) %>%
top_n(n = 50, wt = avg_log2FC) -> top50  #top_n(n = 10, wt = avg_log2FC) -> top10
deg.gene.pos
#head(deg.gene.pos, n = 50) #head(deg.gene.pos, n = 10)


#DEG skin-fetus vs. skin-adult, only search for negative markers and pick top 50 (10)
deg.gene.neg <- FindMarkers(skin_fetus_and_adult_data, ident.1="skin-fetus" , ident.2 ="skin-adult", only.pos = FALSE ) %>%
  top_n(n = 50, wt = -avg_log2FC) -> top50  #top_n(n = 10, wt = avg_log2FC) -> top10
deg.gene.neg
#Reorder deg.gene.neg to get top 10 negative markers
#newdata <- deg.gene.neg[order(deg.gene.neg$avg_log2FC),]
#head(newdata, n = 10)

#write to tsv file
library("xlsx")

#upregulated DEGs
write.xlsx(deg.gene.pos, file = "DEG_for_naivity_signatures.xlsx", sheetName = " 50 upregulated DEGs", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

#downregulated DEGs
write.xlsx(deg.gene.neg, file = "DEG_for_naivity_signatures.xlsx", sheetName = " 50 downregulated DEGs", 
col.names = TRUE, row.names = TRUE, append = TRUE)


# Re-run the standard workflow for visualization and clustering subsetted Seurat object again
skin_fetus_and_adult_data <- FindVariableFeatures(skin_fetus_and_adult_data, selection.method = "vst")
DefaultAssay(skin_fetus_and_adult_data) <- "integrated"

skin_fetus_and_adult_data <- ScaleData(skin_fetus_and_adult_data, verbose = FALSE)
skin_fetus_and_adult_data <- RunPCA(skin_fetus_and_adult_data, npcs = 30, verbose = FALSE)
skin_fetus_and_adult_data <- RunUMAP(skin_fetus_and_adult_data, reduction = "pca", dims = 1:30)
skin_fetus_and_adult_data <- FindNeighbors(skin_fetus_and_adult_data, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
skin_fetus_and_adult_data <- FindClusters(skin_fetus_and_adult_data, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Run Non-linear dimensional reduction (UMAP)
skin.fetus.and.adult.data <- RunUMAP(skin_fetus_and_adult_data, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
DimPlot(skin.fetus.and.adult.data, reduction = "umap", dims = c(1,2), split.by = "orig.ident", label = T, label.size = 6, pt.size = 0.8, cols = NULL) 
DimPlot(skin.fetus.and.adult.data, reduction = "pca", split.by = "orig.ident", label = T, label.size = 6, pt.size = 0.8) 


#Save rds file
saveRDS(skin_fetus_and_adult_data,"./output/Analysis_adult_blood_skin_and_fetus_skin/skin_fetus_and_adult_data.rds") 
#Read in rds file
skin_fetus_and_adult_data <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/skin_fetus_and_adult_data.rds")


#Plot IL1A and IL1B in UMAP
DefaultAssay(skin_fetus_and_adult_data) <- "RNA"
FeaturePlot(skin_fetus_and_adult_data, features = c("IL1A", "IL1B"), split.by = "orig.ident", min.cutoff = 'q10')  & theme(legend.position = c(1,0.6))

# Subset on all cells expressing high abTcell marker genes AND IL1B
testing2 <- subset(x = skin_fetus_and_adult_data, subset = IL1B > 2)
FeaturePlot(testing2, features = c( "CD3E", "TRAC", "TRBC1", "TRBC2"), split.by = "orig.ident") 

#Featureplot to check Tcf7 expression
FeaturePlot(skin_fetus_and_adult_data, features = c("TCF7"), split.by = "orig.ident", min.cutoff = 1)  & theme(legend.position = c(0.7,0.6))


#Featureplot to check for cytokine expression
FeaturePlot(skin_fetus_and_adult_data, features = c("IL2", "IL4", "IL13", "IL17", "TFNB", "IFNG"), split.by = "orig.ident")

#Featureplot for other (from UCell) T cell marker genes
FeaturePlot(skin_fetus_and_adult_data, features = c( "CD3D" , "CD3E",  "CD3G" ,  "CD4" ,  "CD2" ,  "CD7",  "TRAC" ,"TRBC1" ,  "LAT" ), split.by = "orig.ident")  & theme(legend.position = c(1,0.6))

#Dotplot and featureplot to investigate IL1B expressing clusters 13, 16 and 17 closer
markers.to.plot <- c("CD3E", "TRAC", "TRBC1", "TRBC2")
DotPlot(skin_fetus_and_adult_data, features = markers.to.plot, dot.scale = 8, split.by = "orig.ident") + RotatedAxis() 
FeaturePlot(skin_fetus_and_adult_data, features = c( "CD3E", "TRAC", "TRBC1", "TRBC2" ), split.by = "orig.ident")


# Visualize co-expression of two features simultaneously
FeaturePlot(skin_fetus_and_adult_data, features = c("IL1B", "CD3E"), blend = TRUE, cols = c("red", "blue"), split.by = "orig.ident")
FeaturePlot(skin_fetus_and_adult_data, features = c("IL1B", "TRAC"), blend = TRUE, cols = c("red", "blue"), split.by = "orig.ident")
FeaturePlot(skin_fetus_and_adult_data, features = c("IL1B", "TRBC1"), blend = TRUE, cols = c("red", "blue"), split.by = "orig.ident")
FeaturePlot(skin_fetus_and_adult_data, features = c("IL1B", "TRBC2"), blend = TRUE, cols = c("red", "blue"), split.by = "orig.ident")

FeaturePlot(skin_fetus_and_adult_data, features = c("TCF7", "CCR7"),split.by = "orig.ident")
VlnPlot(object = skin_fetus_and_adult_data, features = 'TCF7', split.by = 'orig.ident')

#Subset to remove skin-fetus clusters 15, 17 and 18
subsetted <- subset(x = skin_fetus_and_adult_data, idents = c("15", "17", "18"), invert = TRUE)

#Save rds file
saveRDS(subsetted,"./output/Analysis_adult_blood_skin_and_fetus_skin/abTcell_subsetted_skin_fetus_and_adult_data.rds") 
#Read in rds file
subsetted <- readRDS("./output/Analysis_adult_blood_skin_and_fetus_skin/abTcell_subsetted_skin_fetus_and_adult_data.rds")


# Re-visualize the clusters separated by orig.ident
DimPlot(object = subsetted, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE, split.by = "orig.ident")


#Subset to only clusters 13 and 16 where IL1B expression was high in skin-fetus
IL1B.clusters <- subset(x = subsetted, idents = c("13", "16"))

#Plot dot plot for cell type annotation; Genes from Gustavo's dotplot from him publication
markers.to.plot <- c("CD3E", "TRAC", "TRBC1", "TRBC2", "CD4", "IL7R", "CD8A", "CD8B", "CTLA4", "FOXP3", "ICOS", "TRDC", "TRGC1", "TRGC2", "KLRB1", "GNLY", "GZMB", "MS4A1", "FCER1A", "CST3", "CD14", "LYZ",
                     "PTPRC", "CCR2", "NKG7")
DotPlot(IL1B.clusters, features = markers.to.plot, dot.scale = 8, split.by = "orig.ident") + RotatedAxis()

DotPlot(immune.combined, features = markers.to.plot, cols = c("grey", "light blue", "blue"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()


markers.to.plot <- c()
DotPlot(IL1B.clusters, features = markers.to.plot, dot.scale = 8, split.by = "orig.ident") + RotatedAxis()


#Find DEG for clusters 13 and 16
DEG.IL1B.clusters <- FindAllMarkers(IL1B.clusters, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                                     test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
DEG.IL1B.clusters.top <- DEG.IL1B.clusters %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
View(DEG.IL1B.clusters.top)

#Do enrichment analysis for these top DEGs


#Plot TCf7 expression
FeaturePlot(subsetted, features = c("TCF7"),split.by = "orig.ident")
VlnPlot(object = subsetted, features = 'TCF7', split.by = 'orig.ident')


library(dplyr)
library(tidyverse)
p <- VlnPlot(object = subsetted, features =c("TCF7"), split.by = "orig.ident")
p$data %>% group_by(split) %>% summarize(counts = sum(TCF7, na.rm = TRUE))


#Find percentage of TCF7 expressing cells comparing all cells
#PercentageFeatureSet(object = subsetted, features = "TCF7")
sum(GetAssayData(object = subsetted, slot = "data")["TCF7",]>0) #number of cells expressing IL1B is 821
sum(GetAssayData(object = subsetted, slot = "data")["TCF7",]>0)/nrow(subsetted@meta.data) #percentage is thus 0.1148252


PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(subsetted, "TCF7", group.by = "orig.ident")

#Check IL1B and abT cells features in newly subsetted data
#FeaturePlot(subsetted, features = c("IL1B"), split.by = "orig.ident", min.cutoff = 'q10')
#FeaturePlot(subsetted, features = c( "CD3E", "TRAC", "TRBC1", "TRBC2" ), split.by = "orig.ident")
#markers.to.plot <- c("CD3E", "TRAC", "TRBC1", "TRBC2")
#DotPlot(subsetted, features = markers.to.plot, dot.scale = 8, split.by = "orig.ident") + RotatedAxis() 
#FeaturePlot(subsetted, features = c("TCF7", "CCR7"),split.by = "orig.ident") 


##################### Checking for doublets ##############################
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
suppressMessages(require(DoubletFinder))

#split subsetted object which is integrated data and can thus not be used into two Seurat obects
split.list <- SplitObject(subsetted, split.by = "orig.ident")


skin.adult.subsetted = FindVariableFeatures(split.list$`skin-adult`, verbose = F)
skin.adult.subsetted = ScaleData(skin.adult.subsetted, vars.to.regress = c("nFeature_RNA", "percent.mito"),
                      verbose = F)
skin.adult.subsetted = RunPCA(skin.adult.subsetted, verbose = F, npcs = 20)
skin.adult.subsetted = RunUMAP(skin.adult.subsetted, dims = 1:10, verbose = F)

# Run parameter optimization with paramSweep -> should NOT be run on integrated data
sweep.res <- paramSweep_v3(skin.adult.subsetted)sweep.stats <-
summarizeSweep(sweep.res, GT = FALSE) bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cellscells.
nExp <- round(ncol(subsetted) * 0.04)  # expect 4% doublets
subsetted <- doubletFinder_v3(subsetted, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

##########################################################################

#Check residency markers
#VlnPlot(testing, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
#FeaturePlot(testing, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))




################################################## Enrichment analysis ###################################################
library(enrichR)
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
if (websiteLive) {
  enriched <- enrichr(c("FCER1G","NBEAL1", "MTRNR2L12","SPP1", "AIF1", "ACTB","IFITM3",
                        "HBG2","RACK1", "ATP5F1E","NOP53", "CXCR4", "ATP5MG", "SELENOK", 
                        "TUBA4A", "ZFP36L2", "BTG1", "TSC22D3", "HBA2", "HBA1"), dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2015"]]

#Plot enrichment barplot
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")


###### Run enrichment analysis on all DEGs ########
deg.gene.all <- FindMarkers(skin_fetus_and_adult_data, ident.1="skin-fetus" , ident.2 ="skin-adult", only.pos = FALSE )

#list <- rownames(deg.gene.all)
if (websiteLive) {
  enriched <- enrichr(rownames(deg.gene.all), dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2015"]]

#Plot enrichment barplot
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")


#DOSE
library(DOSE)
#make a single vector containing up -and downregulated top 50 DEG genes
dge_list <- append(rownames(deg.gene.pos), rownames(deg.gene.neg))
dge_list_full <- rbind(deg.gene.pos, deg.gene.neg)


#convert to entrezid from genesymbol
dge_list_full$entrezid <- mapIds(x = org.Hs.eg.db,
                                         keys = rownames(dge_list_full),
                                         column = "ENTREZID",
                                         keytype = "SYMBOL",
                                         multiVals = "first")

dge_list_full <- na.omit(dge_list_full)
dge_list_full$avg_log2FC = sort(dge_list_full$avg_log2FC, decreasing = TRUE)

log_fold_changes_top_100_genes <- dge_list_full$avg_log2FC
names(log_fold_changes_top_100_genes) <- dge_list_full$entrezid



#convert to entrezids
library(org.Hs.eg.db) 
library(clusterProfiler)

gene.df <- bitr(dge_list, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

edo <- enrichDGN(gene.df$ENTREZID)

#barplot
library(enrichplot)
barplot(edo, showCategory=10) 

#dotplot(edo, showCategory=10) 
#heatmap
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

#treeplot
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

#dotplot split by activated and supressed
gene_list = sort(gene.df$SYMBOL, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


################ Universal enrichment analysis ################
library(msigdbr)
msigdbr_species()

#Retrieve all human gene sets
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

gene_list2 = sort(gene.df$ENTREZID, decreasing = TRUE)
em <- enricher(gene_list2, TERM2GENE=m_t2g)
head(em)

C7_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C7_t2g)

em2 <- GSEA(log_fold_changes_top_100_genes, TERM2GENE = C7_t2g)
head(em2)

y <- GSEA(gene_list2, TERM2GENE = cells)
head(y)

#plots
barplot(em, showCategory=10) 
dotplot(em, showCategory=10) + ggtitle("dotplot for ORA")
gseaplot2(em2, geneSetID = 1, title = em2$Description[1])
#gseaplot2(em2, geneSetID = 1, title = em2$Description[1])





################################## Annotate T cell subsets ###################################

#Genes from Gustavo's R-scripts
gn <- c('CD3E','CD4','CD8A', 'CD8B', 'CCR7', 'FCER1A', 'FCGR3A', 'NKG7', 'GZMB', 'GZMA', 
        'GNLY', 'MS4A1', 'CD14', 'LYZ', 'CD34', 'IL7R', 'PDCD1', 'B3GAT1', 'PECAM1', 'ICOS', 
        'IL2RA', 'ITGAE', 'CD69', 'CXCR6', 'ITGA1', 'PTPRC', 'CX3CR1', 'SELL', 'CD19', 'CXCR5', 
        'CCR6', 'FAS', 'KLRB1', 'CD27', 'CD28', 'TRDC', 'TRGC1', 'TRGC2', 'ITGB1', 'SELPLG', 
        'CCR2', 'CST3', 'MS4A7')#, '', '', '', '', '')
names(gn) <- c('T cell', 'CD4 T cell', 'CD8 T cell', 'CD8 T cell', 'Tn', 'DC', 'Monocyte', 'NK cell', 'NK cell', 'NK cell', 
               'NK cell', 'B cell', 'Monocyte', 'Monocyte', 'HSC', 'T cell', 'Tm', 'Tem', 'Tn', 'Treg', 
               'Treg', 'Trm', 'Trm', 'Trm', 'Trm', 'Monocyte', 'trafficking', 'Tn', 'B cell', 'B cell', 
               'B cell', 'Tm', 'RORgT cell', 'Tn and Tcm', 'Tn and Tcm', 'gdT cell', 'gdT cell', 'gdT cell', 'trafficking', 'trafficking', 
               'Monocyte', 'DC', 'Monocyte')

gn

#Dotplot
markers.to.plot <- c("TCF7", "CCR7", "PECAM1", "SELL", "CD27", "CD28", "PDCD1", "FAS", "B3GAT1", "ITGA1", "CXCR6", "CD69", "ITGAE", "LGALS3")
DotPlot(ab_Tcell_subset, features = markers.to.plot, cols = c("red", "blue", "grey"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()

















###### Old

t.cells <- seurat_Tcell_subset
Idents(t.cells) <- "fetal"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

Tn <- subset(seurat_Tcell_subset, idents = "Tn")
Idents(Tn) <- "fetal"
avg.Tn <- as.data.frame(log1p(AverageExpression(Tn, verbose = FALSE)$RNA))
avg.Tn$gene <- rownames(avg.Tn)

#Choose naive T cell sigs
genes.to.label = c("TCF7", "CCR7", "PECAM1", "SELL", "CD27", "CD28")
p1 <- ggplot(avg.t.cells, aes(skin-fetus, skin-adult)) + geom_point() + ggtitle("Tn")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.Tn, aes(fetal, adult)) + geom_point() + ggtitle("Tn and Tcm")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

#"TCF7", "CCR7", "PECAM1", "SELL"
FeaturePlot(seurat_Tcell_subset, features = c("TCF7", "CCR7", "PECAM1", "SELL"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(seurat_Tcell_subset, features = c("TCF7", "CCR7", "PECAM1", "SELL"), split.by = "orig.ident",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)



immune.combined$celltype.con <- paste(Idents(seurat_Tcell_subset), seurat_Tcell_subset$orig.ident, sep = "_")
seurat_Tcell_subset$celltype <- Idents(seurat_Tcell_subset)
Idents(immune.combined) <- "celltype.stim"
response <- FindMarkers(seurat_Tcell_subset, ident.1 = "Tn", verbose = FALSE)
head(response, n = 15)

#Find DEGs
DEG_across_condition <- FindMarkers(seurat_Tcell_subset, ident.1 = "Tn", grouping.var = "orig.ident", verbose = FALSE)
head(DEG_across_condition, n = 15)

#Number of DEGs in total
nrow(DEG_across_condition) # -> 971

#Subset to find up -and downregulated DEGs
downregulated_DEGs <- DEG_across_condition[DEG_across_condition$avg_log2FC < 0, ]
nrow(downregulated_DEGs) # -> 311
#add gene description in extra column
library(org.Hs.eg.db) 

# Add gene full name
downregulated_DEGs$description <- mapIds(x = org.Hs.eg.db,
                                         keys = row.names(downregulated_DEGs),
                                         column = "GENENAME",
                                         keytype = "SYMBOL",
                                         multiVals = "first")
head(downregulated_DEGs, n = 10)

upregulated_DEGs <- DEG_across_condition[DEG_across_condition$avg_log2FC > 0, ]
nrow(upregulated_DEGs) # -> 660
#add gene description in extra column
# Add gene full name
upregulated_DEGs$description <- mapIds(x = org.Hs.eg.db,
                                       keys = row.names(upregulated_DEGs),
                                       column = "GENENAME",
                                       keytype = "SYMBOL",
                                       multiVals = "first")
head(upregulated_DEGs, n = 10)

####################################################################################################

# What are the cell names of specific cells?
#WhichCells(ab_Tcell_subset, idents = "X")

# How can I extract expression matrix for all Tn cells 
#nk.raw.data <- as.matrix(GetAssayData(seurat_Tcell_subset, slot = "counts")[, WhichCells(seurat_Tcell_subset, ident = "Tn")])

# Can I create a Seurat object based on expression of a feature or value in object metadata?
#tcf7_seurat_object <- subset(seurat_Tcell_subset, subset = TCF7 > 1)

#subset(seurat_Tcell_subset, subset = orig.ident == "fetal")

# Can I create a Seurat object of just the Tn cells and Tn and Tcm cells?
#Seurat_only_Tn_and_Tn_Tcm <- subset(seurat_Tcell_subset, idents = c("Tn", "Tn and Tcm"))

################################## Residency markers #########################################

VlnPlot(immune.combined, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
FeaturePlot(immune.combined, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))

VlnPlot(seurat_Tcell_subset, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))
FeaturePlot(seurat_Tcell_subset, features = c("S1PR1", "SELL", "CD69", "LGALS3", "RUNX3", "FABP4", "FABP5", "ITGAE", "SELPLG", "CCR7", "PRDM1"))

################################# Average gene expression within a cluster ###################

# How can I calculate the average expression of all cells within a cluster?
cluster.averages <- AverageExpression(seurat_Tcell_subset)
head(cluster.averages[["RNA"]][, 1:5])

# First, replace spaces with underscores '_' so ggplot2 doesn't fail
orig.levels <- levels(seurat_Tcell_subset)
Idents(seurat_Tcell_subset) <- gsub(pattern = " ", replacement = "_", x = Idents(seurat_Tcell_subset))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(seurat_Tcell_subset) <- orig.levels
cluster.averages <- AverageExpression(seurat_Tcell_subset, return.seurat = TRUE)
cluster.averages


# Plot the average expression of Tn cells vs. Tn and Tcm cells?  Pass do.hover = T for an
# interactive plot to identify gene outliers
CellScatter(cluster.averages, cell1 = "Tn", cell2 = "Tn_and_Tcm")

# How can I calculate expression averages separately for each condition?
#cluster.averages <- AverageExpression(seurat_Tcell_subset, return.seurat = TRUE, add.ident = "condition")
#CellScatter(cluster.averages, cell1 = "Tn_con1", cell2 = "Tn_con2")

# Heatmap of these 'in silico' bulk datasets to visualize agreement between
# conditions
DoHeatmap(cluster.averages, features = unlist(TopFeatures(seurat_Tcell_subset[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)


