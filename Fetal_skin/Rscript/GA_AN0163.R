#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a transplanted
# AML patient (chimeric sample)
# Date: 2019.02.14
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger
# - dimensionality reduction, clustering and biological annotation of clusters
# - differential expression of each cluster
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(openxlsx)
#library(Seurat, lib.loc = paste0(.libPaths(), "/Seurat/Seurat_v3.0"))
library(Seurat)
library(dplyr)
library(SingleR)
library(ggplot2)
library(gridExtra)
#library(stringr)
#library(reshape2)
#library(plyr)
#library(gridExtra)
#library(scales)
#library(VennDiagram)
#library(dplyr)
#library(grid)
#library(readr)
#library(vsn)
#library(edgeR)
#library(biomaRt)
#library(pheatmap)
#library(RColorBrewer)
#library(ggrepel)
#library(ggplotify)

# Define custom theme for ggplot
#theme_set(theme_grey()) # cowplot (loaded along with Seurat v2.3.4 only) sets theme_set(theme_cowplot()) as default
theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = "black"),
                      axis.text.y = element_text(size = 12.8, color = "black"), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = "black"), 
                      axis.ticks.y = element_line(size = 0.4, colour = "black"),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour="black"),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = "black", size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid <- c("EX0008")

# Define Analysis_ID according to ngs_analysis_list
analysisid <- c("GA_AN0163")

# Define samples analyzed (all)
#sa <- c("cd")

# Define reference sample name
ct <- c("blood")
tr <- c("skin")

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c("sc_skin_blood_10x_chimeric")
an.descs <- c("sc_tp")

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx("../../ngs_sample_list.xlsm", sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = "NA", fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], summarise, sum = NA)

# Update table containing the analyses list to include this analysis
# Include info like which read alignment algorithm was used
# Load, update avoiding duplicates (when script is run several times) and save back overwriting old one
nal <- read.xlsx("../../ngs_analysis_list.xlsm", sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = "NA", fillMergedCells = F)
nal <- nal[nal$Analysis_ID == analysisid, ]

# Sample data
metadata <- data.frame(row.names = nsl$Sample_name, 
                       reshape2::colsplit(nsl$Sample_name, "_", c("exp", "celltype", "source")))

#-----------------------------------------------------------------------------------------------------------------#
# Archive old runs of this same analysis
#-----------------------------------------------------------------------------------------------------------------#

## Copy "script" file to "analysis_old/scripts_old" folder adding run date to its name
#file.copy(from = paste0(analysisid, ".R"), 
#          to = paste0("analysis_old/", analysisid, format(Sys.time(), "_%Y%m%d_%H%M%S"), ".R"))
#
## Move "report" file to "analysis_old" folder adding creation date to its name
#file.rename(from = paste0("", analysisid, ".html"), 
#            to = paste0("analysis_old/", analysisid, 
#                        format(file.info(paste0(analysisid, ".html"))$ctime, "_%Y%m%d_%H%M%S"), ".html"))
#
## Move "figures" folder to "analysis_old/figures_old" folder
#if(dir.exists(paste0("figures/", analysisid, "_", an.desc))) {
#  file.copy(from = paste0("figures/", analysisid, "_", an.desc), to = paste0("analysis_old/figures_old/"), 
#            recursive = T, copy.date = T) }
#file.rename(from = paste0("analysis_old/figures_old/", analysisid, "_", an.desc), 
#            to = paste0("analysis_old/figures_old/", analysisid, "_", an.desc, 
#                        format(file.info(paste0("analysis_old/figures_old/", 
#                                                analysisid, "_", an.desc))$ctime, "_%Y%m%d_%H%M%S")))
#unlink(paste0("figures/", analysisid, "_", an.desc), recursive = T)

# Create new folder to store figures from this analysis
dir.create(paste0("figures/", analysisid, "_", an.desc))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/pca"))
dir.figs <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
dir.figs.bd <- paste0("figures/", analysisid, "_", an.desc, "/blood/", analysisid, "_", an.descs, "_blood")
dir.figs.sk <- paste0("figures/", analysisid, "_", an.desc, "/skin/", analysisid, "_", an.descs, "_skin")
#dir.figs.pca <- paste0("figures/", analysisid, "_", an.desc, "/pca/", analysisid, "_", an.descs)

#=================================================================================================================#
# Analysis of processed data for blood sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Clustering and differential expression analysis with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Load UMI counts from 10x data
bd <- Read10X(data.dir = "../../procdata/EX0008/GA_PD0008_cellranger/tp_blood_1/filtered_feature_bc_matrix/")
#bd <- readRDS(paste0(dir.figs.bd, "_seurat.rds"))

# Initialize the Seurat object with the raw (non-normalized data).
bd <- CreateSeuratObject(counts = bd, min.cells = 0, min.features = 0, project = "blood")

# Calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`
mito.features <- grep(pattern = "^MT-", x = rownames(bd), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = bd, slot = 'counts')[mito.features, ]) / 
  Matrix::colSums(GetAssayData(object = bd, slot = 'counts'))
bd[['percent.mito']] <- percent.mito
VlnPlot(bd, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Visualize feature-feature relationships
FeatureScatter(bd, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black", group.by = "orig.ident")
FeatureScatter(bd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black", group.by = "orig.ident")

# Filter out cells that have unique feature counts over 2000 or less than 200
bd <- subset(bd, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 0.2)
#set.seed(1); bd <- subset(bd, cells = sample(colnames(bd), size = 1000)) #length(colnames(bd)); head(colnames(bd),3)
vp.bd <- VlnPlot(bd, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# Employ a global-scaling normalization method that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10000 by default), and log-transforms the result
bd <- NormalizeData(bd, normalization.method = "LogNormalize", scale.factor = 1e4)

# Detection of variable features across the single cells
# Calculate the average expression and dispersion for each feature and z-score for dispersion within each bin
#bd <- FindVariableFeatures(bd, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), 
#                           dispersion.cutoff = c(0.5, Inf))
bd <- FindVariableFeatures(bd, selection.method = "vst", nfeatures = 8000, 
                           mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(bd))
vf.bd <- VariableFeaturePlot(bd, cols = c("black", "red"), pt.size = 1, log = NULL, assay = NULL) + theme_bw() + theme_custom

# Scaling the data and removing unwanted sources of variation (stored in the scale.data slot)
bd <- ScaleData(bd, features = rownames(bd), vars.to.regress = c("nCount_RNA", "percent.mito"))

# Perform linear dimensional reduction
bd <- RunPCA(bd, features = VariableFeatures(bd), verbose = F)

# Examine and visualize PCA results a few different ways (VizDimReduction, DimPlot, and DimHeatmap)
print(bd[['pca']], dims = 1:5, nfeatures = 5, projected = F)
VizDimLoadings(bd, dims = c(1,2))
DimPlot(bd, reduction = "pca")

# Scores each feature in the dataset (including features not included in the PCA) based on their correlation with the PCs
bd <- ProjectDim(bd, do.center = T, verbose = F)
print(bd[['pca']], dims = 1:5, nfeatures = 5, projected = T)
VizDimLoadings(bd, dims = c(1,2), projected = T)
DimPlot(bd, projected = T, reduction = "pca")

# Explore the primary sources of heterogeneity to decide which PCs to include for further downstream analyses
DimHeatmap(bd, dims = 1, cells = NULL, balanced = T)
#DimHeatmap(bd, dims = 1:20, cells = 500, balanced = T)
hm.pc.bd <- DimHeatmap(bd, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle("Heatmap of PC markers from 500 cells (blood)")

# Determine statistically significant principal components
#bd <- JackStraw(bd, dims = 20, prop.freq = 0.01, num.replicate = 100, maxit = 1000)
#bd <- ScoreJackStraw(bd, dims = 1:20)
#JackStrawPlot(bd, dims = 1:20)
ep.bd <- ElbowPlot(bd, ndims = 25)

# Cluster cells (the ‘granularity’ (number of clusters) of the downstream clustering is set to 0.6-1.2 for 3K cells)
bd <- FindNeighbors(bd, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
#bd <- FindClusters(bd, resolution = 1, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)
bd <- FindClusters(bd, resolution = 1.56, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)
bd <- FindClusters(bd, resolution = 1.51, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Run Non-linear dimensional reduction (tSNE and UMAP)
bd <- RunTSNE(bd, reduction = "pca", dims = 1:20, seed.use = 10, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30)
bd <- RunUMAP(bd, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
DimPlot(bd, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.51", split.by = "orig.ident", label = T, label.size = 8, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom
DimPlot(bd, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1.51", split.by = "orig.ident", label = T, label.size = 8, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom
#saveRDS(bd, file = paste0(dir.figs.bd, "_seurat.rds"))
bd[["clust_1"]] <- Idents(object = bd)
ts.sp.bd <- DimPlot(bd, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.51", split.by = "clust_1", 
                    label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_y_continuous(breaks = c(-25,0,25)) + theme_bw() + theme_custom
vp.bdc <- VlnPlot(bd, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "RNA_snn_res.1.51")

# Finding differentially expressed features (cluster biomarkers)
bd.markers <- FindAllMarkers(bd, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                             test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
#bd.markers <- bd.markers[order(bd.markers$cluster, decreasing = F),]
#write.csv(bd.markers, file = paste0(dir.figs.bd, "_cluster_markers.csv"), row.names = T)
#bd.markers <- read.csv(file = paste0(dir.figs.bd, "_cluster_markers.csv"), row.names = 1)
bd.top <- bd.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gn <- c("IL7R", "CD4", "CD8A", "CD8B", "GNLY", "MS4A1", "FCGR3A", "NKG7", "GZMB")

# Distinguish naive and memory
FeaturePlot(object = bd, features = c("S100A4", "CCR7"), cols = c("lightgrey", "blue"), coord.fixed = T)

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
VlnPlot(bd, features = bd.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
VlnPlot(bd, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(bd, features = gn, dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "tsne", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.bd <- DoHeatmap(bd, features = bd.top$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers (blood)")# + NoLegend()

#=================================================================================================================#
# Analysis of processed data for skin sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Clustering and differential expression analysis with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Run workflow described above for skin sample
sk <- Read10X(data.dir = "../../procdata/EX0008/GA_PD0008_cellranger/tp_skin_1/filtered_feature_bc_matrix/")
#sk <- readRDS(paste0(dir.figs.sk, "_seurat.rds"))
sk <- CreateSeuratObject(counts = sk, min.cells = 0, min.features = 0, project = "skin")
sk[['percent.mito']] <- Matrix::colSums(GetAssayData(sk, slot = 'counts')[grep(pattern = "^MT-", rownames(sk), value = T), ]) / 
  Matrix::colSums(GetAssayData(sk, slot = 'counts'))
sk <- subset(sk, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 0.2)
#sk$celltype <- sr.sk$meta.data$celltype
#sk$geno.ident <- sr.sk$meta.data$geno.ident
##sk <- subset(sk, subset = celltype != "Class-switched memory B-cells" & celltype != "Memory B-cells" & celltype != "naive B-cells")
#sk <- subset(sk, subset = RNA_snn_res.1.8 != 13)
sk <- NormalizeData(sk, normalization.method = "LogNormalize", scale.factor = 1e4)
sk <- FindVariableFeatures(sk, selection.method = "vst", nfeatures = 8000, 
                           mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
sk <- ScaleData(sk, features = rownames(sk), vars.to.regress = c("nCount_RNA", "percent.mito"))
sk <- RunPCA(sk, features = VariableFeatures(sk), verbose = F)
sk <- ProjectDim(sk, do.center = T, verbose = F)
#sk <- JackStraw(sk, dims = 20, prop.freq = 0.01, num.replicate = 100, maxit = 1000)
#sk <- ScoreJackStraw(sk, dims = 1:20)
#JackStrawPlot(sk, dims = 1:20)
sk <- FindNeighbors(sk, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
sk <- FindClusters(sk, resolution = 1.2, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)
#sk <- FindClusters(sk, resolution = 1.8, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)
sk@meta.data$clust_1 <- sk@meta.data$RNA_snn_res.1.2
sk <- RunTSNE(sk, reduction = "pca", dims = 1:20, seed.use = 11, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30, reduction.name = "tsne")
#sk <- RunTSNE(sk, reduction = "pca", dims = 1:20, seed.use = 12, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30, reduction.name = "tsne")
sk <- RunUMAP(sk, reduction = "pca", dims = 1:20, seed.use = 1, reduction.name = "umap", spread = 10, min.dist = 0.001)
#sk <- RunUMAP(sk, reduction = "pca", dims = 1:20, seed.use = 1, reduction.name = "umap")
sk.markers <- FindAllMarkers(sk, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                             test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
sk.top <- sk.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(sk.markers, file = paste0(dir.figs.sk, "_skin_cluster_markers.csv"), row.names = T)
#sk.markers <- read.csv(file = paste0(dir.figs.sk, "_skin_cluster_markers.csv"), row.names = 1)
#saveRDS(sk, file = paste0(dir.figs.sk, "_seurat.rds"))

# Overview of the data
FeatureScatter(sk, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black", group.by = "orig.ident")
FeatureScatter(sk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black", group.by = "orig.ident")
vp.skc <- VlnPlot(sk, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "RNA_snn_res.1.2")
vp.sk <- VlnPlot(sk, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")
length(VariableFeatures(sk))
vf.sk <- VariableFeaturePlot(sk, cols = c("black", "red"), pt.size = 1, log = NULL, assay = NULL) + theme_bw() + theme_custom
hm.pc.sk <- DimHeatmap(sk, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle("Heatmap of PC markers from 500 cells (skin)")
#JackStrawPlot(sk, dims = 1:20)
ep.sk <- ElbowPlot(sk, ndims = 25)
VlnPlot(sk, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(sk, features = sk.top$gene[1:6], dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "tsne", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = F)
FeaturePlot(sk, features = gn, dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "tsne", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.sk <- DoHeatmap(sk, features = sk.top$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers (skin)")# + NoLegend()
DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = "orig.ident", label = T, label.size = 8, pt.size = 1, cols = NULL) + theme_bw() + theme_custom
DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = "orig.ident", label = T, label.size = 8, pt.size = 1, cols = NULL) + theme_bw() + theme_custom
ts.sp.sk <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = "clust_1", 
                    label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_x_continuous(breaks = c(-20,0,20)) + theme_bw() + theme_custom

#=================================================================================================================#
# Cluster annotation
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Biological annotation of clusters from blood sample using SingleR
#-----------------------------------------------------------------------------------------------------------------#

# Create the SingleR object using clusters calculated with Seurat
sr.bd <- CreateSinglerObject(GetAssayData(bd, slot = "data"), 
                             annot = NULL, project.name = "tp_blood_ann", min.genes = 200, #500
                             technology = "10X", species = "Human", citation = "",
                             ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                             fine.tune = T, do.signatures = T, clusters = bd@active.ident, do.main.types = T, 
                             reduce.file.size = T, numCores = 4)#SingleR.numCores)
#saveRDS(sr.bd, file = paste0(dir.figs.bd, "_singler.rds"))
#sr.bd <- readRDS(paste0(dir.figs.bd, "_singler.rds"))

# Create the SingleR object using clusters calculated with Seurat using more variable geness
hpca$sd.thres <- sort(matrixStats::rowSds(hpca$data), decreasing = T)[10000]
blueprint_encode$sd.thres <- sort(matrixStats::rowSds(blueprint_encode$data), decreasing = T)[12000]
sr.bd <- CreateSinglerObject(GetAssayData(bd, slot = "data"), 
                             annot = NULL, project.name = "tp_blood_ann", min.genes = 200, #500
                             technology = "10X", species = "Human", citation = "",
                             ref.list = list(), normalize.gene.length = F, variable.genes = "sd",
                             fine.tune = T, do.signatures = T, clusters = bd@active.ident, do.main.types = T, 
                             reduce.file.size = T, numCores = 4)#SingleR.numCores)
#saveRDS(sr.bd, file = paste0(dir.figs.bd, "_singler_sd.rds"))

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
                                                                   "mB cell", "mB cell", "nB cell", "NK cell", "other", "Treg"))

# Correlation of variable genes between the experiment and the reference dataset
SingleR.DrawScatter(sc_data = GetAssayData(bd, slot = "data"), cell_id = 1, ref = hpca, sample_id = 1)

# Correlation between a single cell from the experiment and all reference cells
SingleR.DrawBoxPlot(sc_data = GetAssayData(bd, slot = "data"), cell_id = 1, ref = blueprint_encode, 
                    main_types = T, labels.use = NULL)$p

# Number of cells per annotated cluster
sr.bdn1 <- table(sr.bd$meta.data$orig.ident, sr.bd$meta.data$clusters); knitr::kable(sr.bdn1)
sr.bdn2 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn2)
sr.bdn3 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$clusters); knitr::kable(sr.bdn3)
sr.bdn4 <- table(sr.bd$singler[[2]]$SingleR.single.main$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn4)

# Heatmap of the aggregated scores before fine-tuning for the main cell types:
SingleR.DrawHeatmap(sr.bd$singler[[2]]$SingleR.single, top.n = Inf, clusters = sr.bd$meta.data$cluster)#orig.ident)

# Define colors
levels(factor(sr.bd$meta.data$clusters))
levels(factor(sr.bd$singler[[2]]$SingleR.single$labels[,1]))
#col <- c("darkgoldenrod1", "steelblue1", "wheat2", "firebrick1", "plum3", "darkorchid2", "royalblue1", "darkgoldenrod1", 
#         "plum3", "darkgoldenrod1", "seagreen3", "grey70", "olivedrab")
col <- c("darkorange2", "steelblue1", "wheat2", "firebrick1", "royalblue1", "darkorchid2", "orchid1", "darkgoldenrod1", 
         "grey70", "plum3", "seagreen3", "grey40", "olivedrab")
#col2 <- c("darkorange2", "firebrick1", "darkgoldenrod1", "steelblue1", "royalblue1", "plum3", "seagreen3", 
#          "grey70", "wheat2", "orchid1", "olivedrab", "darkorchid2", "grey40")
#col2 <- c("darkgoldenrod1", "firebrick1", "grey70", "darkorange2", "royalblue1", "plum3", "grey40", 
#          "darkorchid2", "wheat2", "seagreen3", "orchid1", "olivedrab")
#col2 <- c("darkorange2", "firebrick1", "darkgoldenrod1", "steelblue1", "royalblue1", "plum3", "seagreen3", "grey40", 
#          "grey70", "wheat2", "orchid1", "olivedrab", "darkorchid2")
col2 <- c("firebrick1", "grey70", "darkgoldenrod1", "royalblue1", "plum3", "steelblue1", "seagreen3", 
          "wheat2", "darkorchid2", "black", "olivedrab")
col3 <- c("darkgoldenrod1", "firebrick1", "wheat2", "olivedrab", "firebrick1", "grey70", "steelblue1", 
          "darkorchid2", "royalblue1", "plum3", "seagreen3", "darkorchid2", "darkgoldenrod1")

# Draw tSNE plots swith cluster from Seurat object
ts.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col, 
                          labels = sr.bd$meta.data$clusters, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  scale_color_manual("Subset:", values = col, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25))
# Draw tSNE plots with annotation from SingleR at single cell level
ta.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col2, 
                          labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  scale_color_manual("Subset:", values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25))
ta.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = F, do.letters = F, col = col2, 
                          labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  scale_color_manual("Subset:", values = col2) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw tSNE plots with annotation from SingleR at cluster level
tc.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col3, 
                          clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
                          labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  scale_color_manual("Subset:", values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25))
tc.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy, do.labels = F, do.letters = F, col = col3, 
                          clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
                          labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  scale_color_manual("Subset:", values = col3) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw UMAP plots swith cluster from Seurat object
um.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col, 
                          labels = sr.bd$meta.data$clusters, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-50,25,25))
# Draw UMAP plots with annotation from SingleR at single cell level
ua.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col2, 
                          labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-50,25,25))
ua.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = F, do.letters = F, col = col2, 
                             labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col2) + theme_custom + scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-50,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw UMAP plots with annotation from SingleR at cluster level
uc.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col3, 
                          clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
                          labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-50,25,25))
uc.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy.um, do.labels = F, do.letters = F, col = col3, 
                             clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
                             labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = "hsct blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col3) + theme_custom + scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-50,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))

# Visualization of marker gene expression
df.bd <- data.frame(x = sr.bd$meta.data$xy[,1],
                 y = sr.bd$meta.data$xy[,2],
                 t(as.matrix(GetAssayData(bd, slot = "data")[c('CD3E','CD4','CD8A', 'CD8B', 'CCR7', 'FCER1A', 'PPBP', 
                                                               'FCGR3A', 'NKG7', 'GZMB', 'GZMA', 
                                                               'GNLY','MS4A1','CD14', 'LYZ', 'S100A4'),])))
df.bd <- reshape2::melt(df.bd, id.vars = c('x','y'))
ts.mk.bd <- ggplot(df.bd, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
  scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
  scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~variable, ncol = 4) + theme_custom + xlab('') + ylab('')

# Confidence of annotation
sc <- bd # because function has a problem in calling seurat object, and it calls sc instead of the second argument
SingleR.PlotFeature(sr.bd$singler[[2]]$SingleR.single, bd, plot.feature = 'MaxScore', dot.size = 2)
SingleR.PlotFeature(sr.bd$singler[[2]]$SingleR.single, bd,
                    plot.feature = -log10(sr.bd$singler[[2]]$SingleR.single$pval))
bd.cf <- data.frame(nGene = bd[["nFeature_RNA"]], 
                    x = sr.bd$meta.data$xy[,1],
                    y = sr.bd$meta.data$xy[,2], 
                    pval = -log10(sr.bd$singler[[2]]$SingleR.single$pval),
                    Identity = sr.bd$singler[[2]]$SingleR.single$labels)
bd.pva <- ggplot(bd.cf, aes(Identity, y = pval, color = Identity)) + geom_boxplot() + scale_color_manual("Subset:", values = col2) + 
  scale_x_discrete("") + scale_y_continuous(limits = c(0.5,2.5)) + ggtitle("hsct blood") + 
  ylab("-log10(p-value)") + theme_custom + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = NULL)
bd.pvt <- ggplot(bd.cf, aes(x = x, y = y, color = pval)) + geom_point(size = 0.8) + 
  scale_color_gradientn(expression("-"*log[10]*"(p-value)"), colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(50)[c(1:22,29:50)], 
                        limits = c(0.5,2), oob = scales::squish, values = c(0,(1.3-0.5)/2,1)) + #c(0,-log10(0.1)/(max(bd.cf$pval)-min(bd.cf$pval)),1)) +
  #scale_color_gradient2("Subset:", low = "magenta", mid = "black", high = scales::muted("yellow",90,100),  
  #                      midpoint = 1.2, limits = c(0.5,2), oob = scales::squish) + 
  scale_x_continuous("tSNE1", breaks = seq(-20,20,20)) + scale_y_continuous("tSNE2", breaks = seq(-25,25,25)) +
  theme_custom + ggtitle("hsct blood") #+ guides(color = guide_legend(title.position = "left"), title.theme = element_text(angle = 180))

#-----------------------------------------------------------------------------------------------------------------#
# Biological annotation of clusters from skin sample using SingleR
#-----------------------------------------------------------------------------------------------------------------#

# Create the SingleR object using clusters calculated with Seurat
sr.sk <- CreateSinglerObject(GetAssayData(sk, slot = "data"), 
                             annot = NULL, project.name = "tp_skin_ann", min.genes = 200, #500
                             technology = "10X", species = "Human", citation = "",
                             ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                             fine.tune = T, do.signatures = T, clusters = sk@active.ident, do.main.types = T, 
                             reduce.file.size = T, numCores = 4)#SingleR.numCores)
#saveRDS(sr.sk, file = paste0(dir.figs.sk, "_singler.rds"))
#sr.sk <- readRDS(paste0(dir.figs.sk, "_singler.rds"))

# Create the SingleR object using clusters calculated with Seurat using more variable genes
sr.sk <- CreateSinglerObject(GetAssayData(sk, slot = "data"), 
                             annot = NULL, project.name = "tp_skin_ann", min.genes = 200, #500
                             technology = "10X", species = "Human", citation = "",
                             ref.list = list(), normalize.gene.length = F, variable.genes = "sd",
                             fine.tune = T, do.signatures = T, clusters = sk@active.ident, do.main.types = T, 
                             reduce.file.size = T, numCores = 2)#SingleR.numCores)
#saveRDS(sr.sk, file = paste0(dir.figs.sk, "_singler_sd.rds"))

# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
#sr.sk$seurat = seurat.object # (optional)
sr.sk$meta.data$orig.ident <- sk@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.sk$meta.data$xy <- sk@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "tsne"))
sr.sk$meta.data$xy.um <- sk@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "umap"))
sr.sk$meta.data$clusters <- sk@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(sk))
sr.sk$meta.data$celltype <- sr.sk$singler[[2]]$SingleR.single$labels
sr.sk$singler[[2]]$SingleR.single$labels <- plyr::mapvalues(sr.sk$singler[[2]]$SingleR.single$labels, 
                                                            from = levels(factor(sr.sk$singler[[2]]$SingleR.single$labels)), 
                                                            to = c("CD4 Tn", "CD4 Tcm", "CD4 Tem", 
                                                                   "CD8 Tn", "CD8 Tcm", "CD8 Tem", 
                                                                   "other", "other", "other", "other", "other", "other", "NK cell", "Treg"))

# Correlation of variable genes between the experiment and the reference dataset
SingleR.DrawScatter(sc_data = GetAssayData(sk, slot = "data"), cell_id = 1, ref = hpca, sample_id = 1)

# Correlation between a single cell from the experiment and all reference cells
SingleR.DrawBoxPlot(sc_data = GetAssayData(sk, slot = "data"), cell_id = 1, ref = blueprint_encode, 
                    main_types = T, labels.use = NULL)$p

# Number of cells per annotated cluster
sr.skn1 <- table(sr.sk$meta.data$orig.ident, sr.sk$meta.data$clusters); knitr::kable(sr.skn1)
sr.skn2 <- table(sr.sk$singler[[2]]$SingleR.single$labels, sr.sk$meta.data$orig.ident); knitr::kable(sr.skn2)
sr.skn3 <- table(sr.sk$singler[[2]]$SingleR.single$labels, sr.sk$meta.data$clusters); knitr::kable(sr.skn3)
sr.skn4 <- table(sr.sk$singler[[2]]$SingleR.single.main$labels, sr.sk$meta.data$orig.ident); knitr::kable(sr.skn4)

# Heatmap of the aggregated scores before fine-tuning for the main cell types:
SingleR.DrawHeatmap(sr.bd$singler[[2]]$SingleR.single, top.n = Inf, clusters = sr.bd$meta.data$cluster)#orig.ident)

# Define colors
levels(factor(sr.sk$meta.data$clusters))
levels(factor(sr.sk$singler[[2]]$SingleR.single$labels[,1]))
col4 <- c("darkorange2", "firebrick1", "darkgoldenrod1", "steelblue1", "grey70", "plum3", "darkorchid2", 
          "seagreen3", "grey40", "royalblue1", "olivedrab", "orchid1", "wheat2")
col5 <- c("firebrick1", "grey70", "darkgoldenrod1", "royalblue1", "plum3", "steelblue1", "darkorchid2", "black", "olivedrab")
col6 <- c("darkorange2", "darkgoldenrod1", "darkgoldenrod1", "steelblue1", "royalblue1", "darkgoldenrod1", "darkorange2", 
          "royalblue1", "darkgoldenrod1", "steelblue1", "darkorange2", "firebrick1", "olivedrab")

# Draw tSNE plots with cluster from Seurat object
ts.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy, do.label = T, col = col4, do.letters = F, 
                          labels = sr.sk$meta.data$clusters, label.size = 4, dot.size = 1.5,  title = "hsct skin")$p + 
  scale_color_manual("Subset:", values = col4, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20))
# Draw tSNE plots with annotation from SingleR at single cell level
ta.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy, do.label = T, col = col5, 
                          do.letters = F, labels = sr.sk$singler[[2]]$SingleR.single$labels, label.size = 4, 
                          dot.size = 1.5, title = "hsct skin")$p + 
  scale_color_manual("Subset:", values = col5, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20))
ta.sk.nl <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy, do.label = F, col = col5, 
                             do.letters = F, labels = sr.sk$singler[[2]]$SingleR.single$labels, label.size = 4, 
                             dot.size = 1.5, title = "hsct skin")$p + 
  scale_color_manual("Subset:", values = col5) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw tSNE plots with annotation from SingleR at cluster level
tc.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.clusters, sr.sk$meta.data$xy, do.label = T, col = col6, 
                          clusters = sr.sk$meta.data$clusters, 
                          do.letters = F, labels = sr.sk$singler[[2]]$SingleR.clusters$labels, label.size = 4, 
                          dot.size = 1.5, title = "hsct skin")$p + 
  scale_color_manual("Subset:", values = col6, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20))
tc.sk.nl <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.clusters, sr.sk$meta.data$xy, do.label = F, col = col6, 
                             clusters = sr.sk$meta.data$clusters, 
                             do.letters = F, labels = sr.sk$singler[[2]]$SingleR.clusters$labels, label.size = 4, 
                             dot.size = 1.5, title = "hsct skin")$p + 
  scale_color_manual("Subset:", values = col6) + theme_custom + scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw UMAP plots with cluster from Seurat object
um.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy.um, do.label = T, col = col4, 
                          do.letters = F, labels = sr.sk$meta.data$clusters, label.size = 4, dot.size = 1.5, title = "hsct skin")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + #scale_x_continuous(breaks = seq(-4,4,4)) + scale_y_continuous(breaks = seq(-4,8,4)) +  
  scale_color_manual("Subset:", values = col4, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25))
# Draw UMAP plots with annotation from SingleR at single cell level
ua.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy.um, do.label = T, col = col5, 
                          do.letters = F, labels = sr.sk$singler[[2]]$SingleR.single$labels, label.size = 4, 
                          dot.size = 1.5, title = "hsct skin")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col5, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25))
ua.sk.nl <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.single, sr.sk$meta.data$xy.um, do.label = F, col = col5, 
                             do.letters = F, labels = sr.sk$singler[[2]]$SingleR.single$labels, label.size = 4, 
                             dot.size = 1.5, title = "hsct skin")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col5) + theme_custom + scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))
# Draw UMAP plots with annotation from SingleR at cluster level
uc.sk <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.clusters, sr.sk$meta.data$xy.um, do.label = T, col = col6, 
                          clusters = sr.sk$meta.data$clusters, 
                          do.letters = F, labels = sr.sk$singler[[2]]$SingleR.clusters$labels, label.size = 4, 
                          dot.size = 1.5, title = "hsct skin")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col6, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom + scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25))
uc.sk.nl <- SingleR.PlotTsne(sr.sk$singler[[2]]$SingleR.clusters, sr.sk$meta.data$xy.um, do.label = F, col = col6, 
                             clusters = sr.sk$meta.data$clusters, 
                             do.letters = F, labels = sr.sk$singler[[2]]$SingleR.clusters$labels, label.size = 4, 
                             dot.size = 1.5, title = "hsct skin")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual("Subset:", values = col6) + theme_custom + scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8), title.theme = element_text(hjust = 0)))

# Visualization of marker gene expression
df.sk <- data.frame(x = sr.sk$meta.data$xy[,1],
                    y = sr.sk$meta.data$xy[,2],
                    t(as.matrix(GetAssayData(sk, slot = "data")[c('CD3E','CD4','CD8A', 'CD8B', 'CCR7', 'FCER1A', 'PPBP', 
                                                                  'FCGR3A', 'NKG7', 'GZMB', 'GZMA', 
                                                                  'GNLY','MS4A1','CD14', 'LYZ', 'S100A4'),])))
df.sk <- reshape2::melt(df.sk, id.vars = c('x','y'))
ts.mk.sk <- ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
  scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
  scale_x_continuous(breaks = c(-20,0,20)) + facet_wrap(~variable, ncol = 4) + theme_custom + xlab('') + ylab('')

# Confidence of annotation
sc <- sk # because function has a problem in calling seurat object, and it calls sc instead of the second argument
SingleR.PlotFeature(sr.sk$singler[[2]]$SingleR.single, sk, plot.feature = 'MaxScore', dot.size = 2)
SingleR.PlotFeature(sr.sk$singler[[2]]$SingleR.single, sk,
                    plot.feature = -log10(sr.sk$singler[[2]]$SingleR.single$pval))
sk.cf <- data.frame(nGene = sk[["nFeature_RNA"]], 
                    x = sr.sk$meta.data$xy[,1],
                    y = sr.sk$meta.data$xy[,2], 
                    pval = -log10(sr.sk$singler[[2]]$SingleR.single$pval),
                    Identity = sr.sk$singler[[2]]$SingleR.single$labels)
sk.pva <- ggplot(sk.cf, aes(Identity, y = pval, color = Identity)) + geom_boxplot() + scale_color_manual("Subset:", values = col5) + 
  scale_x_discrete("") + scale_y_continuous(limits = c(0.5,2.5)) + ggtitle("hsct skin") + 
  ylab("-log10(p-value)") + theme_custom + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = NULL)
sk.pvt <- ggplot(sk.cf, aes(x = x, y = y, color = pval)) + geom_point(size = 0.8) + 
  scale_color_gradientn(expression("-"*log[10]*"(p-value)"), colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(50)[c(1:21,29:50)], 
                        limits = c(0.5,2), oob = scales::squish, values = c(0,(1.3-0.5)/2,1)) + #c(0,-log10(0.05)/(max(sk.cf$pval)-min(max(sk.cf$pval))),1)) +
  #scale_color_gradient2("Subset:", low = "magenta", mid = "black", high = scales::muted("yellow",90,100),  
  #                      midpoint = 1.2, limits = c(0.5,2), oob = scales::squish) + 
  scale_x_continuous("tSNE1", breaks = seq(-20,20,20)) + scale_y_continuous("tSNE2", breaks = seq(-20,20,20)) +
  ylab("-log10(p-value)") + theme_custom + ggtitle("hsct skin")

#-----------------------------------------------------------------------------------------------------------------#
# Distinction of host and donor cells
#-----------------------------------------------------------------------------------------------------------------#

# Load barcodes from host genotype
bc <- read.delim("figures/GA_AN0162_sc_skin_blood_genotyping/GA_AN0162_sc_skin_genotyping_barcodes_host.tsv",as.is = T,header = F)[,1]
# Keep only the barcodes present in seurat object after fitlering steps
bc <- gsub("-1", "", bc)
bc <- bc[bc %in% names(Idents(sk))]

# Classify idents according to barcodes from host and donor genotypes
#Idents(sk) <- names(Idents(sk))
#Idents(sk) <- ifelse(names(Idents(sk)) %in% bc, "host", "donor")
sk$geno.ident <- names(Idents(sk))
sk$geno.ident <- ifelse(sk$geno.ident %in% bc, "host", "donor")
col8 <- c("grey70", "royalblue2")

# Plot genotypes sepparately
#ts.dn <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = NULL, label = T, label.size = 5, pt.size = 1, cols = NULL, cells = colnames(sk)[!(colnames(sk) %in% bc)]) + 
#  xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("skin (donor)") + theme_bw() + theme_custom
#col7 <- unique(ggplot_build(ts.dn)$data[[1]][order(ggplot_build(ts.dn)$data[[1]]$group, decreasing = F), ]$colour)
col7 <- c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC")
ts.ht <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = NULL, label = T, label.size = 5, pt.size = 1, cols = col7, cells = bc) + 
  scale_y_continuous(limits = c(-35,35), expand = c(0,0)) + scale_x_continuous(limits = c(-31,32), expand = c(0,0)) + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("hsct skin (host)") + theme_bw() + theme_custom
um.ht <- DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = NULL, label = T, label.size = 5, pt.size = 1, cols = col7, cells = bc) + 
  scale_y_continuous(limits = c(-48,43), expand = c(0,0)) + scale_x_continuous(limits = c(-60,40), expand = c(0,0)) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("hsct skin (host)") + theme_bw() + theme_custom
ts.dn <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = NULL, label = T, label.size = 5, pt.size = 1, cols = col7, cells = colnames(sk)[!(colnames(sk) %in% bc)]) + 
  scale_y_continuous(limits = c(-35,35), expand = c(0,0)) + scale_x_continuous(limits = c(-31,32), expand = c(0,0)) + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("hsct skin (donor)") + theme_bw() + theme_custom
um.dn <- DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1.2", split.by = NULL, label = T, label.size = 5, pt.size = 1, cols = col7, cells = colnames(sk)[!(colnames(sk) %in% bc)]) + 
  scale_y_continuous(limits = c(-48,43), expand = c(0,0)) + scale_x_continuous(limits = c(-60,40), expand = c(0,0)) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("hsct skin (donor)") + theme_bw() + theme_custom
ts.ge <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "geno.ident", split.by = NULL, label = F, label.size = 5, pt.size = 1, cols = col8) + 
  scale_y_continuous(limits = c(-35,35), expand = c(0,0)) + scale_x_continuous(limits = c(-31,32), expand = c(0,0)) + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("hsct skin (genotyped)") + theme_bw() + theme_custom
um.ge <- DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "geno.ident", split.by = NULL, label = F, label.size = 5, pt.size = 1, cols = col8) + 
  scale_y_continuous(limits = c(-48,43), expand = c(0,0)) + scale_x_continuous(limits = c(-60,40), expand = c(0,0)) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("hsct skin (genotyped)") + theme_bw() + theme_custom

# Plot with ggplot2
df.sk2 <- data.frame(tx = sk@reductions$tsne@cell.embeddings[,1], 
                     ty = sk@reductions$tsne@cell.embeddings[,2], 
                     ux = sk@reductions$umap@cell.embeddings[,1], 
                     uy = sk@reductions$umap@cell.embeddings[,2], 
                     cells = colnames(sk), 
                     donor = sk@meta.data$geno.ident, 
                     celltype = sr.sk$singler[[2]]$SingleR.single$labels[,1])#, t(as.matrix(GetAssayData(sk, slot = "data")[genelist,])))
df.sk2 <- reshape2::melt(df.sk2, id.vars = c("tx", "ty", "ux", "uy", "cells", "donor", "celltype"))
#df.mk[[k]]$value <- df.mk[[k]]$value + rnorm(n = length(df.mk[[k]]$value)) / 100000 #add noise
ts.ge2 <- ggplot(df.sk2, aes(x = tx, y = ty, color = donor)) + geom_point(shape = 16, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey60", "royalblue2"), guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("hsct skin (genotyped through SNPs)"))
ts.ge2
ts.ge3 <- ggplot(df.sk2 %>% arrange(desc(donor)), aes(x = tx, y = ty, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey40", "black"), guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
ts.ge3
ts.ge4 <- ggplot(df.sk2 %>% arrange(desc(donor)), aes(x = tx, y = ty, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 1) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
ts.ge4
ts.ge5 <- ggplot(df.sk2 %>% arrange(donor), aes(x = tx, y = ty, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 1) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
ts.ge5
ts.ge6 <- ggplot(df.sk2 %>% arrange(donor), aes(x = tx, y = ty, color = donor, fill = donor)) + geom_point(shape = 21, size = 0.8, alpha = 1, stroke = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
ts.ge6
um.ge2 <- ggplot(df.sk2, aes(x = ux, y = uy, color = donor)) + geom_point(shape = 16, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey60", "royalblue2"), guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("hsct skin (genotyped through SNPs)"))
um.ge2
um.ge3 <- ggplot(df.sk2 %>% arrange(desc(donor)), aes(x = ux, y = uy, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey40", "black"), guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
um.ge3
um.ge4 <- ggplot(df.sk2 %>% arrange(desc(donor)), aes(x = ux, y = uy, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 1) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
um.ge4
um.ge5 <- ggplot(df.sk2 %>% arrange(donor), aes(x = ux, y = uy, color = donor, fill = donor)) + geom_point(shape = 21, size = 1.5, alpha = 1) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
um.ge5
um.ge6 <- ggplot(df.sk2 %>% arrange(donor), aes(x = ux, y = uy, color = donor, fill = donor)) + geom_point(shape = 21, size = 0.8, alpha = 1, stroke = 0.5) + 
  scale_color_manual("Genotype:", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
  scale_fill_manual("Genotype:", values = c("white", "black")) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("HSCT skin (genotyped through SNPs)"))
um.ge6
ta.sk.nl.sp <- ggplot(df.sk2, aes(x = tx, y = ty, color = celltype)) + geom_point(shape = 16, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col5, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("hsct skin (genotyped through SNPs)"))
ta.sk.nl.sp
ua.sk.nl.sp <- ggplot(df.sk2, aes(x = ux, y = uy, color = celltype)) + geom_point(shape = 16, size = 1.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col5, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_x_continuous(breaks = seq(-50,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
  facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("hsct skin (genotyped through SNPs)"))
ua.sk.nl.sp

# Tranfer barcodes to singler object
sr.sk$meta.data$geno.ident <- sk$geno.ident

# Number of cells per annotated cluster
sr.skn5 <- table(sr.sk$singler[[2]]$SingleR.single$labels, sr.sk$meta.data$geno.ident); knitr::kable(sr.skn5)
sr.skn6 <- table(sr.sk$singler[[2]]$SingleR.single.main$labels, sr.sk$meta.data$geno.ident); knitr::kable(sr.skn6)

# Tranfer barcodes to singler object
#sr.sk$meta.data$geno.ident <- sk@active.ident
#sr.sk$meta.data$clusters <- sk@meta.data$RNA_snn_res.1.2

# Create singler object for each genotype subset
#sr.ht <- SingleR.Subset(sr.sk, bc)
#sr.ht$meta.data$xy.um <- subset(sr.ht$meta.data$xy.um, row.names(sr.ht$meta.data$xy.um) %in% bc)
#sr.dn <- SingleR.Subset(sr.sk, colnames(sk)[!(colnames(sk) %in% bc)])
#sr.dn$meta.data$xy.um <- subset(sr.dn$meta.data$xy.um, row.names(sr.dn$meta.data$xy.um) %in% colnames(sk)[!(colnames(sk) %in% bc)])

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Prepare plots for saving
sv <- list()
sv[["ts.bd"]] <- grid.arrange(grobs = list(ts.bd + guides(color = "none"), cowplot::get_legend(ts.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ta.bd"]] <- grid.arrange(grobs = list(ta.bd + guides(color = "none"), cowplot::get_legend(ta.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ta.bd.nl"]] <- grid.arrange(grobs = list(ta.bd.nl + guides(color = "none"), cowplot::get_legend(ta.bd.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["tc.bd"]] <- grid.arrange(grobs = list(tc.bd + guides(color = "none"), cowplot::get_legend(tc.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["tc.bd.nl"]] <- grid.arrange(grobs = list(tc.bd.nl + guides(color = "none"), cowplot::get_legend(tc.bd.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["um.bd"]] <- grid.arrange(grobs = list(um.bd + guides(color = "none"), cowplot::get_legend(um.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ua.bd"]] <- grid.arrange(grobs = list(ua.bd + guides(color = "none"), cowplot::get_legend(ua.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ua.bd.nl"]] <- grid.arrange(grobs = list(ua.bd.nl + guides(color = "none"), cowplot::get_legend(ua.bd.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["uc.bd"]] <- grid.arrange(grobs = list(uc.bd + guides(color = "none"), cowplot::get_legend(uc.bd)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["uc.bd.nl"]] <- grid.arrange(grobs = list(uc.bd.nl + guides(color = "none"), cowplot::get_legend(uc.bd.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ts.sk"]] <- grid.arrange(grobs = list(ts.sk + guides(color = "none"), cowplot::get_legend(ts.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ta.sk"]] <- grid.arrange(grobs = list(ta.sk + guides(color = "none"), cowplot::get_legend(ta.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ta.sk.nl"]] <- grid.arrange(grobs = list(ta.sk.nl + guides(color = "none"), cowplot::get_legend(ta.sk.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["tc.sk"]] <- grid.arrange(grobs = list(tc.sk + guides(color = "none"), cowplot::get_legend(tc.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["tc.sk.nl"]] <- grid.arrange(grobs = list(tc.sk.nl + guides(color = "none"), cowplot::get_legend(tc.sk.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["um.sk"]] <- grid.arrange(grobs = list(um.sk + guides(color = "none"), cowplot::get_legend(um.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ua.sk"]] <- grid.arrange(grobs = list(ua.sk + guides(color = "none"), cowplot::get_legend(ua.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["ua.sk.nl"]] <- grid.arrange(grobs = list(ua.sk.nl + guides(color = "none"), cowplot::get_legend(ua.sk.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["uc.sk"]] <- grid.arrange(grobs = list(uc.sk + guides(color = "none"), cowplot::get_legend(uc.sk)), ncol = 2, 
                         top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
sv[["uc.sk.nl"]] <- grid.arrange(grobs = list(uc.sk.nl + guides(color = "none"), cowplot::get_legend(uc.sk.nl)), ncol = 2, 
                              top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))

# Save plots for blood from Seurat analysis
ggsave(plot = vp.bdc, file = paste0(dir.figs.bd, "_violin_expression_cluster.pdf"), height = 5, width = 15)
ggsave(plot = vp.bd, file = paste0(dir.figs.bd, "_violin_expression.pdf"), height = 5, width = 7)
ggsave(plot = vf.bd, file = paste0(dir.figs.bd, "_var_features.pdf"), height = 4, width = 6)
ggsave(plot = hm.pc.bd, file = paste0(dir.figs.bd, "_heatmap_pca.pdf"), height = 36, width = 20)
ggsave(plot = ep.bd, file = paste0(dir.figs.bd, "_elbowplot.pdf"), height = 4, width = 4)
ggsave(plot = hm.mk.bd, file = paste0(dir.figs.bd, "_heatmap_markers.pdf"), height = 12, width = 20)
ggsave(plot = ts.sp.bd, file = paste0(dir.figs.bd, "_tsne_split.pdf"), height = 8, width = 8)
ggsave(plot = bp.ac.bd, file = paste0(dir.figs.bd, "_boxplot_ann_confidence.pdf"), height = 8, width = 8)
ggsave(plot = bd.pvt, file = paste0(dir.figs.bd, "_blood_annotation_pval_tsne.pdf"), height = 4, width = 5)
ggsave(plot = bd.pva, file = paste0(dir.figs.bd, "_blood_annotation_pval.pdf"), height = 4, width = 5)

# Save plots for blood from SingleR analysis
ggsave(plot = sv[["ts.bd"]], file = paste0(dir.figs.bd, "_tsne.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ta.bd"]], file = paste0(dir.figs.bd, "_tsne_sc_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ta.bd.nl"]], file = paste0(dir.figs.bd, "_tsne_sc_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["tc.bd"]], file = paste0(dir.figs.bd, "_tsne_cl_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["tc.bd.nl"]], file = paste0(dir.figs.bd, "_tsne_cl_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["um.bd"]], file = paste0(dir.figs.bd, "_umap.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ua.bd"]], file = paste0(dir.figs.bd, "_umap_sc_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ua.bd.nl"]], file = paste0(dir.figs.bd, "_umap_sc_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["uc.bd"]], file = paste0(dir.figs.bd, "_umap_cl_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["uc.bd.nl"]], file = paste0(dir.figs.bd, "_umap_cl_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = ts.mk.bd, file = paste0(dir.figs.bd, "_markers_expression.pdf"), height = 8, width = 8)
tt <- ttheme_minimal(rowhead=list(fg_params=list(fontface="bold")))
pdf(paste0(dir.figs.bd, "_cell_number.pdf"), height = 18, width = 10)
grid.arrange(tableGrob(sr.bdn1, theme=tt), tableGrob(sr.bdn2, theme=tt), 
             tableGrob(sr.bdn3, theme=tt), tableGrob(sr.bdn4, theme=tt), nrow = 4, ncol = 1)
dev.off()

# Save plots for skin from Seurat analysis
ggsave(plot = vp.skc, file = paste0(dir.figs.sk, "_violin_expression_cluster.pdf"), height = 5, width = 15)
ggsave(plot = vp.sk, file = paste0(dir.figs.sk, "_violin_expression.pdf"), height = 5, width = 7)
ggsave(plot = vf.sk, file = paste0(dir.figs.sk, "_var_features.pdf"), height = 4, width = 6)
ggsave(plot = hm.pc.sk, file = paste0(dir.figs.sk, "_heatmap_pca.pdf"), height = 36, width = 20)
ggsave(plot = ep.sk, file = paste0(dir.figs.sk, "_elbowplot.pdf"), height = 4, width = 4)
ggsave(plot = hm.mk.sk, file = paste0(dir.figs.sk, "_heatmap_markers.pdf"), height = 12, width = 20)
ggsave(plot = ts.sp.sk, file = paste0(dir.figs.sk, "_tsne_split.pdf"), height = 8, width = 8)
ggsave(plot = bp.ac.sk, file = paste0(dir.figs.sk, "_boxplot_ann_confidence.pdf"), height = 8, width = 8)
ggsave(plot = sk.pvt, file = paste0(dir.figs.sk, "_skin_annotation_pval_tsne.pdf"), height = 4, width = 5)
ggsave(plot = sk.pva, file = paste0(dir.figs.sk, "_skin_annotation_pval.pdf"), height = 4, width = 5)

# Save plots for skin from SingleR analysis
ggsave(plot = sv[["ts.sk"]], file = paste0(dir.figs.sk, "_tsne.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ta.sk"]], file = paste0(dir.figs.sk, "_tsne_sc_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ta.sk.nl"]], file = paste0(dir.figs.sk, "_tsne_sc_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["tc.sk"]], file = paste0(dir.figs.sk, "_tsne_cl_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["tc.sk.nl"]], file = paste0(dir.figs.sk, "_tsne_cl_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["um.sk"]], file = paste0(dir.figs.sk, "_umap.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ua.sk"]], file = paste0(dir.figs.sk, "_umap_sc_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["ua.sk.nl"]], file = paste0(dir.figs.sk, "_umap_sc_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = sv[["uc.sk"]], file = paste0(dir.figs.sk, "_umap_cl_annotated.pdf"), height = 4, width = 8)
ggsave(plot = sv[["uc.sk.nl"]], file = paste0(dir.figs.sk, "_umap_cl_annotated_no_label.pdf"), height = 4, width = 8)
ggsave(plot = ts.mk.sk, file = paste0(dir.figs.sk, "_markers_expression.pdf"), height = 8, width = 8)
pdf(paste0(dir.figs.sk, "_cell_number.pdf"), height = 27, width = 10)
grid.arrange(tableGrob(sr.skn1, theme=tt), tableGrob(sr.skn2, theme=tt), 
             tableGrob(sr.skn3, theme=tt), tableGrob(sr.skn4, theme=tt), 
             tableGrob(sr.skn5, theme=tt), tableGrob(sr.skn6, theme=tt), nrow = 6, ncol = 1)
dev.off()

# Save plots for skin from genotyping
ggsave(plot = ts.ht, file = paste0(dir.figs.sk, "_host_tsne.pdf"), height = 3.7, width = 4.5)
ggsave(plot = um.ht, file = paste0(dir.figs.sk, "_host_umap.pdf"), height = 3.7, width = 4.5)
ggsave(plot = ts.dn, file = paste0(dir.figs.sk, "_donor_tsne.pdf"), height = 3.7, width = 4.5)
ggsave(plot = um.dn, file = paste0(dir.figs.sk, "_donor_umap.pdf"), height = 3.7, width = 4.5)
ggsave(plot = ts.ge, file = paste0(dir.figs.sk, "_genotyped_tsne2.pdf"), height = 3.7, width = 4.5)
ggsave(plot = um.ge, file = paste0(dir.figs.sk, "_genotyped_umap2.pdf"), height = 3.7, width = 4.5)
ggsave(plot = ts.ge2, file = paste0(dir.figs.sk, "_genotyped_tsne.pdf"), height = 4, width = 4.5)
ggsave(plot = ts.ge3, file = paste0(dir.figs.sk, "_genotyped_tsne_bw.pdf"), height = 4, width = 4.5)
ggsave(plot = ts.ge4, file = paste0(dir.figs.sk, "_genotyped_tsne_bw2.pdf"), height = 4, width = 4.5)
ggsave(plot = ts.ge5, file = paste0(dir.figs.sk, "_genotyped_tsne_bw3.pdf"), height = 4, width = 4.5)
ggsave(plot = ts.ge6, file = paste0(dir.figs.sk, "_genotyped_tsne_bw4.pdf"), height = 4, width = 4.5)
ggsave(plot = um.ge2, file = paste0(dir.figs.sk, "_genotyped_umap.pdf"), height = 4, width = 4.5)
ggsave(plot = um.ge3, file = paste0(dir.figs.sk, "_genotyped_umap_bw.pdf"), height = 4, width = 4.5)
ggsave(plot = um.ge4, file = paste0(dir.figs.sk, "_genotyped_umap_bw2.pdf"), height = 4, width = 4.5)
ggsave(plot = um.ge5, file = paste0(dir.figs.sk, "_genotyped_umap_bw3.pdf"), height = 4, width = 4.5)
ggsave(plot = um.ge6, file = paste0(dir.figs.sk, "_genotyped_umap_bw4.pdf"), height = 4, width = 4.5)
ggsave(plot = ta.sk.nl.sp, file = paste0(dir.figs.sk, "_tsne_sc_annotated_no_label_split.pdf"), height = 4, width = 9)
ggsave(plot = ua.sk.nl.sp, file = paste0(dir.figs.sk, "_umap_sc_annotated_no_label_split.pdf"), height = 4, width = 9)

#-----------------------------------------------------------------------------------------------------------------#
# End of script
#-----------------------------------------------------------------------------------------------------------------#

#=================================================================================================================#
# R session info
#=================================================================================================================#

#> devtools::session_info()
#─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────
#setting  value                       
#version  R version 3.5.1 (2018-07-02)
#os       macOS  10.14.2              
#system   x86_64, darwin15.6.0        
#ui       RStudio                     
#language (EN)                        
#collate  en_US.UTF-8                 
#ctype    en_US.UTF-8                 
#tz       Europe/Berlin               
#date     2018-12-15                  
#
#─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────
#package       * version    date       lib source                           
#AnnotationDbi   1.44.0     2018-10-30 [1] Bioconductor                     
#assertthat      0.2.0      2017-04-11 [1] CRAN (R 3.5.0)                   
#backports       1.1.2      2017-12-13 [1] CRAN (R 3.5.0)                   
#bibtex          0.4.2      2017-06-30 [1] CRAN (R 3.5.0)                   
#bindr           0.1.1      2018-03-13 [1] CRAN (R 3.5.0)                   
#bindrcpp        0.2.2      2018-03-29 [1] CRAN (R 3.5.0)                   
#Biobase         2.42.0     2018-10-30 [1] Bioconductor                     
#BiocGenerics    0.28.0     2018-10-30 [1] Bioconductor                     
#BiocParallel    1.16.2     2018-11-28 [1] Bioconductor                     
#bit             1.1-14     2018-05-29 [1] CRAN (R 3.5.0)                   
#bit64           0.9-7      2017-05-08 [1] CRAN (R 3.5.0)                   
#bitops          1.0-6      2013-08-17 [1] CRAN (R 3.5.0)                   
#blob            1.1.1      2018-03-25 [1] CRAN (R 3.5.0)                   
#callr           3.1.0      2018-12-10 [1] CRAN (R 3.5.0)                   
#caTools         1.17.1.1   2018-07-20 [1] CRAN (R 3.5.0)                   
#cli             1.0.1      2018-09-25 [1] CRAN (R 3.5.0)                   
#cluster         2.0.7-1    2018-04-13 [1] CRAN (R 3.5.1)                   
#codetools       0.2-15     2016-10-05 [1] CRAN (R 3.5.1)                   
#colorspace      1.3-2      2016-12-14 [1] CRAN (R 3.5.0)                   
#cowplot         0.9.3      2018-07-15 [1] CRAN (R 3.5.0)                   
#crayon          1.3.4      2017-09-16 [1] CRAN (R 3.5.0)                   
#data.table      1.11.8     2018-09-30 [1] CRAN (R 3.5.0)                   
#DBI             1.0.0      2018-05-02 [1] CRAN (R 3.5.0)                   
#desc            1.2.0      2018-05-01 [1] CRAN (R 3.5.0)                   
#devtools        2.0.1      2018-10-26 [1] CRAN (R 3.5.1)                   
#digest          0.6.18     2018-10-10 [1] CRAN (R 3.5.0)                   
#DO.db           2.9        2018-11-02 [1] Bioconductor                     
#DOSE            3.8.0      2018-10-30 [1] Bioconductor                     
#dplyr         * 0.7.8      2018-11-10 [1] CRAN (R 3.5.0)                   
#fastmatch       1.1-0      2017-01-28 [1] CRAN (R 3.5.0)                   
#fgsea           1.8.0      2018-10-30 [1] Bioconductor                     
#fitdistrplus    1.0-11     2018-09-10 [1] CRAN (R 3.5.0)                   
#fs              1.2.6      2018-08-23 [1] CRAN (R 3.5.0)                   
#future          1.10.0     2018-10-17 [1] CRAN (R 3.5.0)                   
#future.apply    1.0.1      2018-08-26 [1] CRAN (R 3.5.0)                   
#gbRd            0.4-11     2012-10-01 [1] CRAN (R 3.5.0)                   
#gdata           2.18.0     2017-06-06 [1] CRAN (R 3.5.0)                   
#ggplot2         3.1.0      2018-10-25 [1] CRAN (R 3.5.0)                   
#ggrepel         0.8.0      2018-05-09 [1] CRAN (R 3.5.0)                   
#ggridges        0.5.1      2018-09-27 [1] CRAN (R 3.5.0)                   
#globals         0.12.4     2018-10-11 [1] CRAN (R 3.5.0)                   
#glue            1.3.0      2018-07-17 [1] CRAN (R 3.5.0)                   
#GO.db           3.7.0      2018-11-02 [1] Bioconductor                     
#GOSemSim        2.8.0      2018-10-30 [1] Bioconductor                     
#gplots          3.0.1      2016-03-30 [1] CRAN (R 3.5.0)                   
#gridExtra       2.3        2017-09-09 [1] CRAN (R 3.5.0)                   
#gtable          0.2.0      2016-02-26 [1] CRAN (R 3.5.0)                   
#gtools          3.8.1      2018-06-26 [1] CRAN (R 3.5.0)                   
#htmltools       0.3.6      2017-04-28 [1] CRAN (R 3.5.0)                   
#htmlwidgets     1.3        2018-09-30 [1] CRAN (R 3.5.0)                   
#httr            1.4.0      2018-12-11 [1] CRAN (R 3.5.1)                   
#ica             1.0-2      2018-05-24 [1] CRAN (R 3.5.0)                   
#igraph          1.2.2      2018-07-27 [1] CRAN (R 3.5.0)                   
#IRanges         2.16.0     2018-10-30 [1] Bioconductor                     
#irlba           2.3.2      2018-01-11 [1] CRAN (R 3.5.0)                   
#jsonlite        1.6        2018-12-07 [1] CRAN (R 3.5.0)                   
#KernSmooth      2.23-15    2015-06-29 [1] CRAN (R 3.5.1)                   
#lattice         0.20-38    2018-11-04 [1] CRAN (R 3.5.0)                   
#lazyeval        0.2.1      2017-10-29 [1] CRAN (R 3.5.0)                   
#listenv         0.7.0      2018-01-21 [1] CRAN (R 3.5.0)                   
#lmtest          0.9-36     2018-04-04 [1] CRAN (R 3.5.0)                   
#lsei            1.2-0      2017-10-23 [1] CRAN (R 3.5.0)                   
#magrittr        1.5        2014-11-22 [1] CRAN (R 3.5.0)                   
#MASS            7.3-51.1   2018-11-01 [1] CRAN (R 3.5.1)                   
#Matrix          1.2-15     2018-11-01 [1] CRAN (R 3.5.1)                   
#memoise         1.1.0      2017-04-21 [1] CRAN (R 3.5.0)                   
#metap           1.0        2018-07-25 [1] CRAN (R 3.5.0)                   
#munsell         0.5.0      2018-06-12 [1] CRAN (R 3.5.0)                   
#npsurv          0.4-0      2017-10-14 [1] CRAN (R 3.5.0)                   
#openxlsx      * 4.1.0      2018-05-26 [1] CRAN (R 3.5.0)                   
#pbapply         1.3-4      2018-01-10 [1] CRAN (R 3.5.0)                   
#pillar          1.3.0      2018-07-14 [1] CRAN (R 3.5.0)                   
#pkgbuild        1.0.2      2018-10-16 [1] CRAN (R 3.5.0)                   
#pkgconfig       2.0.2      2018-08-16 [1] CRAN (R 3.5.0)                   
#pkgload         1.0.2      2018-10-29 [1] CRAN (R 3.5.0)                   
#plotly          4.8.0      2018-07-20 [1] CRAN (R 3.5.1)                   
#plyr            1.8.4      2016-06-08 [1] CRAN (R 3.5.0)                   
#png             0.1-7      2013-12-03 [1] CRAN (R 3.5.0)                   
#prettyunits     1.0.2      2015-07-13 [1] CRAN (R 3.5.0)                   
#processx        3.2.1      2018-12-05 [1] CRAN (R 3.5.1)                   
#ps              1.2.1      2018-11-06 [1] CRAN (R 3.5.1)                   
#purrr           0.2.5      2018-05-29 [1] CRAN (R 3.5.0)                   
#qvalue          2.14.0     2018-10-30 [1] Bioconductor                     
#R.methodsS3     1.7.1      2016-02-16 [1] CRAN (R 3.5.0)                   
#R.oo            1.22.0     2018-04-22 [1] CRAN (R 3.5.0)                   
#R.utils         2.7.0      2018-08-27 [1] CRAN (R 3.5.0)                   
#R6              2.3.0      2018-10-04 [1] CRAN (R 3.5.0)                   
#RANN            2.6        2018-07-16 [1] CRAN (R 3.5.0)                   
#RColorBrewer    1.1-2      2014-12-07 [1] CRAN (R 3.5.0)                   
#Rcpp            1.0.0      2018-11-07 [1] CRAN (R 3.5.0)                   
#Rdpack          0.10-1     2018-10-04 [1] CRAN (R 3.5.0)                   
#remotes         2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
#reshape2        1.4.3      2017-12-11 [1] CRAN (R 3.5.0)                   
#reticulate      1.10       2018-08-05 [1] CRAN (R 3.5.0)                   
#rlang           0.3.0.1    2018-10-25 [1] CRAN (R 3.5.0)                   
#ROCR            1.0-7      2015-03-26 [1] CRAN (R 3.5.0)                   
#rprojroot       1.3-2      2018-01-03 [1] CRAN (R 3.5.0)                   
#RSQLite         2.1.1      2018-05-06 [1] CRAN (R 3.5.0)                   
#rstudioapi      0.8        2018-10-02 [1] CRAN (R 3.5.0)                   
#rsvd            1.0.0      2018-11-06 [1] CRAN (R 3.5.0)                   
#Rtsne           0.15       2018-11-10 [1] CRAN (R 3.5.0)                   
#S4Vectors       0.20.1     2018-11-09 [1] Bioconductor                     
#scales          1.0.0      2018-08-09 [1] CRAN (R 3.5.0)                   
#SDMTools        1.1-221    2014-08-05 [1] CRAN (R 3.5.0)                   
#sessioninfo     1.1.1      2018-11-05 [1] CRAN (R 3.5.1)                   
#Seurat        * 3.0.0.9000 2018-12-12 [1] Github (satijalab/seurat@aa9ffda)
#stringi         1.2.4      2018-07-20 [1] CRAN (R 3.5.0)                   
#stringr         1.3.1      2018-05-10 [1] CRAN (R 3.5.0)                   
#survival        2.43-3     2018-11-26 [1] CRAN (R 3.5.0)                   
#tibble          1.4.2      2018-01-22 [1] CRAN (R 3.5.0)                   
#tidyr           0.8.2      2018-10-28 [1] CRAN (R 3.5.0)                   
#tidyselect      0.2.5      2018-10-11 [1] CRAN (R 3.5.0)                   
#tsne            0.1-3      2016-07-15 [1] CRAN (R 3.5.0)                   
#usethis         1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
#viridisLite     0.3.0      2018-02-01 [1] CRAN (R 3.5.0)                   
#withr           2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
#yaml            2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
#zip             1.0.0      2017-04-25 [1] CRAN (R 3.5.0)                   
#zoo             1.8-4      2018-09-19 [1] CRAN (R 3.5.0)                   