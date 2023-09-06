#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a transplanted
# AML patient (chimeric sample)
# Date: 2019.04.05
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://satijalab.org/seurat/pancreas_integration_label_transfer.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger
# - transpose skin data onto blood data healthy donor (with seurat CCA) 
# using 1000 genes instead of the default 2000 (to increase homogeneity of the projection)
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
experimentid <- c("EX0004")

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub(".R", "", stringr::str_split(rstudioapi::getSourceEditorContext()$path, "/", simplify = T)[,5])

# Define samples analyzed (all)
#sa <- c("cd")

# Define reference sample name
ct <- c("blood")
tr <- c("skin")

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c("sc_hd_10x_skin_vs_blood")
an.descs <- c("sc_hd_skin_vs_blood")

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
dir.create(paste0("figures/", analysisid, "_", an.desc, "/seurat"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
dir.figs <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
dir.figss <- paste0("figures/", analysisid, "_", an.desc, "/seurat/", analysisid, "_", an.descs)
#dir.figs.bd <- paste0("figures/", analysisid, "_", an.desc, "/blood/", analysisid, "_", an.descs, "_blood")
#dir.figs.sk <- paste0("figures/", analysisid, "_", an.desc, "/skin/", analysisid, "_", an.descs, "_skin")

#=================================================================================================================#
# Analysis of processed data for blood sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Data set integration with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Load blood Seurat object from healthy donor and transplanted patient
hd <- list()
hd[["bd"]] <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_seurat_blood.rds"))
hd$bd@meta.data[6:19] <- NULL
hd$bd@meta.data[["donor"]] <- c("blood (healthy)")
hd[["sk"]] <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_seurat_skin.rds"))
#hd$sk@meta.data[c(5,7)] <- NULL
hd$sk@meta.data[["donor"]] <- c("skin (healthy)")

# Load blood singleR object from healthy donor and transplanted patient
sr <- list()
sr[["bd"]] <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_singler_blood.rds"))
sr[["sk"]] <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_singler_skin.rds"))

# Transfer singleR cell type annotation to Seurat objects
hd$bd@meta.data[["celltype"]] <- sr$bd[["singler"]][[2]][["SingleR.single"]][["labels"]]
hd$sk@meta.data[["celltype"]] <- sr$sk[["singler"]][[2]][["SingleR.single"]][["labels"]]
#saveRDS(hd, file = paste0(dir.figs, "_merged_seurat.rds"))
#hd <- readRDS(paste0(dir.figs, "_merged_seurat.rds"))

# Save original clustering
hd$bd@meta.data$RNA_snn_res.bd <- hd$bd@meta.data$RNA_snn_res.1
hd$sk@meta.data$RNA_snn_res.sk <- hd$sk@meta.data$RNA_snn_res.1
hd$bd@meta.data[["RNA_snn_res.1"]] <- NULL
hd$sk@meta.data[["RNA_snn_res.1"]] <- NULL

# Identify anchors between the individual datasets
hd.anc <- FindIntegrationAnchors(object.list = hd, dims = 1:20, anchor.features = 1000, scale = F)

# Integrate data
hd.int <- IntegrateData(anchorset = hd.anc, dims = 1:20, features.to.integrate = rownames(hd$bd))

# Switch to integrated assay (the variable features of this assay are automatically set during IntegrateData)
DefaultAssay(hd.int) <- "integrated"

# Run the standard workflow for visualization
#hd.int <- ScaleData(hd.int, features = rownames(hd.int))
hd.int <- ScaleData(hd.int, features = rownames(hd.int), vars.to.regress = c("nCount_RNA", "percent.mito"))
hd.int <- RunPCA(hd.int, verbose = F)
ep.pc <- ElbowPlot(hd.int, ndims = 35)
hd.int <- RunTSNE(hd.int, reduction = "pca", dims = 1:20, seed.use = 10, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30)
hd.int <- RunUMAP(hd.int, reduction = "pca", dims = 1:20, seed.use = 1, spread = 10, min.dist = 0.001, reduction.name = "umap")
hd.int@meta.data[["orig.ident_old"]] <- hd.int@meta.data[["orig.ident"]]
hd.int@meta.data[["orig.ident"]] <- c("healthy")
#saveRDS(hd.int, file = paste0(dir.figs, "_integrated_seurat.rds"))
hd.int <- readRDS(paste0(dir.figs, "_integrated_seurat.rds"))

#-----------------------------------------------------------------------------------------------------------------#
# Clustering and differential expression analysis with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Find clusters
hd.int <- FindNeighbors(hd.int, reduction = "pca", k.param = 30, dims = 1:20, do.plot = F)
hd.int <- FindClusters(hd.int, resolution = 0.7, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)

# Define colors
col1 <- c("darkorange2", "firebrick1", "darkgoldenrod1", "steelblue1", "royalblue1", "plum3", "grey40", 
          "grey70", "black", "black", "seagreen3", "orchid1", "wheat2", "darkorchid2", "black", "olivedrab")
col2 <- c("darkorange2", "firebrick1", "darkgoldenrod1", "steelblue1", "royalblue1", "plum3", 
          "grey70", "seagreen3", "wheat2", "orchid1", "darkorchid2", "black", "olivedrab")

# Draw plots
sts.hd <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "integrated_snn_res.0.7", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("tSNE 1") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + ylab("tSNE 2") #+ NoLegend()
sts.hd.ss <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "integrated_snn_res.0.7", split.by = "donor", label = T, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("tSNE 1") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + ylab("tSNE 2") #+ NoLegend()
sta.hd <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "celltype", split.by = "donor", label = T, label.size = 4, pt.size = 0.5, cols = col1) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sta.hd.nl <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "celltype", split.by = "donor", label = F, label.size = 4, pt.size = 0.5, cols = col1) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25))#+ NoLegend()
sts.hd.sc <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "donor", split.by = "integrated_snn_res.0.7", label = F, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + scale_x_continuous(breaks = seq(-25,25,25))#+ NoLegend()
sts.hd.sa <- DimPlot(hd.int, reduction = "tsne", dims = c(1,2), group.by = "donor", split.by = "celltype", label = F, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()

sum.hd <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "integrated_snn_res.0.7", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sum.hd.ss <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "integrated_snn_res.0.7", split.by = "donor", label = T, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sua.hd <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "celltype", split.by = "donor", label = T, label.size = 4, pt.size = 0.5, cols = col1) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sua.hd.nl <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "celltype", split.by = "donor", label = F, label.size = 4, pt.size = 0.5, cols = col1) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sum.hd.sc <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "donor", split.by = "integrated_snn_res.0.7", label = F, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
sum.hd.sa <- DimPlot(hd.int, reduction = "umap", dims = c(1,2), group.by = "donor", split.by = "celltype", label = F, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()

# Number of overlapping cells per annotated cluster
ov.hd <- table(hd.int@meta.data$celltype, hd.int@meta.data$donor); knitr::kable(ov.hd)

#-----------------------------------------------------------------------------------------------------------------#
# Data vizualization
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data to visualize tSNE embeddings
df <- data.frame(tx = hd.int@reductions$tsne@cell.embeddings[,1],
                    ty = hd.int@reductions$tsne@cell.embeddings[,2],
                    ux = hd.int@reductions$umap@cell.embeddings[,1],
                    uy = hd.int@reductions$umap@cell.embeddings[,2],
                    cell = colnames(hd.int), 
                    donor = hd.int@meta.data$donor, 
                    celltype = hd.int@meta.data$celltype, 
                    cluster = hd.int@meta.data$integrated_snn_res.0.7)
#df <- reshape2::melt(df, id.vars = c('tx','ty','ux','uy'))
df <- reshape2::melt(df, id.vars = c('tx','ty','ux','uy','donor'))

# Plot tSNE
ts.dn <- ggplot(df[df$variable == "cell",], aes(x = tx, y = ty, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), 
                     guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.dn.sp <- ggplot(df[df$variable == "cell",], aes(x = tx, y = ty, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), 
                     guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.ct <- ggplot(df[df$variable == "celltype",], aes(x = tx, y = ty, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.ct.sp <- ggplot(df[df$variable == "celltype",], aes(x = tx, y = ty, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~value, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.ct.dn <- ggplot(df[df$variable == "celltype",], aes(x = tx, y = ty, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.ct.sp.dn <- ggplot(df[df$variable == "celltype",], aes(x = tx, y = ty, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~value, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.cl <- ggplot(df[df$variable == "cluster",], aes(x = tx, y = ty, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.cl.sp <- ggplot(df[df$variable == "cluster",], aes(x = tx, y = ty, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~factor(value, levels = 0:max(as.numeric(value))), ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.cl.dn <- ggplot(df[df$variable == "cluster",], aes(x = tx, y = ty, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.cl.sp.dn <- ggplot(df[df$variable == "cluster",], aes(x = tx, y = ty, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~factor(value, levels = 0:max(as.numeric(value))), ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")
ts.ct.dn.cn <- ggplot(df[df$variable == "celltype",], aes(x = tx, y = ty)) + 
  geom_point(data = df[df$variable == "celltype" & df$donor == "skin (healthy)",], aes(color = value), shape = 16, size = 1.2, alpha = 0.5) + 
  geom_density_2d(data = df[df$variable == "celltype" & df$donor == "blood (healthy)",], color = "grey30", size = 0.2) + 
  scale_color_manual("Skin:", values = col2, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25), limits = c(-42,35)) + scale_y_continuous(breaks = c(-25,0,25), limits = c(-45,45)) + #facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle("Healthy donor")

# Plot UMAP
um.dn <- ggplot(df[df$variable == "cell",], aes(x = ux, y = uy, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), 
                     guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.dn.sp <- ggplot(df[df$variable == "cell",], aes(x = ux, y = uy, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Sample:", values = c("black", "firebrick2"), 
                     guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.ct <- ggplot(df[df$variable == "celltype",], aes(x = ux, y = uy, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.ct.sp <- ggplot(df[df$variable == "celltype",], aes(x = ux, y = uy, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~value, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.ct.dn <- ggplot(df[df$variable == "celltype",], aes(x = ux, y = uy, color = value)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.ct.sp.dn <- ggplot(df[df$variable == "celltype",], aes(x = ux, y = uy, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cell type:", values = c("black", "firebrick2"), guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~value, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.cl <- ggplot(df[df$variable == "cluster",], aes(x = ux, y = uy, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.cl.sp <- ggplot(df[df$variable == "cluster",], aes(x = ux, y = uy, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~factor(value, levels = 0:max(as.numeric(value))), ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.cl.dn <- ggplot(df[df$variable == "cluster",], aes(x = ux, y = uy, color = factor(value, levels = 0:max(as.numeric(value))))) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = col1, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.cl.sp.dn <- ggplot(df[df$variable == "cluster",], aes(x = ux, y = uy, color = donor)) + 
  geom_point(shape = 16, size = 0.5, alpha = 0.5) + 
  scale_color_manual("Cluster:", values = c("black", "firebrick2"), guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + facet_wrap(~factor(value, levels = 0:max(as.numeric(value))), ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")
um.ct.dn.cn <- ggplot(df[df$variable == "celltype",], aes(x = ux, y = uy)) + 
  geom_point(data = df[df$variable == "celltype" & df$donor == "skin (healthy)",], aes(color = value), shape = 16, size = 1.2, alpha = 0.5) + 
  geom_density_2d(data = df[df$variable == "celltype" & df$donor == "blood (healthy)",], color = "grey30", size = 0.2) + 
  scale_color_manual("Skin:", values = col2, guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) + 
  scale_x_continuous(breaks = c(-25,0,25), limits = c(-51,33)) + scale_y_continuous(breaks = c(-25,0,25), limits = c(-50,44)) + #facet_wrap(~donor, ncol = 4) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("Healthy donor")

#-----------------------------------------------------------------------------------------------------------------#
# Projecting data set onto a reference
#-----------------------------------------------------------------------------------------------------------------#

## Identify transfer anchors between the individual datasets
#hd.trs <- FindTransferAnchors(reference = hd$bd, query = hd$sk, dims = 1:20)
#
## Add predicted IDs to query data set
#hd.prd <- TransferData(anchorset = hd.trs, refdata = hd$bd$celltype, dims = 1:20)
#hd.sk <- AddMetaData(object = hd$sk, metadata = hd.prd)
#
## Evaluate predicted cell type annotations
#hd$sk$prediction.match <- hd$sk$predicted.id == hd$sk$celltype
#table(hd$sk$predicted.id, hd$sk$celltype)
#table(hd$sk$prediction.match)
#table(hd$sk$predicted.id)
#table(hd$sk$celltype)

#-----------------------------------------------------------------------------------------------------------------#
# Projecting data onto an existing umap embedding
#-----------------------------------------------------------------------------------------------------------------#

## Store original umap embeddings
#hd$bd@reductions$umap_orig <- hd$bd@reductions$umap
#hd$sk@reductions$umap_orig <- hd$sk@reductions$umap
#
## Extract pca coordinates from healthy donor and transplanted patient
#hd.bd.pca <- Embeddings(hd$bd, reduction = "pca")[, paste0("PC_", 1:20), drop = F]
#hd.sk.pca <- Embeddings(hd$sk, reduction = "pca")[, paste0("PC_", 1:20), drop = F]
#
## Run umap on the PCA values (calling umap-learn directly)
#bd.um <- umap::umap(hd.bd.pca, method = "umap-learn", n_neighbors = 30L, 
#                    min_dist = 0.001, random_state = 1L, 
#                    spread = 10, n_components = 2L, metric = "correlation", #n_epochs = 200, 
#                    set.op.mix.ratio = 1, local_connectivity = 1, negative.sample.rate = 5L, 
#                    a = umap:::find.ab.params(10,0.001)["a"], 
#                    b = umap:::find.ab.params(10,0.001)["b"], verbose = T)
#hd$bd@reductions$umap@cell.embeddings[,1:2] <- bd.um$layout
#
## Plot umap
#DimPlot(hd$bd, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.bd", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
#
## Project transplanted patient data set
##tp.um <- umap:::umap.learn.predict(hd.um, bd.tp.pca)
#sk.um <- bd.um$UMAP$transform(hd.sk.pca)
#hd$sk@reductions$umap@cell.embeddings[,1:2] <- sk.um
#
## Plot umap
#DimPlot(hd$sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.sk", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
#
## Store umap embeddings
#bd$hd@reductions$umap_pc20 <- bd$hd@reductions$umap
#bd$tp@reductions$umap_pc20 <- bd$tp@reductions$umap
#
## Extract scaled count data from healthy donor and transplanted patient filtered for merged top variable features
#var.feat <- union(rownames(HVFInfo(bd$hd)[rownames(HVFInfo(bd$hd)) %in% VariableFeatures(bd$hd)[1:8000],]), 
#                  rownames(HVFInfo(bd$tp)[rownames(HVFInfo(bd$tp)) %in% VariableFeatures(bd$tp)[1:8000],]))
#bd.hd.um <- t(as.matrix(GetAssayData(bd$hd, slot = "scale.data", assay = "RNA")[var.feat,]))
#bd.tp.um <- t(as.matrix(GetAssayData(bd$tp, slot = "scale.data", assay = "RNA")[var.feat,]))
#
## Use umap-learn
#hd.um <- umap::umap(bd.hd.um, method = "umap-learn", n_neighbors = 30L, 
#                    min_dist = 0.001, random_state = 1L, 
#                    spread = 10, n_components = 2L, metric = "correlation", #n_epochs = 200, 
#                    set.op.mix.ratio = 1, local_connectivity = 1, negative.sample.rate = 5L, 
#                    a = umap:::find.ab.params(10,0.001)["a"], 
#                    b = umap:::find.ab.params(10,0.001)["b"], verbose = T)
#bd$hd@reductions$umap@cell.embeddings[,1:2] <- hd.um$layout
#
## Plot umap
#DimPlot(bd$hd, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
#
## Project transplanted patient data set
#tp.um <- hd.um$UMAP$transform(bd.tp.um)
#bd$tp@reductions$umap@cell.embeddings[,1:2] <- tp.um
#
## Plot umap
#DimPlot(bd$tp, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1.51", split.by = "orig.ident", label = T, label.size = 4, pt.size = 0.8, cols = NULL) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
#
## Store original umap embeddings
##saveRDS(bd, file = paste0(dir.figs, "_merged_seurat_umap8000.rds"))
#bd$hd@reductions$umap_8000 <- bd$hd@reductions$umap
#bd$tp@reductions$umap_8000 <- bd$tp@reductions$umap
##bd$hd@reductions$umap <- bd$hd@reductions$umap_8000
##bd$tp@reductions$umap <- bd$tp@reductions$umap_8000

#=================================================================================================================#
# Cluster annotation
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Biological annotation of clusters from blood sample using SingleR
#-----------------------------------------------------------------------------------------------------------------#

# Create the SingleR object using clusters calculated with Seurat
hd.pd <- GetAssayData(hd.int, slot = "data", assay = "integrated")
hd.pd@Dimnames[[2]] <- paste0(1:ncol(hd.pd), "_", hd.pd@Dimnames[[2]])
sr.hd <- CreateSinglerObject(hd.pd, 
                             annot = NULL, project.name = "integrated_blood_ann", min.genes = 200, #500
                             technology = "10X", species = "Human", citation = "",
                             ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                             fine.tune = T, do.signatures = F, clusters = NULL, do.main.types = T, 
                             reduce.file.size = T, numCores = 4)#SingleR.numCores)
saveRDS(sr.hd, file = paste0(dir.figs, "_singler.rds"))
#sr.hd <- readRDS(paste0(dir.figs.bd, "_singler.rds"))

# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
#sr.hd$seurat = seurat.object # (optional)
sr.hd$meta.data$orig.ident <- hd.int@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.hd$meta.data$donor <- hd.int@meta.data$donor
sr.hd$meta.data$celltype <- hd.int@meta.data$celltype
sr.hd$meta.data$RNA_snn_res.bd <- hd.int@meta.data$RNA_snn_res.bd
sr.hd$meta.data$RNA_snn_res.sk <- hd.int@meta.data$RNA_snn_res.sk
sr.hd$meta.data$xy <- hd.int@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = "tsne"))
sr.hd$meta.data$xy.um <- hd.int@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = "umap"))
sr.hd$meta.data$clusters <- hd.int@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(bd))

# Number of cells per annotated cluster
sr.hdn1 <- table(sr.hd$meta.data$orig.ident, sr.hd$meta.data$clusters); knitr::kable(sr.hdn1)
sr.hdn2 <- table(sr.hd$singler[[2]]$SingleR.single$labels, sr.hd$meta.data$orig.ident); knitr::kable(sr.hdn2)
sr.hdn3 <- table(sr.hd$singler[[2]]$SingleR.single$labels, sr.hd$meta.data$clusters); knitr::kable(sr.hdn3)
sr.hdn4 <- table(sr.hd$singler[[2]]$SingleR.single.main$labels, sr.hd$meta.data$orig.ident); knitr::kable(sr.hdn4)
sr.hdn5 <- table(sr.hd$meta.data$celltype, sr.hd$meta.data$clusters); knitr::kable(sr.hdn5)
sr.hdn6 <- table(sr.hd$meta.data$celltype, sr.hd$meta.data$donor); knitr::kable(sr.hdn6)

# Heatmap of the aggregated scores before fine-tuning for the main cell types:
#SingleR.DrawHeatmap(sr.hd$singler[[2]]$SingleR.single, top.n = Inf, clusters = sr.hd$meta.data$cluster)#orig.ident)

# Draw tSNE plots swith cluster from Seurat object
ts.bd <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy, do.labels = T, do.letters = F, col = col1, 
                          labels = sr.hd$meta.data$clusters, label.size = 4, dot.size = 1, title = "blood")$p + 
  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
# Draw tSNE plots with annotation from SingleR at single cell level from individual objects
ta.bd <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy, do.labels = T, do.letters = F, col = col2, 
                          labels = sr.hd$meta.data$celltype, label.size = 4, dot.size = 1, title = "blood")$p + 
  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
ta.bd.nl <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy, do.labels = F, do.letters = F, col = col2, 
                             labels = sr.hd$meta.data$celltype, label.size = 4, dot.size = 1, title = "blood")$p + 
  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
# Draw UMAP plots swith cluster from Seurat object
um.bd <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy.um, do.labels = T, do.letters = F, col = col1, 
                          labels = sr.hd$meta.data$clusters, label.size = 4, dot.size = 1, title = "blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
# Draw UMAP plots with annotation from SingleR at single cell level
ua.bd <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy.um, do.labels = T, do.letters = F, col = col2, 
                          labels = sr.hd$meta.data$celltype, label.size = 4, dot.size = 1, title = "blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
ua.bd.nl <- SingleR.PlotTsne(sr.hd$singler[[2]]$SingleR.single, sr.hd$meta.data$xy.um, do.labels = F, do.letters = F, col = col2, 
                             labels = sr.hd$meta.data$celltype, label.size = 4, dot.size = 1, title = "blood")$p + 
  xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Funciton to save plots with same size
grid <- function(data, layout=rbind(c(1,2))) {
  grid.arrange(grobs = list(data + guides(color = "none"), cowplot::get_legend(data)), ncol = 2, layout_matrix = layout, 
               top = grid::textGrob("", gp = grid::gpar(fontsize = 16, font = 2)))
}

# Save plots for blood from Seurat analysis
ggsave(plot = sts.hd.sc, file = paste0(dir.figss, "_tsne_split_cluster.pdf"), height = 8, width = 8)
ggsave(plot = sts.hd.sa, file = paste0(dir.figss, "_tsne_split_ann.pdf"), height = 8, width = 8)
ggsave(plot = sts.hd.ss, file = paste0(dir.figss, "_tsne_split_sample.pdf"), height = 4, width = 8)
ggsave(plot = sum.hd.sc, file = paste0(dir.figss, "_umap_split_cluster.pdf"), height = 8, width = 8)
ggsave(plot = sum.hd.sa, file = paste0(dir.figss, "_umap_split_ann.pdf"), height = 8, width = 8)
ggsave(plot = sum.hd.ss, file = paste0(dir.figss, "_umap_split_sample.pdf"), height = 4, width = 8)

ggsave(plot = grid(sts.hd), file = paste0(dir.figss, "_tsne.pdf"), height = 4, width = 8)
ggsave(plot = grid(sta.hd), file = paste0(dir.figss, "_tsne_sc_annotated.pdf"), height = 4, width = 12)
ggsave(plot = grid(sta.hd.nl), file = paste0(dir.figss, "_tsne_sc_annotated_no_label.pdf"), height = 4, width = 12)
ggsave(plot = grid(sum.hd), file = paste0(dir.figss, "_umap.pdf"), height = 4, width = 8)
ggsave(plot = grid(sua.hd), file = paste0(dir.figss, "_umap_sc_annotated.pdf"), height = 4, width = 12)
ggsave(plot = grid(sua.hd.nl), file = paste0(dir.figss, "_umap_sc_annotated_no_label.pdf"), height = 4, width = 12)

# Save plots for blood from SingleR analysis
ggsave(plot = grid(ts.dn, rbind(c(1,2))), file = paste0(dir.figs, "_tsne.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.dn.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_split.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.ct, rbind(c(1,2))), file = paste0(dir.figs, "_tsne_celltype.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.ct.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_celltype_split.pdf"), height = 10, width = 12)
ggsave(plot = grid(ts.ct.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_celltype_donor.pdf"), height = 4, width = 9)
ggsave(plot = grid(ts.ct.dn.cn, rbind(c(1,2))), file = paste0(dir.figs, "_tsne_celltype_donor_contour.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.ct.sp.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_celltype_split_donor.pdf"), height = 10, width = 12)
ggsave(plot = grid(ts.cl, rbind(c(1,2))), file = paste0(dir.figs, "_tsne_cluster.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.cl.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_cluster_split.pdf"), height = 10, width = 12)
ggsave(plot = grid(ts.cl.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_cluster_donor.pdf"), height = 4, width = 8)
ggsave(plot = grid(ts.cl.sp.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_tsne_cluster_split_donor.pdf"), height = 10, width = 12)

ggsave(plot = grid(um.dn, rbind(c(1,2))), file = paste0(dir.figs, "_umap.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.dn.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_split.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.ct, rbind(c(1,2))), file = paste0(dir.figs, "_umap_celltype.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.ct.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_celltype_split.pdf"), height = 10, width = 12)
ggsave(plot = grid(um.ct.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_celltype_donor.pdf"), height = 4, width = 9)
ggsave(plot = grid(um.ct.dn.cn, rbind(c(1,2))), file = paste0(dir.figs, "_umap_celltype_donor_contour.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.ct.sp.dn, rbind(c(1,2))), file = paste0(dir.figs, "_umap_celltype_split_donor.pdf"), height = 10, width = 12)
ggsave(plot = grid(um.cl, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_cluster.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.cl.sp, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_cluster_split.pdf"), height = 10, width = 12)
ggsave(plot = grid(um.cl.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_cluster_donor.pdf"), height = 4, width = 8)
ggsave(plot = grid(um.cl.sp.dn, rbind(c(1,1,2))), file = paste0(dir.figs, "_umap_cluster_split_donor.pdf"), height = 10, width = 12)

# Save table with cell type overlap
tt <- ttheme_minimal(rowhead=list(fg_params=list(fontface="bold")))
pdf(paste0(dir.figs, "_celltype_overlap.pdf"), height = 8, width = 8)
grid.arrange(tableGrob(ov.hd, theme=tt), nrow = 1, ncol = 1)
dev.off()

#-----------------------------------------------------------------------------------------------------------------#
# End of script
#-----------------------------------------------------------------------------------------------------------------#

#=================================================================================================================#
# R session info
#=================================================================================================================#

#> devtools::session_info()
#─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#setting  value                       
#version  R version 3.5.1 (2018-07-02)
#os       macOS  10.14.2              
#system   x86_64, darwin15.6.0        
#ui       RStudio                     
#language (EN)                        
#collate  en_US.UTF-8                 
#ctype    en_US.UTF-8                 
#tz       Europe/Berlin               
#date     2019-03-12                  
#
#─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#package       * version    date       lib source                           
#annotate        1.60.0     2018-10-30 [1] Bioconductor                     
#AnnotationDbi   1.44.0     2018-10-30 [1] Bioconductor                     
#assertthat      0.2.0      2017-04-11 [1] CRAN (R 3.5.0)                   
#backports       1.1.2      2017-12-13 [1] CRAN (R 3.5.0)                   
#bibtex          0.4.2      2017-06-30 [1] CRAN (R 3.5.0)                   
#bindr           0.1.1      2018-03-13 [1] CRAN (R 3.5.0)                   
#bindrcpp      * 0.2.2      2018-03-29 [1] CRAN (R 3.5.0)                   
#Biobase         2.42.0     2018-10-30 [1] Bioconductor                     
#BiocGenerics    0.28.0     2018-10-30 [1] Bioconductor                     
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
#dplyr         * 0.7.8      2018-11-10 [1] CRAN (R 3.5.0)                   
#evaluate        0.12       2018-10-09 [1] CRAN (R 3.5.0)                   
#fitdistrplus    1.0-11     2018-09-10 [1] CRAN (R 3.5.0)                   
#fs              1.2.6      2018-08-23 [1] CRAN (R 3.5.0)                   
#future          1.10.0     2018-10-17 [1] CRAN (R 3.5.0)                   
#future.apply    1.0.1      2018-08-26 [1] CRAN (R 3.5.0)                   
#gbRd            0.4-11     2012-10-01 [1] CRAN (R 3.5.0)                   
#gdata           2.18.0     2017-06-06 [1] CRAN (R 3.5.0)                   
#geneplotter     1.60.0     2018-10-30 [1] Bioconductor                     
#ggplot2       * 3.1.0      2018-10-25 [1] CRAN (R 3.5.0)                   
#ggrepel         0.8.0      2018-05-09 [1] CRAN (R 3.5.0)                   
#ggridges        0.5.1      2018-09-27 [1] CRAN (R 3.5.0)                   
#globals         0.12.4     2018-10-11 [1] CRAN (R 3.5.0)                   
#glue            1.3.0      2018-07-17 [1] CRAN (R 3.5.0)                   
#gplots          3.0.1      2016-03-30 [1] CRAN (R 3.5.0)                   
#graph           1.60.0     2018-10-30 [1] Bioconductor                     
#gridExtra     * 2.3        2017-09-09 [1] CRAN (R 3.5.0)                   
#GSEABase        1.44.0     2018-10-30 [1] Bioconductor                     
#GSVA            1.30.0     2018-10-30 [1] Bioconductor                     
#gtable          0.2.0      2016-02-26 [1] CRAN (R 3.5.0)                   
#gtools          3.8.1      2018-06-26 [1] CRAN (R 3.5.0)                   
#highr           0.7        2018-06-09 [1] CRAN (R 3.5.0)                   
#htmltools       0.3.6      2017-04-28 [1] CRAN (R 3.5.0)                   
#htmlwidgets     1.3        2018-09-30 [1] CRAN (R 3.5.0)                   
#httpuv          1.4.5      2018-07-19 [1] CRAN (R 3.5.0)                   
#httr            1.4.0      2018-12-11 [1] CRAN (R 3.5.1)                   
#ica             1.0-2      2018-05-24 [1] CRAN (R 3.5.0)                   
#igraph          1.2.2      2018-07-27 [1] CRAN (R 3.5.0)                   
#IRanges         2.16.0     2018-10-30 [1] Bioconductor                     
#irlba           2.3.2      2018-01-11 [1] CRAN (R 3.5.0)                   
#jsonlite        1.6        2018-12-07 [1] CRAN (R 3.5.0)                   
#KernSmooth      2.23-15    2015-06-29 [1] CRAN (R 3.5.1)                   
#knitr           1.20       2018-02-20 [1] CRAN (R 3.5.0)                   
#labeling        0.3        2014-08-23 [1] CRAN (R 3.5.0)                   
#later           0.7.5      2018-09-18 [1] CRAN (R 3.5.0)                   
#lattice         0.20-38    2018-11-04 [1] CRAN (R 3.5.0)                   
#lazyeval        0.2.1      2017-10-29 [1] CRAN (R 3.5.0)                   
#listenv         0.7.0      2018-01-21 [1] CRAN (R 3.5.0)                   
#lmtest          0.9-36     2018-04-04 [1] CRAN (R 3.5.0)                   
#lsei            1.2-0      2017-10-23 [1] CRAN (R 3.5.0)                   
#magrittr        1.5        2014-11-22 [1] CRAN (R 3.5.0)                   
#MASS            7.3-51.1   2018-11-01 [1] CRAN (R 3.5.1)                   
#Matrix          1.2-15     2018-11-01 [1] CRAN (R 3.5.1)                   
#matrixStats     0.54.0     2018-07-23 [1] CRAN (R 3.5.0)                   
#memoise         1.1.0      2017-04-21 [1] CRAN (R 3.5.0)                   
#metap           1.0        2018-07-25 [1] CRAN (R 3.5.0)                   
#mime            0.6        2018-10-05 [1] CRAN (R 3.5.0)                   
#munsell         0.5.0      2018-06-12 [1] CRAN (R 3.5.0)                   
#npsurv          0.4-0      2017-10-14 [1] CRAN (R 3.5.0)                   
#openxlsx      * 4.1.0      2018-05-26 [1] CRAN (R 3.5.0)                   
#outliers        0.14       2011-01-24 [1] CRAN (R 3.5.0)                   
#pbapply         1.3-4      2018-01-10 [1] CRAN (R 3.5.0)                   
#pbmcapply       1.3.0      2018-11-05 [1] CRAN (R 3.5.0)                   
#pheatmap        1.0.12     2019-01-04 [1] CRAN (R 3.5.2)                   
#pillar          1.3.0      2018-07-14 [1] CRAN (R 3.5.0)                   
#pkgbuild        1.0.2      2018-10-16 [1] CRAN (R 3.5.0)                   
#pkgconfig       2.0.2      2018-08-16 [1] CRAN (R 3.5.0)                   
#pkgload         1.0.2      2018-10-29 [1] CRAN (R 3.5.0)                   
#plotly          4.8.0      2018-07-20 [1] CRAN (R 3.5.1)                   
#plyr            1.8.4      2016-06-08 [1] CRAN (R 3.5.0)                   
#png             0.1-7      2013-12-03 [1] CRAN (R 3.5.0)                   
#prettyunits     1.0.2      2015-07-13 [1] CRAN (R 3.5.0)                   
#processx        3.2.1      2018-12-05 [1] CRAN (R 3.5.1)                   
#promises        1.0.1      2018-04-13 [1] CRAN (R 3.5.0)                   
#ps              1.2.1      2018-11-06 [1] CRAN (R 3.5.1)                   
#purrr           0.2.5      2018-05-29 [1] CRAN (R 3.5.0)                   
#R.methodsS3     1.7.1      2016-02-16 [1] CRAN (R 3.5.0)                   
#R.oo            1.22.0     2018-04-22 [1] CRAN (R 3.5.0)                   
#R.utils         2.7.0      2018-08-27 [1] CRAN (R 3.5.0)                   
#R6              2.3.0      2018-10-04 [1] CRAN (R 3.5.0)                   
#RANN            2.6        2018-07-16 [1] CRAN (R 3.5.0)                   
#RColorBrewer    1.1-2      2014-12-07 [1] CRAN (R 3.5.0)                   
#Rcpp            1.0.0      2018-11-07 [1] CRAN (R 3.5.0)                   
#RCurl           1.95-4.11  2018-07-15 [1] CRAN (R 3.5.0)                   
#Rdpack          0.10-1     2018-10-04 [1] CRAN (R 3.5.0)                   
#remotes         2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
#reshape2        1.4.3      2017-12-11 [1] CRAN (R 3.5.0)                   
#reticulate      1.10       2018-08-05 [1] CRAN (R 3.5.0)                   
#rlang           0.3.0.1    2018-10-25 [1] CRAN (R 3.5.0)                   
#rmarkdown       1.10       2018-06-11 [1] CRAN (R 3.5.0)                   
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
#shiny           1.2.0      2018-11-02 [1] CRAN (R 3.5.0)                   
#shinythemes     1.1.2      2018-11-06 [1] CRAN (R 3.5.0)                   
#SingleR       * 0.2.1      2018-12-12 [1] local                            
#stringi         1.2.4      2018-07-20 [1] CRAN (R 3.5.0)                   
#stringr         1.3.1      2018-05-10 [1] CRAN (R 3.5.0)                   
#survival        2.43-3     2018-11-26 [1] CRAN (R 3.5.0)                   
#tibble          1.4.2      2018-01-22 [1] CRAN (R 3.5.0)                   
#tidyr           0.8.2      2018-10-28 [1] CRAN (R 3.5.0)                   
#tidyselect      0.2.5      2018-10-11 [1] CRAN (R 3.5.0)                   
#tsne            0.1-3      2016-07-15 [1] CRAN (R 3.5.0)                   
#umap            0.2.0.0    2018-07-25 [1] CRAN (R 3.5.1)                   
#usethis         1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
#viridisLite     0.3.0      2018-02-01 [1] CRAN (R 3.5.0)                   
#withr           2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
#XML             3.98-1.16  2018-08-19 [1] CRAN (R 3.5.0)                   
#xtable          1.8-3      2018-08-29 [1] CRAN (R 3.5.0)                   
#yaml            2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
#zip             1.0.0      2017-04-25 [1] CRAN (R 3.5.0)                   
#zoo             1.8-4      2018-09-19 [1] CRAN (R 3.5.0)                   