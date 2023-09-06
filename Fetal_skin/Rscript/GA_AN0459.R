#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a healthy donor
# Date: 2020.07.22
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger
# - dimensionality reduction, clustering and biological annotation for clusters using customized references
# - exclude non-T cells from correlation
# Gattinoni L, Lugli E, Ji Y, Pos Z et al. A human memory T cell subset with stem cell-like properties. 
# Nat Med 2011 Sep 18;17(10):1290-7. PMID: 21926977 (BioProject: PRJNA131239)
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(openxlsx)
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

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid <- c("EX0008")

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
an.desc <- c("sc_tp_skin_blood_10x_tscm_tcf1_expression")
an.descs <- c("sc_tp_skin_blood")

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
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/refs"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/pca"))
dir.figs <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
#dir.figs.bd <- paste0("figures/", analysisid, "_", an.desc, "/blood/", analysisid, "_", an.desc)
#dir.figs.sk <- paste0("figures/", analysisid, "_", an.desc, "/skin/", analysisid, "_", an.descs)
#dir.figs.rf <- paste0("figures/", analysisid, "_", an.desc, "/refs/", analysisid, "_", an.descs)
#dir.figs.pca <- paste0("figures/", analysisid, "_", an.desc, "/pca/", analysisid, "_", an.descs)

#=================================================================================================================#
# Analysis of processed data for skin sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Load datasets
#-----------------------------------------------------------------------------------------------------------------#

# Load blood Seurat object from healthy donor and transplanted patient
pd <- list()
#pd[['hd.bd']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_seurat_blood.rds'))
#pd[['hd.sk']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_seurat_skin.rds'))
#pd$hd.bd@meta.data[6:19] <- NULL
#pd$hd.bd@meta.data[['sample']] <- c('Healthy blood')
#pd$hd.sk@meta.data[['sample']] <- c('Healthy skin')
#pd$hd.bd@meta.data$cluster <- pd$hd.bd@meta.data$RNA_snn_res.1
#pd$hd.sk@meta.data$cluster <- pd$hd.sk@meta.data$RNA_snn_res.1
#pd$hd.bd@meta.data[['RNA_snn_res.1']] <- NULL
#pd$hd.sk@meta.data[['RNA_snn_res.1']] <- NULL

# Load blood Seurat object from healthy donor and transplanted patient
pd[['tp.bd']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/blood/GA_AN0163_sc_tp_blood_seurat.rds'))
pd[['tp.sk']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/skin/GA_AN0163_sc_tp_skin_seurat.rds'))
pd$tp.bd@meta.data[c(5,7)] <- NULL
pd$tp.sk@meta.data[c(6)] <- NULL
pd$tp.bd@meta.data[['sample']] <- c('HSCT blood')
pd$tp.sk@meta.data[['sample']] <- c('HSCT skin')
pd$tp.bd@meta.data$cluster <- pd$tp.bd@meta.data$RNA_snn_res.1.51
pd$tp.sk@meta.data$cluster <- pd$tp.sk@meta.data$RNA_snn_res.1.2
pd$tp.bd@meta.data[['RNA_snn_res.1.51']] <- NULL
pd$tp.sk@meta.data[['RNA_snn_res.1.2']] <- NULL

# Load singleR object from healthy donor and transplanted patient
sr <- list()
#sr[['hd.bd']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_singler_blood.rds'))
#sr[['hd.sk']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_singler_skin.rds'))
sr[['tp.bd']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/blood/GA_AN0163_sc_tp_blood_singler.rds'))
sr[['tp.sk']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/skin/GA_AN0163_sc_tp_skin_singler.rds'))
# Transfer singleR cell type annotation to Seurat objects
#pd$hd.bd@meta.data[['celltype.sr']] <- sr$hd.bd[['singler']][[2]][['SingleR.single']][['labels']]
#pd$hd.sk@meta.data[['celltype.sr']] <- sr$hd.sk[['singler']][[2]][['SingleR.single']][['labels']]
pd$tp.bd@meta.data[['celltype.sr']] <- sr$tp.bd[['singler']][[2]][['SingleR.single']][['labels']]
pd$tp.sk@meta.data[['celltype.sr']] <- sr$tp.sk[['singler']][[2]][['SingleR.single']][['labels']]

# Annotate cell types to general cell types
#pd$hd.bd$subset <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), 
#                                   to = c('CD4 Tn','CD4 Tm','CD4 Tm','CD8 Tn','NK cell','CD8 Tm','NKT',
#                                          'B cell','gdT','CD4 Tn','B cell','Treg','n.d.'))
#pd$hd.bd$celltype <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), 
#                                     to = c('CD4 T','CD4 T','CD4 T','CD8 T','NK cell','CD8 T','NKT',
#                                            'B cell','gdT','CD4 T','B cell','Treg','n.d.'))
#pd$hd.bd$celltype <- factor(pd$hd.bd$celltype, levels = levels(pd$hd.bd$celltype)[c(1,2,6,7,4,3,5,8)])
pd$tp.bd$subset <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), 
                                   to = c('CD4 Tn','CD8 Tn','B cell','CD4 Tm','CD8 Tm','NK cell','gdT',
                                          'CD4 Tn','gdT','CD8 Tn','B cell','B cell','Treg'))
pd$tp.bd$celltype <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), 
                                     to = c('CD4 T','CD8 T','B cell','CD4 T','CD8 T','NK cell','gdT',
                                            'CD4 T','gdT','CD8 T','B cell','B cell','Treg'))
pd$tp.bd$celltype <- factor(pd$tp.bd$celltype, levels = levels(pd$tp.bd$celltype)[c(1,2,5,6,4,3)])
#pd$hd.sk$subset <- plyr::mapvalues(pd$hd.sk$cluster, from = levels(factor(pd$hd.sk$cluster)), 
#                                   to = c('CD4 T', 'CD4 T', 'CD4 T', 'CD4 T', 'CD8 T', 
#                                          'CD4 T', 'gdT', 'gdT', 'CD4 T', 'Treg', 'NK cell', 'CD8 T', 'CD8 T'))
#pd$hd.sk$subset <- factor(pd$hd.sk$subset, levels = levels(pd$hd.sk$subset)[c(1,2,3,4,5)])
#pd$hd.sk$celltype <- pd$hd.sk$subset
pd$tp.sk$subset <- plyr::mapvalues(pd$tp.sk$cluster, from = levels(factor(pd$tp.sk$cluster)), 
                                   to = c('gdT', 'CD4 T', 'CD4 T', 'CD4 T', 'gdT', 'CD4 T', 'CD8 T', 'CD8 T', 
                                          'CD8 T', 'NK cell','Treg', 'NKT', 'B cell'))
pd$tp.sk$subset <- factor(pd$tp.sk$subset, levels = levels(pd$tp.sk$subset)[c(2,3,1,5,6,4,7)])
pd$tp.sk$celltype <- pd$tp.sk$subset

# Load barcodes from host genotype in skin and bld
bc <- list(tp.sk = gsub('-1', '', read.delim('../201902/figures/GA_AN0162_sc_skin_blood_genotyping/GA_AN0162_sc_skin_genotyping_barcodes_host.tsv',as.is = T,header = F)[,1]), 
           tp.bd = gsub('-1', '', readLines('../201907/figures/GA_AN0260_sc_tp_blood_genotyping/GA_AN0260_sc_tp_blood_genotyping_barcodes_blood_host_cells.txt')))

# Keep only the barcodes present in seurat object after fitlering steps
bc$tp.sk <- bc$tp.sk[bc$tp.sk %in% names(Idents(pd$tp.sk))]
bc$tp.bd <- bc$tp.bd[bc$tp.bd %in% names(Idents(pd$tp.bd))]

# Classify idents according to barcodes from host and donor genotypes
pd$tp.sk[['genotype']] <- names(Idents(pd$tp.sk))
pd$tp.sk[['genotype']] <- ifelse(pd$tp.sk$genotype %in% bc$tp.sk, 'host', 'donor')
pd$tp.bd[['genotype']] <- names(Idents(pd$tp.bd))
pd$tp.bd[['genotype']] <- ifelse(pd$tp.bd$genotype %in% bc$tp.bd, 'host', 'donor')
#pd$hd.sk[['genotype']] <- c('host')
#pd$hd.bd[['genotype']] <- c('host')

# # Merge data sets
# tp.int <- merge(x = pd$tp.bd, y = pd$tp.sk, add.cell.ids = c("bd", "sk"), project = "HSCT patient")
# #pd$tp.mg$geno.ident[is.na(pd$tp.mg$geno.ident)] <- c("healthy (skin)")
# #pd$tp.mg <- ScaleData(pd$tp.mg, features = rownames(pd$tp.mg), vars.to.regress = c("nCount_RNA", "percent.mito"), assay = 'RNA')
# #saveRDS(pd$tp.mg, file = paste0(dir.figs, "_seurat_merged.rds"))
# #pd$tp.mg <- readRDS(paste0(dir.figs, "_seurat_merged.rds"))
# 
# # Annotate blood and skin subsets
# tp.int$celltypef <- paste0(ifelse(tp.int$orig.ident=='blood','b','s'), tp.int$celltype)
# #tp.int$subsetf <- factor(tp.int$subsetf, levels = levels(factor(tp.int$subsetf))[c(3,2,5,4,6,9,8,7,1)])
# tp.int$celltypef <- factor(tp.int$celltypef, levels = levels(factor(tp.int$celltypef))[c(2,3,4,6,5,1,8,9,10,13,12,11,7)])
# #VlnPlot(tp.int, features = c('GZMB'), ncol = 2, group.by = 'sample')
# table(tp.int$celltypef)
# 
# # Annotate blood and skin subsets split between memory and naive
# tp.int$subsetf <- paste0(ifelse(tp.int$orig.ident=='blood','b','s'), gsub('T','T',tp.int$subset))
# #tp.int$subsetf <- factor(tp.int$subsetf, levels = levels(factor(tp.int$subsetf))[c(3,2,5,4,6,9,8,7,1)])
# tp.int$subsetf <- factor(tp.int$subsetf, levels = levels(factor(tp.int$subsetf))[c(3,2,5,4,6,8,7,1,10,11,12,15,14,13,9)])
# #VlnPlot(tp.int, features = c('GZMB'), ncol = 2, group.by = 'sample')
# table(tp.int$subsetf)
# tp.int$genotype.subset <- paste(tp.int$genotype, tp.int$subset)
# 
# # Set clusters names
# tp.int[['cluster_split']] <- tp.int$cluster
# tp.int$cluster_split <- c(paste0('blood_cl', tp.int$cluster_split[tp.int$orig.ident == "blood"]), 
#                           paste0('skin_cl', tp.int$cluster_split[tp.int$orig.ident == "skin"]))
# tp.int$cluster_split <- factor(tp.int$cluster_split, levels = levels(factor(tp.int$cluster_split))[c(1,2,6:13,3,4,5,14,15,19:26,16,17,18)])
# tp.int$donor <- factor('HSCT')
# 
# # Extract tSNE coordenates
# tp.int <- FindVariableFeatures(tp.int, selection.method = "vst", nfeatures = 8000, 
#                                mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
# tp.int <- ScaleData(tp.int, features = rownames(tp.int))#, vars.to.regress = c("nCount_RNA", "percent.mito"))
# tp.int <- RunPCA(tp.int, features = VariableFeatures(tp.int), verbose = F)
# tp.int <- RunTSNE(tp.int, reduction = "pca", dims = 1:20, seed.use = 11, tsne.method = "Rtsne", dim.embed = 2, perplexity = 30, reduction.name = "tsne")
# tp.int@reductions$tsne@cell.embeddings <- rbind(pd$tp.bd@reductions$tsne@cell.embeddings, pd$tp.sk@reductions$tsne@cell.embeddings)
# rownames(tp.int@reductions$tsne@cell.embeddings) <- colnames(tp.int)

# Change cell names
#rownames(pd$tp.bd@meta.data) <- paste0('bd_',rownames(pd$tp.bd@meta.data))
#rownames(pd$tp.sk@meta.data) <- paste0('sk_',rownames(pd$tp.sk@meta.data))
#names(Idents(pd$tp.bd)) <- paste0('bd_',names(Idents(pd$tp.bd)))
#saveRDS(pd$tp.mg, file = paste0(dir.figs, "_seurat_merged.rds"))
#saveRDS(tp.int, file = paste0(dir.figs, "_seurat_merged.rds"))
pd$tp.mg <- readRDS(file = paste0("figures/GA_AN0453_sc_tp_skin_blood_10x_correlation_tcomp/GA_AN0453_sc_tp_skin_blood_seurat_merged.rds"))

#-----------------------------------------------------------------------------------------------------------------#
# Single cell correlation with bulk transcriptomes
#-----------------------------------------------------------------------------------------------------------------#

# # Load bulk transcirptomes
# pathdir <- list.dirs("../202007/figures/GA_AN0452_Tcompartments_PRJNA245319")
# files <- list.files(pathdir, ".csv", full.names = T)
# ref <- Map(x=files, function(x)  read.csv(x, stringsAsFactors = F) )
# names(ref) <- stringr::str_split(names(ref), "expression_|\\.csv", simplify = T)[,2]
# ref <- ref$TCM_vs_TN[,c(20,23:34)]
# ref <- ref[!is.na(ref$genename), ]
# ref <- ref[!(ref$genename == ''), ]
# ref <- ref[order(rowSums(ref[,2:13]), decreasing = T), ]
# ref <- ref[!duplicated(ref$genename), ]
# row.names(ref) <- ref$genename
# ref <- ref[,-1]
# 
# # Create SingleR reference object
# main.1 <- paste(stringr::str_split(colnames(ref), "_\\d", n = 2, simplify = T)[, 1])
# type.1 <- stringr::str_split(colnames(ref), "_\\d", n = 2, simplify = T)[, 1]
# #ref.1 <- list(name = "cheuk", data = as.matrix(log2(ref+1)), types = type.1, main_types = main.1)
# ref.1 <- list(name = "gattinoni", data = as.matrix(ref), types = type.1, main_types = main.1)
# #ref.1 <- list(name = "zielinski_skin", data = as.matrix(log2(ref.1[-25]+1)), types = type, main_types = main)
# # If using the de method, predefine the variable genes
# ref.1$de.genes <- CreateVariableGeneSet(ref.1$data, type.1, 500)#200)
# ref.1$de.genes.main <- CreateVariableGeneSet(ref.1$data, main.1, 500)#300)
# # If using the sd method, define an sd threshold
# ref.1$sd.thres <- sort(matrixStats::rowSds(ref.1$data), decreasing = T)[10000]#[4000] # or any other threshold
# #saveRDS(ref, file = paste0(dir.figs.rf, "_singler_ref_cheuk.rds"))
# 
# # Correlation of variable genes between the experiment and the reference dataset
# SingleR.DrawScatter(sc_data = GetAssayData(tp.int, slot = "data"), cell_id = 1, ref = ref.1, sample_id = 1)
# # Correlation between a single cell from the experiment and all reference cells
# SingleR.DrawBoxPlot(sc_data = GetAssayData(tp.int, slot = "data"), cell_id = 1, ref = ref.1, main_types = F, labels.use = NULL)$p
# SingleR(sc_data = as.matrix(GetAssayData(tp.int, slot = "data")[, c(1, 1)]), 
#                 ref_data = ref.1$data, types = type.1, fine.tune = F, 
#                 sd.thres = ref.1$sd.thres, do.pvals = F)

# Create the SingleR object using clusters calculated with Seurat
#sr.tp <- CreateSinglerObject(GetAssayData(tp.int, slot = "data"), 
#                             annot = NULL, project.name = "ann_gattinoni", min.genes = 200, #500
#                             technology = "10X", species = "Human", citation = "", normalize.gene.length = F, 
#                             variable.genes = "de", ref.list = list(ref.1), #list(hpca, blueprint_encode, ref), 
#                             fine.tune = T, do.signatures = F, do.main.types = F, clusters = tp.int$cluster_split, 
#                             reduce.file.size = T, numCores = 4)#SingleR.numCores)
#saveRDS(sr.tp, file = paste0(dir.figs, "_singler_gattinoni.rds"))
sr.tp <- readRDS(file = paste0("figures/GA_AN0453_sc_tp_skin_blood_10x_correlation_tcomp/GA_AN0453_sc_tp_skin_blood_singler_gattinoni.rds"))
table(sr.tp[["singler"]][[1]][["SingleR.single"]][["labels"]], sr.tp[["singler"]][[1]][["SingleR.single"]][["labels1"]])
# sr.tp$meta.data$orig.ident <- tp.int@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.tp$meta.data$xy <- rbind(pd$tp.bd@reductions$tsne@cell.embeddings, pd$tp.sk@reductions$tsne@cell.embeddings)#tp.int@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "tsne"))
# identical(sub('bd_|sk_','',colnames(tp.int)), rownames(sr.tp$meta.data$xy))
rownames(sr.tp$meta.data$xy) <- colnames(tp.int)
sr.tp$meta.data$xy.um <- rbind(pd$tp.bd@reductions$tsne@cell.embeddings, pd$tp.sk@reductions$tsne@cell.embeddings)#tp.int@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "umap"))
rownames(sr.tp$meta.data$xy.um) <- colnames(tp.int)
# #sr.tp$meta.data$clusters <- sk@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(sk))
# sr.tp$meta.data$clusters <- tp.int$cluster_split
# sr.tp$meta.data$celltypef <- tp.int$celltypef
# #sr.tp$singler[[1]]$SingleR.single$labels[,1] <- factor(sr.tp$singler[[1]]$SingleR.single$labels[,1], 
# #                                                       levels = levels(factor(sort(sr.tp$singler[[1]]$SingleR.single$labels[,1])))[c(3,4,1,2)])

# Filter labels to remove non-T cells
#sr.tp$singler[[1]]$SingleR.single$labels[,1] <- ifelse(tp.int@meta.data$celltype %in% c('B cell','NK cell'), 'other', sr.tp$singler[[1]]$SingleR.single$labels[,1])
sr.tp$singler[[1]]$SingleR.single$labels[,1] <- ifelse(pd$tp.mg@meta.data$celltype %in% c('B cell','NK cell'), 'other', sr.tp$singler[[1]]$SingleR.single$labels[,1])

# Transfer singleR cell type annotation to Seurat objects
pd$tp.mg@meta.data[['Compartment']] <- sr.tp[['singler']][[1]][['SingleR.single']][['labels']]
pd$tp.mg@meta.data$Compartment <- factor(pd$tp.mg@meta.data$Compartment, levels = sort(as.character(unique(pd$tp.mg@meta.data$Compartment)))[c(4,5,2,3,1)])
table(pd$tp.mg@meta.data$Compartment)

#-----------------------------------------------------------------------------------------------------------------#
# Visualize gene expression
#-----------------------------------------------------------------------------------------------------------------#

# Define gene list
genelist <- sort(c('TCF7','TOX')) #,'CD127', 'CD95
all(genelist %in% rownames(tp.int)); genelist[!(genelist %in% rownames(tp.int))]

# Define colors
col1 <- c("darkseagreen2",'grey20', "firebrick3",'royalblue3', 'grey80')
#pal <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')
#blues <- RColorBrewer::brewer.pal(n = 9, name = 'Blues')
#reds <- RColorBrewer::brewer.pal(n = 9, name = 'Reds')
#col.ct <- c('black', pal[c(1,2,5,6,8,4)],'darkorchid3','grey80', 
#            #c(blues[c(3,5,7,9)],pal[c(5,6)],reds[c(8)],pal[c(4)],pal[c(7,8)],'#CD96CD','darkorchid3','grey40'))
#            pal[c(2,6,8,4)],'#CD96CD','darkorchid3','grey80')
#col.ss <- col.ct[-c(1,2,4)]
#col.hm2 <- grDevices::colorRampPalette(c(rev(RColorBrewer::brewer.pal(n = 9, name = 'YlGnBu')), (RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[2:9])))(50)
col <- list(cluster = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(1:10,1:10)], 
            subset = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(1:10,1:10)])
#col$hm <- colorRampPalette(c('#440154', '#31688E', '#35B779', '#E8E419'))(50)
col$hm <- colorRampPalette(c('grey80', 'skyblue2', 'royalblue1', 'royalblue4'))(50)
col$tissue <- setNames(c('white','grey60'), levels(pd$tp.mg$orig.ident))
col$genotype <- setNames(c('white','grey60'), levels(pd$tp.mg$genotype))

# Extract data
df.mk <- Map(x=names(pd)[3], function(x) { 
  df <- data.frame(tx = c(pd$tp.bd@reductions$tsne@cell.embeddings[,1], pd$tp.sk@reductions$tsne@cell.embeddings[,1]),
                   ty = c(pd$tp.bd@reductions$tsne@cell.embeddings[,2], pd$tp.sk@reductions$tsne@cell.embeddings[,2]),
                   #tx = sr.tp$meta.data$xy[,1],
                   #ty = sr.tp$meta.data$xy[,2], 
                   #tx = pd[[x]]@reductions$tsne@cell.embeddings[,1], 
                   #ty = pd[[x]]@reductions$tsne@cell.embeddings[,2], 
                   sample = pd[[x]]@meta.data$orig.ident, 
                   celltype = pd[[x]]@meta.data$celltype, 
                   cells = colnames(pd[[x]]), 
                   genotype = pd[[x]]@meta.data$genotype, 
                   subset = pd[[x]]@meta.data$subsetf, 
                   #subset.genotype = pd[[x]]$subset.genotype,
                   compartment = pd[[x]]@meta.data$Compartment, 
                   t(as.matrix(GetAssayData(pd[[x]], slot = "data")[genelist,,drop=F])))
  df <- reshape2::melt(df, id.vars = c("tx", "ty", "sample", "celltype", "cells", "genotype", "subset",'compartment')) 
  df$value_noise <- df$value + rnorm(n = length(df$value)) / 100000 #add noise 
  return(df) })
df.mk$tp.mg$celltype <- factor(df.mk$tp.mg$celltype, levels = levels(factor(df.mk$tp.mg$celltype))[c(2,3,4,7,6,5,1)])
#df.mk$tp.mg$compartment <- factor(df.mk$tp.mg$compartment, levels = levels(factor(df.mk$tp.mg$compartment))[c()])

# Plot expression over tSNE
ts.mk <- Map(x=names(pd)[3], function(x) { lapply(setNames(genelist, nm = genelist), function(y) { 
  lapply(setNames(c('blood', 'skin'), nm = c('blood', 'skin')), function(z) {
    ggplot(df.mk[[x]][df.mk[[x]]$variable == y & df.mk[[x]]$sample == z,] %>% arrange((value)), aes(x = tx, y = ty, color = value)) + 
      geom_point(shape = 16, size = ifelse(z=='blood',0.8,1.2), alpha = 1) + 
      #scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(min(df.mk[[k]][x,]$value),3)) + 
      scale_color_gradientn('Expression', colours = col$hm, oob = scales::squish, 
                            limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
      scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
      #facet_wrap(~variable+genotype, ncol = 4) + 
      theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
      ggtitle(paste0('HSCT ',z, ' (lymphocytes)'), subtitle = paste0("Gene expression: ", y)) })})})
#grid.arrange(grobs = ts.mk$sk, ncol = 4)
ts.mkt <- Map(x=names(pd)[3], function(x) { lapply(setNames(genelist, nm = genelist), function(y) {
    ggplot(df.mk[[x]][df.mk[[x]]$variable == y,] %>% arrange((value)), aes(x = tx, y = ty, color = value)) + 
      geom_point(shape = 16, size = 0.8, alpha = 1) + 
      #scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(min(df.mk[[k]][x,]$value),3)) + 
      scale_color_gradientn('Expression', colours = col$hm, oob = scales::squish, 
                            limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
      scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + facet_wrap(~sample, ncol = 2) + 
      theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
      ggtitle(paste0('HSCT patient', ' (lymphocytes)'), subtitle = paste0("Gene expression: ", y)) })})

# Plot expression over tSNE (different colors for host and donor)
col$hm2 <- colorRampPalette(c('grey80', 'lightsalmon1', 'firebrick1', 'firebrick4'))(50)
ts.mks <- Map(x=names(pd)[3], function(x) { lapply(setNames(genelist, nm = genelist), function(y) { 
  lapply(setNames(c('blood', 'skin'), nm = c('blood', 'skin')), function(z) {
    ggplot(mapping = aes(x = tx, y = ty, color = value, fill = value)) + 
      geom_point(data = df.mk[[x]][df.mk[[x]]$variable == y & df.mk[[x]]$sample == z & df.mk[[x]]$genotype == 'donor',] %>% arrange((value)), 
                 shape = 21, size = ifelse(z=='blood',1.2,1.6), alpha = 1, stroke = 0, color='grey80') + 
      scale_fill_gradientn('Donor', colours = col$hm2, oob = scales::squish, limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
      #new_scale_color() + 
      geom_point(data = df.mk[[x]][df.mk[[x]]$variable == y & df.mk[[x]]$sample == z & df.mk[[x]]$genotype == 'host',] %>% arrange((value)), 
                 shape = 16, size = ifelse(z=='blood',0.8,1.2), alpha = 1) + 
      #scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(min(df.mk[[k]][x,]$value),3)) + 
      scale_color_gradientn('Host', colours = col$hm, oob = scales::squish, limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
      scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
      #facet_wrap(~variable+genotype, ncol = 4) + 
      theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
      ggtitle(paste0('HSCT ',z, ' (lymphocytes)'), subtitle = paste0("Gene expression: ", y)) })})})

#-----------------------------------------------------------------------------------------------------------------#
# Differential gene expression between blood and skin subsets
#-----------------------------------------------------------------------------------------------------------------#

# Split compartment annotation between tissues
#pd$tp.mg$Compartmentt <- as.character(pd$tp.mg$Compartment)
#pd$tp.mg$Compartmentt[pd$tp.mg$orig.ident %in% 'blood'] <- sub('^','b',pd$tp.mg$Compartmentt[pd$tp.mg$orig.ident %in% 'blood'])
#pd$tp.mg$Compartmentt[pd$tp.mg$orig.ident %in% 'skin'] <- sub('^','s',pd$tp.mg$Compartmentt[pd$tp.mg$orig.ident %in% 'skin'])
#table(pd$tp.mg$Compartmentt)
pd$tp.bd@meta.data$Compartment <- pd$tp.mg$Compartment[pd$tp.mg$orig.ident %in% 'blood']
pd$tp.sk@meta.data$Compartment <- pd$tp.mg$Compartment[pd$tp.mg$orig.ident %in% 'skin']

# Find differentially expressed features
de <- Map(x=names(pd)[1:2], function(x) { #lapply(setNames(unique(pd[[x]]$orig.ident), nm = unique(pd[[x]]$orig.ident)), function(w) { 
  celltypes <- levels(pd[[x]]$Compartment)[-grep('other',levels(pd[[x]]$Compartment))]
  #celltypes <- unique(pd[[x]]$Compartmentt)[-grep('other',unique(pd[[x]]$Compartmentt))]
  lapply(setNames(celltypes, nm = celltypes), function(y) { 
    lapply(setNames(celltypes[!(celltypes %in% y)], nm = celltypes[!(celltypes %in% y)]), function(z) { #print(paste0(x,': ',y,' vs ',z))
      Idents(pd[[x]]) <- pd[[x]]@meta.data$Compartment
      de.int <- FindMarkers(pd[[x]], assay = 'RNA', ident.1 = y, ident.2 = z, 
                            #min.pct = 0, logfc.threshold = 0, test.use = 'wilcox', min.diff.pct = -Inf, random.seed = 1)
                            min.pct = 0.1, logfc.threshold = 0.25, test.use = 'wilcox', min.diff.pct = -Inf, random.seed = 1)
      de.int[, c('gene','target','ref','vs')] <- list(row.names(de.int), y, z, paste0(y, '_vs_', z))
      de.int %>% group_by(dn=avg_logFC<=0) %>% group_by(sig=p_val_adj<0.05) %>% #filter(p_val_adj < 0.05) %>% 
        arrange(p_val) %>% arrange(dn) }) }) }) #})
#saveRDS(de, file = paste0(dir.figs, "_diff_genes_compartments.rds"))
de <- readRDS(paste0(dir.figs, "_diff_genes_compartments.rds"))
#deb <- Map(x=names(de), function(x) bind_rows(lapply(setNames(names(de[[x]][[y]]), nm = names(de[[x]][[y]]))), .id = 'tissue') )
deb <- bind_rows(lapply(setNames(names(de), nm = c('blood','skin')), function(x) bind_rows(lapply(setNames(names(de[[1]]), nm = names(de[[1]])), function(y) 
  bind_rows(de[[x]][[y]]) )) ), .id = 'sample' )

# Plot expression as violin plots
vp.mk <- Map(x=names(pd)[3], function(x) { lapply(setNames(genelist, nm = genelist), function(y) { 
  ya <- max(df.mk[[x]][df.mk[[x]]$variable == y & !(df.mk[[x]]$compartment %in% 'other'),]$value)
  sig <- deb[deb$gene %in% y & deb$dn == FALSE,]
  #sig$p_val_adj <- as.numeric(formatC(sig$p_val_adj,format = 'e', digits = 1))
  #sig$sig <- ifelse(sig$value<0.05, '*', 'ns')
  sig$xs <- as.numeric(factor(sig$target, levels = levels(pd[[3]]$Compartment)[-5]))
  sig$xe <- as.numeric(factor(sig$ref, levels = levels(pd[[3]]$Compartment)[-5]))
  sig$yp <- c(if(sum(sig$sample == 'blood')>0) 1:nrow(sig[sig$sample == 'blood',]), if(sum(sig$sample == 'skin')>0) 1:nrow(sig[sig$sample == 'skin',]) )
  sig <- sig[(sig$sig),]
  #max <- ceiling(max(df.mk[[x]][df.mk[[x]]$variable == y,]$value))
  p<-ggplot(df.mk[[x]][df.mk[[x]]$variable == y & !(df.mk[[x]]$compartment %in% 'other'),],  aes(x = compartment, y = value, fill = compartment)) + 
    geom_jitter(shape = 16, size = 0.5, alpha = 0.5, position=position_jitter(width=0.2, height=0, seed=1), color = "black", show.legend = F) + 
    #geom_boxplot(outlier.colour = NA, color='grey40', lwd = 0.3, width=0.6, alpha = 0.4, show.legend = F) + 
    geom_violin(scale = "area", color = "black", size = 0.2, show.legend = F, alpha = 0.6) + 
    #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
    scale_x_discrete('') + 
    #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.1, label=paste0('*')), inherit.aes = F, size = 5) + 
    #scale_y_continuous(breaks = seq(0,10,1), limits = c(-0.01,ifelse(max<4,4,max))) + 
    scale_y_continuous(breaks = seq(0,5,1), limits = c(-0.01,5.01)) + 
    scale_fill_manual(values = col1) + xlab('') + ylab(expression(''*log[2]*'(norm. expression)')) + 
    ggtitle('HSCT patient - CD45+ cells', subtitle = paste0("Gene expression: ", y)) + facet_wrap(~sample, ncol = 4) + 
    #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
    theme_custom + theme(aspect.ratio = NULL, axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1))
  if(nrow(sig)>0) p + geom_text(data = sig, aes(y = 2.05 - 0.45*yp+(0.9*ya), x=(xs+xe-1)/2, label=paste0('',formatC(p_val_adj,format = 'e', digits = 1))), 
                                nudge_x = 0.5, size = 2.5, color = "black", inherit.aes = F) + 
    geom_segment(data = sig, colour = "black", show.legend = F, size=0.3, 
                 aes(x = xs-0.1, xend = xe+0.1, y = (1.83-0.45*yp+(0.9*ya)), yend = (1.83-0.45*yp+(0.9*ya))), inherit.aes = F) else p
  #if(dim(de$tp.mg[de$tp.mg$gene %in% y,])[1]!=0) 
  #  p+geom_text(data=de$tp.mg[de$tp.mg$gene %in% y & de$tp.mg$sig == T,], 
  #              aes(x=1.5,y=5.95, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size=3) else p
  ##annotate('text', x=1.5, y=5.05, label=de$tp.mg[de$tp.mg$gene %in% genelist,]$p_val_adj) + 
}) })

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Save tSNE plots for expression
lapply(names(ts.mk$tp.mg), function(x) lapply(names(ts.mk$tp.mg[[x]]), function(y) ggsave(plot = ts.mk$tp.mg[[x]][[y]], file = paste0(dir.figs, "_tsne_", x,'_',y, ".pdf"), height = 4, width = 4.5)) )
lapply(names(ts.mks$tp.mg), function(x) lapply(names(ts.mks$tp.mg[[x]]), function(y) ggsave(plot = ts.mks$tp.mg[[x]][[y]], file = paste0(dir.figs, "_tsne_", x,'_',y, "_split.pdf"), height = 4, width = 4.5)) )
lapply(names(ts.mkt$tp.mg), function(x) ggsave(plot = ts.mkt$tp.mg[[x]], file = paste0(dir.figs, "_tsne_", x, ".pdf"), height = 4.2, width = 8))

# Save violin plots
lapply(names(vp.mk$tp.mg), function(x) { ggsave(paste0(dir.figs, "_violin_", x, ".pdf"), vp.mk$tp.mg[[x]], height = 4.3, width = 5) })

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