#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a transplanted
# AML patient (chimeric sample)
# Date: 2020.10.15
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger
# - annotate clusters through expression of cell type markers
# - modify analysis GA_AN0305 to match cluster numbers and subsets for consistency between all samples
# - modify analysis GA_AN0353 to have a unique color for n.d. (instead of grey, that is also used for B cells)
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
library(grid)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
#library(stringr)
#library(reshape2)
#library(plyr)
#library(scales)
#library(VennDiagram)
#library(readr)
#library(vsn)
#library(edgeR)
#library(biomaRt)
#library(RColorBrewer)
#library(ggrepel)
#library(ggplotify)

# Define custom theme for ggplot
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
                      strip.text = element_text(size = 14, face = "plain", margin = margin(0,0,5,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid1 <- c('EX0004')
experimentid2 <- c('EX0008')

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])

# Define samples analyzed (all)
#sa <- c('cd')

# Define reference sample name
ct <- c('blood')
tr <- c('skin')

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c("sc_10x_hsct_consistent_annotation")
an.descs <- c("sc_hsct")

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx('../../ngs_sample_list.xlsm', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid1 | nsl$Experiment_ID == experimentid2, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], summarise, sum = NA)

# Update table containing the analyses list to include this analysis
# Include info like which read alignment algorithm was used
# Load, update avoiding duplicates (when script is run several times) and save back overwriting old one
nal <- read.xlsx('../../ngs_analysis_list.xlsm', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nal <- nal[nal$Analysis_ID == analysisid, ]

# Sample data
metadata <- data.frame(row.names = nsl$Sample_name, 
                       reshape2::colsplit(nsl$Sample_name, '_', c('exp', 'celltype', 'source')))

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
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/healthy_blood'))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/healthy_skin'))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/hsct_blood'))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/hsct_skin'))
dir.figs <- list(all = paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs) )
dir.figs$hd.bd <- paste0('figures/', analysisid, '_', an.desc, '/healthy_blood/', analysisid, '_', an.descs, '_healthy_blood')
dir.figs$hd.sk <- paste0('figures/', analysisid, '_', an.desc, '/healthy_skin/', analysisid, '_', an.descs, '_healthy_skin')
dir.figs$tp.bd <- paste0('figures/', analysisid, '_', an.desc, '/hsct_blood/', analysisid, '_', an.descs, '_hsct_blood')
dir.figs$tp.sk <- paste0('figures/', analysisid, '_', an.desc, '/hsct_skin/', analysisid, '_', an.descs, '_hsct_skin')
#dir.figs$hd.bd <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs, '_healthy_blood')
#dir.figs$hd.sk <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs, '_healthy_skin')
#dir.figs$tp.bd <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs, '_hsct_blood')
#dir.figs$tp.sk <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs, '_hsct_skin')

#=================================================================================================================#
# Analysis of processed data from scRNA-Seq
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Prepare data
#-----------------------------------------------------------------------------------------------------------------#

# Load blood Seurat object from healthy donor and transplanted patient
pd <- list()
pd[['hd.bd']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_seurat_blood.rds'))
pd[['hd.sk']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_seurat_skin.rds'))
pd$hd.bd@meta.data[6:19] <- NULL
pd$hd.bd@meta.data[['sample']] <- c('Healthy blood')
pd$hd.sk@meta.data[['sample']] <- c('Healthy skin')
pd$hd.bd@meta.data$cluster <- pd$hd.bd@meta.data$RNA_snn_res.1
pd$hd.sk@meta.data$cluster <- pd$hd.sk@meta.data$RNA_snn_res.1
pd$hd.bd@meta.data[['RNA_snn_res.1']] <- NULL
pd$hd.sk@meta.data[['RNA_snn_res.1']] <- NULL

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
sr[['hd.bd']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/blood/GA_AN0146_sc_skin_blood_10x_singler_blood.rds'))
sr[['hd.sk']] <- readRDS(paste0('../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_singler_skin.rds'))
sr[['tp.bd']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/blood/GA_AN0163_sc_tp_blood_singler.rds'))
sr[['tp.sk']] <- readRDS(paste0('../201902/figures/GA_AN0163_sc_skin_blood_10x_chimeric/skin/GA_AN0163_sc_tp_skin_singler.rds'))
# Transfer singleR cell type annotation to Seurat objects
pd$hd.bd@meta.data[['celltype.sr']] <- sr$hd.bd[['singler']][[2]][['SingleR.single']][['labels']]
pd$hd.sk@meta.data[['celltype.sr']] <- sr$hd.sk[['singler']][[2]][['SingleR.single']][['labels']]
pd$tp.bd@meta.data[['celltype.sr']] <- sr$tp.bd[['singler']][[2]][['SingleR.single']][['labels']]
pd$tp.sk@meta.data[['celltype.sr']] <- sr$tp.sk[['singler']][[2]][['SingleR.single']][['labels']]

# Annotate cell types to general cell types
pd$hd.bd$subset <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), 
                                   to = c('CD4 Tn','CD4 Tm','CD4 Tm','CD8 Tn','NK cell','CD8 Tm','NKT',
                                          'B cell','gdT','CD4 Tn','B cell','Treg','n.d.'))
pd$hd.bd$celltype <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), 
                                     to = c('CD4 T','CD4 T','CD4 T','CD8 T','NK cell','CD8 T','NKT',
                                            'B cell','gdT','CD4 T','B cell','Treg','n.d.'))
pd$hd.bd$celltype <- factor(pd$hd.bd$celltype, levels = levels(pd$hd.bd$celltype)[c(1,2,6,7,4,3,5,8)])
pd$tp.bd$subset <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), 
                                   to = c('CD4 Tn','CD8 Tn','B cell','CD4 Tm','CD8 Tm','NK cell','gdT',
                                          'CD4 Tn','gdT','CD8 Tn','B cell','B cell','Treg'))
pd$tp.bd$celltype <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), 
                                     to = c('CD4 T','CD8 T','B cell','CD4 T','CD8 T','NK cell','gdT',
                                            'CD4 T','gdT','CD8 T','B cell','B cell','Treg'))
pd$tp.bd$celltype <- factor(pd$tp.bd$celltype, levels = levels(pd$tp.bd$celltype)[c(1,2,5,6,4,3)])
pd$hd.sk$subset <- plyr::mapvalues(pd$hd.sk$cluster, from = levels(factor(pd$hd.sk$cluster)), 
                                   to = c('CD4 T', 'CD4 T', 'CD4 T', 'CD4 T', 'CD8 T', 
                                          'CD4 T', 'gdT', 'gdT', 'CD4 T', 'Treg', 'NK cell', 'CD8 T', 'CD8 T'))
pd$hd.sk$subset <- factor(pd$hd.sk$subset, levels = levels(pd$hd.sk$subset)[c(1,2,3,4,5)])
pd$hd.sk$celltype <- pd$hd.sk$subset
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
pd$hd.sk[['genotype']] <- c('host')
pd$hd.bd[['genotype']] <- c('host')

# # Get read counts data for each object
# ct <- list(hd.bd = GetAssayData(pd$hd.bd, slot = 'data'), 
#            hd.sk = GetAssayData(pd$hd.sk, slot = 'data'), 
#            tp.bd = GetAssayData(pd$tp.bd, slot = 'data'), 
#            tp.sk = GetAssayData(pd$tp.sk, slot = 'data'))
#ct <- Map(x=names(pd), function(x) GetAssayData(pd[[x]], slot = 'counts'))

#-----------------------------------------------------------------------------------------------------------------#
# Manual annotation
#-----------------------------------------------------------------------------------------------------------------#

# Rearrange cluster numbers to match subsets in all four samples
for(x in names(pd)) pd[[x]]$cluster_old <- pd[[x]]$cluster
#pd$hd.bd$cluster_new <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), to = c(2,3,4,5,9,6,8,12,10,1,11,7,13))
pd$hd.bd$cluster_new <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), to = c(2,3,4,5,10,6,9,12,8,1,11,7,13))
pd$hd.bd$cluster_new <- factor(pd$hd.bd$cluster_new, levels = sort(as.numeric(levels(factor(pd$hd.bd$cluster_new)))))
pd$hd.sk$cluster_new <- plyr::mapvalues(pd$hd.sk$cluster, from = levels(factor(pd$hd.sk$cluster)), to = c(2,1,3,4,7,5,12,11,6,10,13,9,8))
pd$hd.sk$cluster_new <- factor(pd$hd.sk$cluster_new, levels = sort(as.numeric(levels(factor(pd$hd.sk$cluster_new)))))
pd$tp.bd$cluster_new <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), to = c(2,4,11,3,6,10,8,1,9,5,12,13,7))
pd$tp.bd$cluster_new <- factor(pd$tp.bd$cluster_new, levels = sort(as.numeric(levels(factor(pd$tp.bd$cluster_new)))))
pd$tp.sk$cluster_new <- plyr::mapvalues(pd$tp.sk$cluster, from = levels(factor(pd$tp.sk$cluster)), to = c(9,1,4,2,10,3,5,6,7,12,8,11,13))
pd$tp.sk$cluster_new <- factor(pd$tp.sk$cluster_new, levels = sort(as.numeric(levels(factor(pd$tp.sk$cluster_new)))))
#for(x in names(pd)) pd[[x]]$subset_cluster <- paste0(pd[[x]]@meta.data$cluster_new,' - ',pd[[x]]@meta.data$celltype)
#pd$hd.bd$subset_cluster <- factor(pd$hd.bd$subset_cluster, levels = levels(factor(pd$hd.bd$subset_cluster))[c(1,6:13,2:5)])
#pd$hd.sk$subset_cluster <- factor(pd$hd.sk$subset_cluster, levels = levels(factor(pd$hd.sk$subset_cluster))[c(1,6:13,2:5)])
#pd$tp.bd$subset_cluster <- factor(pd$tp.bd$subset_cluster, levels = levels(factor(pd$tp.bd$subset_cluster))[c(1,6:13,2:5)])
#pd$tp.sk$subset_cluster <- factor(pd$tp.sk$subset_cluster, levels = levels(factor(pd$tp.sk$subset_cluster))[c(1,6:13,2:5)])
for(x in names(pd)) pd[[x]]$subset_cluster <- paste0(ifelse(as.numeric(pd[[x]]@meta.data$cluster)<10, paste0('  ',pd[[x]]@meta.data$cluster), paste0(pd[[x]]@meta.data$cluster)),' - ',pd[[x]]@meta.data$celltype)
for(x in names(pd)) pd[[x]]$subset_cluster_new <- paste0(ifelse(as.numeric(pd[[x]]@meta.data$cluster_new)<10, paste0('  ',pd[[x]]@meta.data$cluster_new), paste0(pd[[x]]@meta.data$cluster_new)),' - ',pd[[x]]@meta.data$celltype)

# Extract data
df <- Map(x=names(pd), function(x) { 
  data.frame(cells = colnames(pd[[x]]), 
             tx = pd[[x]]@reductions$tsne@cell.embeddings[,1], 
             ty = pd[[x]]@reductions$tsne@cell.embeddings[,2], 
             ux = pd[[x]]@reductions$umap@cell.embeddings[,1], 
             uy = pd[[x]]@reductions$umap@cell.embeddings[,2], 
             sample = pd[[x]]@meta.data$sample, 
             cluster = pd[[x]]@meta.data$cluster,
             #cluster_new = factor(pd[[x]]$cluster_new, levels = sort(as.numeric(levels(factor(pd[[x]]$cluster_new))))),
             cluster_new = pd[[x]]@meta.data$cluster_new,
             compartment = pd[[x]]@meta.data$subset, 
             subset = pd[[x]]@meta.data$celltype, 
             #subset_cluster = paste0(pd[[x]]@meta.data$cluster_new,' - ',pd[[x]]@meta.data$celltype),
             subset_cluster = pd[[x]]@meta.data$subset_cluster,
             subset_cluster_new = pd[[x]]@meta.data$subset_cluster_new, 
             genotype = pd[[x]]@meta.data$genotype, 
             sex = ifelse(as.matrix(GetAssayData(pd[[x]])["RPS4Y1",])>0, 'male', 'n.d.'), 
             row.names = 1:length(colnames(pd[[x]]))) })
lapply(df, dim)

# Position for cluster labels in the tSNE plot
#labels.loc <- Map(x=names(pd), function(x) { lapply(X = levels(factor(pd$hd.bd$cluster)), FUN = function(group) {
#  data.use <- df[[x]][df[[x]][, 'cluster'] == group, c('tx','ty')]
#  return(apply(X = data.use, MARGIN = 2, FUN = median, na.rm = TRUE)) }) })
labels.loc <- Map(x=names(pd), function(x) { 
  data.use <- lapply(X = levels(factor(pd[[x]]$cluster)), FUN = function(group) { 
    data.use <- df[[x]][df[[x]][, 'cluster'] == group, c('tx','ty')]
    return(apply(X = data.use, MARGIN = 2, FUN = median, na.rm = TRUE)) }) 
  names(data.use) <- levels(factor(pd[[x]]$cluster)) 
  data.frame(as.data.frame(x = t(x = as.data.frame(x = data.use))), cluster = levels(factor(pd[[x]]$cluster))) })
labels.loc_new <- Map(x=names(pd), function(x) { 
  data.use <- lapply(X = levels(factor(pd[[x]]$cluster_new)), FUN = function(group) { 
    data.use <- df[[x]][df[[x]][, 'cluster_new'] == group, c('tx','ty')]
    return(apply(X = data.use, MARGIN = 2, FUN = median, na.rm = TRUE)) }) 
  names(data.use) <- levels(factor(pd[[x]]$cluster_new)) 
  data.frame(as.data.frame(x = t(x = as.data.frame(x = data.use))), cluster = levels(factor(pd[[x]]$cluster_new))) })
#geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
#plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames$x, y = xynames$y, label = id), ...)

# Define subset colors
lapply(names(pd), function(x) table(pd[[x]]$celltype))
pal <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')
blues <- RColorBrewer::brewer.pal(n = 9, name = 'Blues')
reds <- RColorBrewer::brewer.pal(n = 9, name = 'Reds')
cols <- list(subset = list(hd.bd =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40', 'wheat'),#'grey80'),  
                           hd.sk =  c(pal[c(2,6,8,4)], 'darkorchid3'),
                           tp.bd =  c(pal[c(2,6,8,4)], 'darkorchid3', 'grey40'), 
                           tp.sk =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40')), 
             cluster = list(#hd.bd = c(pal[c(1,2)],'royalblue4',pal[c(5)], 'darkorchid3', pal[c(6,8)], 'grey60', '#CD96CD', 'cyan3', 'grey40', pal[c(4)], 'grey80'),
                            #hd.bd = c(pal[c(1,2)],'royalblue4',pal[c(5)], 'darkorchid3', pal[c(6,8)], 'grey60', '#CD96CD', 'cyan3', 'grey40', pal[c(4)], 'grey80'),
                            hd.bd = c(pal[c(1,2)],'royalblue4',pal[c(5)], 'darkorchid3', pal[c(6)],'#CD96CD', 'grey60', pal[c(8)], 'cyan3', 'grey40', pal[c(4)], 'wheat'),#'grey80'),  
                            #hd.sk = c('cyan3',pal[c(1)],'royalblue1',pal[c(2)],pal[c(5)],'royalblue4', 'darkorchid3', pal[c(6)], 'black', '#CD96CD', 'cyan3', 'grey40', pal[c(4)], 'grey80'),
                            hd.sk = c(blues[c(4,3,5,7)],pal[c(5)],blues[c(9)], pal[c(7,8)],'cyan3',pal[c(4)],'darkorchid3',reds[c(8)],pal[c(6)],'#CD96CD', 'grey40', 'grey60', 'grey80'),
                            #hd.sk = c(blues[c(4,3)],'cyan3',blues[c(5)],pal[c(5)],blues[c(7)], pal[c(7,8)],blues[c(9)],pal[c(4)],'darkorchid3',reds[c(8)],pal[c(6)],'#CD96CD', 'grey40', 'grey60', 'grey80'),
                            tp.bd = c(pal[c(1)],pal[c(5)],'grey40','royalblue4',pal[c(6)],'darkorchid3',pal[c(7)],pal[c(2)],pal[c(8)],reds[c(8)],'grey60','grey80',pal[c(4)]), 
                            tp.sk = c(pal[c(7)],blues[c(3,9,5)],pal[c(8)],blues[c(7)],pal[c(5)],pal[c(6)],reds[c(8)],'darkorchid3',pal[c(4)],'#CD96CD','grey40')) )
coln <- list(subset = list(hd.bd =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40', 'wheat'),#'grey80'),  
                           hd.sk =  c(pal[c(2,6,8,4)], 'darkorchid3'),
                           tp.bd =  c(pal[c(2,6,8,4)], 'darkorchid3', 'grey40'), 
                           tp.sk =  c(pal[c(2,6,8,4)], "#CD96CD", 'darkorchid3', 'grey40')), 
             cluster = list(#hd.bd = c('cyan3', pal[c(1,2)],'royalblue4',pal[c(5)], pal[c(6,4,8)], 'darkorchid3', '#CD96CD', 'grey40', 'grey60', 'grey80'),
                            hd.bd = c('cyan3', pal[c(1,2)],'royalblue4',pal[c(5)], pal[c(6,4,8)], '#CD96CD', 'darkorchid3', 'grey40', 'grey60', 'wheat'),#'grey80'),  
                            hd.sk = c(blues[c(3,4,5,7,9)],'cyan3',pal[c(5,6)],reds[c(8)], pal[c(4,7,8)],'darkorchid3'),
                            tp.bd = c(pal[c(2,1)],'royalblue4',pal[c(5)],reds[c(8)],pal[c(6)],pal[c(4,7,8)],'darkorchid3','grey40','grey60','grey80'), 
                            tp.sk = c(blues[c(3,5,7,9)],pal[c(5,6)],reds[c(8)],pal[c(4)],pal[c(7,8)],'#CD96CD','darkorchid3','grey40')) )

# Draw tSNE plots
ts.ss <- Map(x=names(pd), function(x) { 
  ggplot(df[[x]], aes(x = tx, y = ty, color = subset)) + 
    geom_point(shape = 16, size = 0.5, alpha = 1) + #size = 1, alpha = 0.6) +  
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
    scale_color_manual('Subset', values = cols$subset[[x]], guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + #alpha = 0.6))) +  
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample) })
ts.cl <- Map(x=names(pd), function(x) { 
  ggplot(df[[x]], aes(x = tx, y = ty, color = cluster)) + 
    geom_point(shape = 16, size = 0.5, alpha = 1) + 
    geom_point(data = labels.loc[[x]], aes(x = tx, y = ty), size = 6, shape = 16, color = 'white', alpha = 0.7) + 
    geom_text(data = labels.loc[[x]], aes(x = tx, y = ty, label = cluster), size = 4, color = 'black') +
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
    scale_color_manual('Cluster', values = cols$cluster[[x]], guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample) })
ts.cn <- Map(x=names(pd), function(x) { 
  #ggplot(df[[x]], aes(x = tx, y = ty, color = cluster_new)) + 
  ggplot(df[[x]], aes(x = tx, y = ty, color = subset_cluster_new)) + 
    geom_point(shape = 16, size = 0.5, alpha = 1) + 
    geom_point(data = labels.loc_new[[x]], aes(x = tx, y = ty), size = 6, shape = 16, color = 'white', alpha = 0.7) + 
    geom_text(data = labels.loc_new[[x]], aes(x = tx, y = ty, label = cluster), size = 4, color = 'black') + 
    #geom_blank(aes(color = subset), show.legend = T) + facet_wrap(~genotype, ncol = 2) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
    scale_color_manual('Cluster', values = coln$cluster[[x]], guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + #, labels = levels(df[[x]]$subset_cluster)) + 
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample) })
#gridExtra::marrangeGrob(grobs = ts.cl, nrow = 2,ncol=2)



# Define position for labels
pos <- list(hd.bd = data.frame(x=c(-30,25,31,25,-14,20,-15,31),y=c(-15,-37,30,14,40,42,-40,-5),s=levels(df$hd.bd$subset)), 
            #um.cl = data.frame(x=c(24,-28,-24,25,18,16,-36,27),y=c(35,30,-40,-24,-34,-42,16,21),s=levels(df$hd.bd$subset)))
            #hd.bd = data.frame(x=c(-10,-30,25,20,31,25,-14,20,-15,31),y=c(-28,-15,-37,21,30,14,40,42,-40,-5),s=levels(df$hd.bd$subset)), 
            #um.cl = data.frame(x=c(24,-17,-28,-25,-24,25,18,16,-36,27),y=c(35,-4,30,-20,-40,-24,-34,-42,16,21),s=levels(df$hd.bd$subset)))
            hd.sk = data.frame(x=c(19,9,28,-25,25),y=c(-28,30,-15,-20,25),s=levels(df$hd.sk$subset)), 
            #um.cl = data.frame(x=c(6,6,-4,3,3),y=c(-2,5,8,-5,8),s=levels(df$hd.sk$subset)))
            tp.bd = data.frame(x=c(-23,-23,-5,30.5,27,10),y=c(30,-13,-11,2,-35,-37),s=levels(df$tp.bd$subset)), 
            #um.cl = data.frame(x=c(-20,35,23,4,13,-35),y=c(25,32,-40,-20,-55,-20),s=levels(df$tp.bd$subset)))
            #tp.bd = data.frame(x=c(-23,27.5,-23,27,-5,30.5,27,10),y=c(30,23,-13,-17,-11,2,-35,-37),s=levels(df$tp.bd$subset)), 
            #um.cl = data.frame(x=c(-20,34,35,-24,23,4,13,-35),y=c(25,-20,32,-32,-40,-20,-55,-20),s=levels(df$tp.bd$subset)))
            tp.sk = data.frame(x=c(-22,25,-9,2,-9,14,7),y=c(-21,20,30,-25,-23,30,-31),s=levels(df$tp.sk$subset)) )
            #um.cl = data.frame(x=c(-8,-30,23,10,-35,-26,-48),y=c(-28,37,-41,10,-7,-40,0),s=levels(df$tp.sk$subset)))

# Plot data with labels inside
ts.sl <- Map(x=names(pd), function(x) { 
  ggplot(df[[x]], aes(x = tx, y = ty, color = subset)) + 
    geom_point(shape = 16, size = 0.5, alpha = 1) + #size = 1, alpha = 0.6) +  
    annotate("text", x = pos[[x]]$x, y = pos[[x]]$y, size = 5, label = pos[[x]]$s, color = cols$subset[[x]]) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
    scale_color_manual('Subset', values = cols$subset[[x]], guide = guide_legend(override.aes = list(size = 4, alpha = 1))) + #alpha = 0.6))) +  
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample) })

#-----------------------------------------------------------------------------------------------------------------#
# Percentage expressing cells and expression levels per cluster of cell type marker genes
#-----------------------------------------------------------------------------------------------------------------#

# Define cell type markers
gn <- c(`T cell`='CD3E', `CD4 T cell`='CD4', `CD4 T cell`='IL7R', `CD8 T cell`='CD8A', `CD8 T cell`='CD8B', 
        `Treg`='CTLA4', `Treg`='ICOS', `Treg`='FOXP3', #`Treg`='IL2RA', `Treg`='ITGAE', `RORgT cell`='KLRB1', 
        `NK and NKT cell`='KLRB1', `NK cell`='GNLY', `NK cell`='GZMB', #`NK and NKT cell`='GZMK',#`NK cell`='NKG7', `NK cell`='GZMA', 
        `Tn`='CCR7', `Tn`='SELL', `Tm`='S100A4', #`Tn`='PECAM1', `Tm`='FAS', `Tem`='B3GAT1', `HSC`='CD34', 
        #`Tn and Tcm`='CD27', `Tn and Tcm`='CD28', `Trm`='CD69', `Trm`='PDCD1', `Trm`='CXCR6', `Trm`='ITGA1', 
        `B cell`='MS4A1', #`B cell`='CD19', `B cell`='CXCR5', `B cell`='CCR6', `DC`='FCER1A', `DC`='CST3', 
        #`Monocyte`='FCGR3A', `Monocyte`='CD14', `Monocyte`='LYZ', `Monocyte`='PTPRC', `Monocyte`='CCR2', `Monocyte`='MS4A7', 
        #`trafficking`='CX3CR1', `trafficking`='ITGB1', `trafficking`='SELPLG', 
        `gdT cell`='TRDC', `gdT cell`='TRGC1', `gdT cell`='TRGC2', `abT cell`='TRAC', `abT cell`='TRBC1', `abT cell`='TRBC2')

# Extract expression data
df.mk <- Map(x=names(pd), function(x) { 
  df.int <- data.frame(cells = colnames(pd[[x]]), 
                       tx = pd[[x]]@reductions$tsne@cell.embeddings[,1], 
                       ty = pd[[x]]@reductions$tsne@cell.embeddings[,2], 
                       ux = pd[[x]]@reductions$umap@cell.embeddings[,1], 
                       uy = pd[[x]]@reductions$umap@cell.embeddings[,2], 
                       sample = pd[[x]]@meta.data$sample, 
                       cluster = pd[[x]]@meta.data$cluster,
                       #cluster_new = factor(pd[[x]]$cluster_new, levels = sort(as.numeric(levels(factor(pd[[x]]$cluster_new))))),
                       cluster_new = pd[[x]]@meta.data$cluster_new,
                       compartment = pd[[x]]@meta.data$subset, 
                       subset = pd[[x]]@meta.data$celltype, 
                       #subset_cluster = paste0(pd[[x]]@meta.data$cluster_new,' - ',pd[[x]]@meta.data$celltype),
                       subset_cluster = pd[[x]]@meta.data$subset_cluster,
                       subset_cluster_new = pd[[x]]@meta.data$subset_cluster_new, 
                       genotype = pd[[x]]@meta.data$genotype, 
                       sex = ifelse(as.matrix(GetAssayData(pd[[x]])["RPS4Y1",])>0, 'male', 'n.d.'), 
                       t(as.matrix(GetAssayData(pd[[x]], slot = "data")[gn,])) )
  reshape2::melt(df.int, id.vars = names(df.int)[1:14]) })
lapply(df.mk, dim)

# Calculate percentage of expressing cells for each gene
#sort(unique(names(gn), decreasing = F))[c(8,1,5,10,9,3,4,11,7,6,2)]
sort(unique(names(gn), decreasing = F))[c(8,1,10,9,3,4,11,5,6,7,2)]
df.pct <- Map(x=names(pd), function(x) df.mk[[x]][,c('cells','subset_cluster_new','genotype','variable','value')] )
df.pct.cl <- Map(x=names(pd), function(x) { 
  df.int <- plyr::ddply(df.pct[[x]], c('subset_cluster_new','variable'), summarise, 
                        percentage = 100*(sum(value>0)/length(value)), value = mean(value))
  #df.int$variable <- as.character(df.int$variable)
  df.int$variable <- factor(df.int$variable, levels = sort(levels(df.int$variable), decreasing = T))
  df.int$celltype <- names(gn[match(df.int$variable, gn)])
  df.int$celltype <- factor(df.int$celltype, levels = sort(unique(df.int$celltype), decreasing = F)[c(8,1,10,9,3,4,11,5,6,7,2)])
  return(df.int) })

# Plot data
#col.hm <- viridis::viridis_pal()(52)[1:50]
col.hm <- colorRampPalette(c('grey80', 'skyblue2', 'royalblue1', 'royalblue4'))(50)
dot.pct.cl <- Map(x=names(pd), function(x) 
  ggplot(df.pct.cl[[x]], aes(subset_cluster_new, variable, color = value, size = percentage)) + 
    geom_point(shape = 16) + 
    scale_x_discrete('clusters', expand = c(0,0.6), position = 'top') + 
    scale_y_discrete('', expand = c(0,0.6)) + 
    scale_color_gradientn('Expression', colours = col.hm, oob = scales::squish, #expression(''*log[2]*'(exp)'), 
                          limits = c(0,2), breaks = seq(0,2, length.out = 3)) + 
    scale_size_continuous('Percentage', limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
    ggtitle('HSCT skin') + theme_custom + facet_grid(celltype~., scales = 'free', space = 'free') + 
    theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
          plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
          axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 0), #size = 12.8, 
          strip.text.y = element_text(size = 12, face = 'bold', margin = margin(0,0,0,5), angle=0,vjust=0.5,hjust=0)) )

#-----------------------------------------------------------------------------------------------------------------#
# Percentage expressing cells and expression levels per cluster (remove memory and naive markers)
#-----------------------------------------------------------------------------------------------------------------#

# Filter data table
df.pct.clf <- Map(x=names(pd), function(x) df.pct.cl[[x]][df.pct.cl[[x]]$celltype %in% levels(df.pct.cl[[x]]$celltype)[-c(3,4)],] )

# Plot data
dot.pct.clf <- Map(x=names(pd), function(x) 
  ggplot(df.pct.clf[[x]], aes(subset_cluster_new, variable, color = value, size = percentage)) + 
    geom_point(shape = 16) + 
    scale_x_discrete('clusters', expand = c(0,0.6), position = 'top') + 
    scale_y_discrete('', expand = c(0,0.6)) + 
    scale_color_gradientn('Expression', colours = col.hm, oob = scales::squish, #expression(''*log[2]*'(exp)'), 
                          limits = c(0,2), breaks = seq(0,2, length.out = 3)) + 
    scale_size_continuous('Percentage', limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
    ggtitle('Healthy skin') + theme_custom + facet_grid(celltype~., scales = 'free', space = 'free') + 
    theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
          plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
          axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 0), #size = 12.8, 
          strip.text.y = element_text(size = 12, face = 'bold', margin = margin(0,0,0,5), angle=0,vjust=0.5,hjust=0)) )

#-----------------------------------------------------------------------------------------------------------------#
# Percentage expressing cells and expression levels per cluster and genotype
#-----------------------------------------------------------------------------------------------------------------#

## Prepare data table
#df.pct.gt <- plyr::ddply(df.pct, c('cluster','genotype','variable'), summarise, 
#                         percentage = 100*(sum(value>0)/length(value)), value = mean(value))
#df.pct.gt$variable <- as.character(df.pct.gt$variable)
#df.pct.gt$celltype <- names(gn[match(df.pct.gt$variable, gn)])
#df.pct.gt$variable <- factor(df.pct.gt$variable, levels = sort(unique(df.pct.gt$variable), decreasing = T))
#df.pct.gt$celltype <- factor(df.pct.gt$celltype, levels = sort(unique(df.pct.gt$celltype), decreasing = F)[c(6,1,4,2,3,7,5)])
#
## Plot data
#dot.pct.gt <- ggplot(df.pct.gt, aes(cluster, variable, color = value, size = percentage)) + 
#  geom_point(shape = 16) + 
#  scale_x_discrete('clusters', expand = c(0,0.6), position = 'top') + 
#  scale_y_discrete('', expand = c(0,0.6)) + 
#  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col.hm, oob = scales::squish,
#                        limits = c(0,2), breaks = seq(0,2, length.out = 3)) + 
#  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
#  ggtitle('HSCT skin') + theme_custom + facet_grid(celltype~genotype, scales = 'free_y', space = 'free', switch = 'x') + 
#  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
#        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
#        strip.text.y = element_text(size = 12, face = 'bold', margin = margin(0,0,0,5), angle=0,vjust=0.5,hjust=0), 
#        strip.text.x = element_text(size = 12, face = 'bold', margin = margin(5,0,0,0), angle=0,vjust=0.5,hjust=0.5))

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Funciton to save plots with same size
grid <- function(data, layout=rbind(c(1,1,1,1,2,3,3,3))) {
  grid.arrange(grobs = list(data + guides(color = 'none', fill = 'none', size = 'none', shape = 'none'), 
                            cowplot::get_legend(data)), ncol = 2, layout_matrix = layout, 
               top = grid::textGrob('', gp = grid::gpar(fontsize = 16, font = 2)))
}

# Save plots for subset annotation
#lapply(names(pd), function(x) ggsave(plot = grid(ts.ss[[x]]), file = paste0(dir.figs[[x]], "_tsne_annotation.pdf"), height = 4, width = 8) )
#lapply(names(pd), function(x) ggsave(plot = grid(ts.cl[[x]]), file = paste0(dir.figs[[x]], "_tsne_clusters.pdf"), height = 4, width = 8) )
#lapply(names(pd), function(x) ggsave(plot = grid(ts.cn[[x]]), file = paste0(dir.figs[[x]], "_tsne_clusters_new.pdf"), height = 4, width = 8) )
lapply(names(pd), function(x) ggsave(plot = ts.ss[[x]], file = paste0(dir.figs[[x]], "_tsne_annotation.pdf"), height = 3.8, width = 6) )
lapply(names(pd), function(x) ggsave(plot = ts.cl[[x]], file = paste0(dir.figs[[x]], "_tsne_clusters.pdf"), height = 3.8, width = 6) )
lapply(names(pd), function(x) ggsave(plot = ts.cn[[x]], file = paste0(dir.figs[[x]], "_tsne_clusters_new.pdf"), height = 3.8, width = 6) )
lapply(names(pd), function(x) ggsave(plot = ts.sl[[x]], file = paste0(dir.figs[[x]], "_tsne_annotation_label.pdf"), height = 3.8, width = 6) )

# Draw plots for saving
sv.list <- list(ts.mn = grid(ts.mn), um.mn = grid(um.mn), ts.mnl = grid(ts.mnl), um.mnl = grid(um.mnl),  
                ts.mn.gt = grid(ts.mn.gt, rbind(c(1,1,1,2))), um.mn.gt = grid(um.mn.gt, rbind(c(1,1,1,2))), 
                ts.tcr = grid(ts.tcr, rbind(c(1,1,1,1,2))), ts.cd8 = grid(ts.cd8, rbind(c(1,1,1,1,2))),
                um.tcr = grid(um.tcr, rbind(c(1,1,1,1,2))), um.cd8 = grid(um.cd8, rbind(c(1,1,1,1,2))), 
                dot.pct.cl = grid(dot.pct.cl, rbind(c(1,1,2))), dot.pct.clf = grid(dot.pct.clf, rbind(c(1,1,2))), dot.pct.gt = grid(dot.pct.gt, rbind(c(1,1,2))),  
                dot.pct.tcr.cl = grid(dot.pct.tcr.cl, rbind(c(1,2))), dot.pct.tcr.gt = grid(dot.pct.tcr.gt, rbind(c(1,1,2))))

# Save plots
lapply(names(pd), function(x) ggsave(plot = dot.pct.cl[[x]], file = paste0(dir.figs[[x]], '_pct_expression_clusters.pdf'), height = 8.5, width = 7.5) )
lapply(names(pd), function(x) ggsave(plot = dot.pct.clf[[x]], file = paste0(dir.figs[[x]], '_pct_expression_clusters_filtered.pdf'), height = 7.7, width = 7.5) )
#ggsave(plot = sv.list$dot.pct.gt, file = paste0(dir.figs.sk, '_pct_expression_genotype.pdf'), height = 6.8, width = 11)

#-----------------------------------------------------------------------------------------------------------------#
# End of script
#-----------------------------------------------------------------------------------------------------------------#

#=================================================================================================================#
# R session info
#=================================================================================================================#

#> devtools::session_info()
#─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#setting  value                       
#version  R version 3.5.1 (2018-07-02)
#os       macOS  10.14.6              
#system   x86_64, darwin15.6.0        
#ui       RStudio                     
#language (EN)                        
#collate  en_US.UTF-8                 
#ctype    en_US.UTF-8                 
#tz       Europe/Berlin               
#date     2020-03-11                  
#
#─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
#pheatmap      * 1.0.12     2019-01-04 [1] CRAN (R 3.5.2)                   
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
#RColorBrewer  * 1.1-2      2014-12-07 [1] CRAN (R 3.5.0)                   
#Rcpp            1.0.0      2018-11-07 [1] CRAN (R 3.5.0)                   
#RCurl           1.95-4.11  2018-07-15 [1] CRAN (R 3.5.0)                   
#Rdpack          0.10-1     2018-10-04 [1] CRAN (R 3.5.0)                   
#remotes         2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
#reshape2        1.4.3      2017-12-11 [1] CRAN (R 3.5.0)                   
#reticulate      1.10       2018-08-05 [1] CRAN (R 3.5.0)                   
#rlang           0.4.0      2019-06-25 [1] CRAN (R 3.5.2)                   
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
#usethis         1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
#viridis         0.5.1      2018-03-29 [1] CRAN (R 3.5.0)                   
#viridisLite     0.3.0      2018-02-01 [1] CRAN (R 3.5.0)                   
#withr           2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
#XML             3.98-1.16  2018-08-19 [1] CRAN (R 3.5.0)                   
#xtable          1.8-3      2018-08-29 [1] CRAN (R 3.5.0)                   
#yaml            2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
#zip             1.0.0      2017-04-25 [1] CRAN (R 3.5.0)                   
#zoo             1.8-4      2018-09-19 [1] CRAN (R 3.5.0)                   
#
#[1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library