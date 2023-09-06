#=================================================================================================================#
# Exploratory analysis of BD Rhapsody targeted scRNA-Seq from CD3+ T cells from healthy skin and blood (matched)
# AbSeq: CD45RA, CD45RO, CD4, CD103, CD27, CD28, CD95, CCR7, PD-1, CD69
# Date: 2020.07.09
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze BD Rhapsody data processed with Seven Bridges Genomics platform
# - gene set expression scores for Tscm signatures (overlaping differentially expressed genes when compared to 
# Tn, Tem and Tcm from Gatinoni PMID21926977)
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
.libPaths(paste0(.libPaths(), '/Seurat_v3.1.1'))#, paste0(.libPaths(), '/Seurat_v3.1.1'))
library(openxlsx)
library(Seurat, lib.loc = paste0(.libPaths(), '/Seurat_v3.1.1'))
library(dplyr)
library(SingleR)
library(ggplot2)
library(gridExtra)
library(destiny)
library(Biobase)
library(ComplexHeatmap)
library(GSEABase)
#library(pheatmap)
#library(stringr)
#library(reshape2)
#library(plyr)
#library(gridExtra)
#library(scales)
#library(VennDiagram)
#library(grid)
#library(readr)
#library(vsn)
#library(edgeR)
#library(biomaRt)
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
                      strip.text = element_text(size = 16, face = "plain", margin = margin(0,0,0,5)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid <- c("EX0010")

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
an.desc <- c("sc_bd_hd_skin_vs_blood_gsea_tscm_gatinoni")
an.descs <- c("sc_hd_skin_vs_blood_tscm")

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
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/tissue"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/tissue_genotype"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/tissue_subset"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/tissue_subset_genotype"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood_memory_genotype"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/tsne"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/CD4"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/CD8"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/seurat"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/dpt"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/dpt_not_integrated"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
dir.figs <- list()
dir.figs$all <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
#dir.figs$t <- paste0("figures/", analysisid, "_", an.desc, "/tissue/", analysisid, "_", an.descs, '_tissue')
#dir.figs$tg <- paste0("figures/", analysisid, "_", an.desc, "/tissue_genotype/", analysisid, "_", an.descs, '_tissue_genotype')
#dir.figs$ts <- paste0("figures/", analysisid, "_", an.desc, "/tissue_subset/", analysisid, "_", an.descs, '_tissue_subset')
#dir.figs$tsg <- paste0("figures/", analysisid, "_", an.desc, "/tissue_subset_genotype/", analysisid, "_", an.descs, '_tissue_subset_genotype')
#dir.figs$bmg <- paste0("figures/", analysisid, "_", an.desc, "/blood_memory_genotype/", analysisid, "_", an.descs, '_blood_memory_genotype')
dir.figs$tsne <- paste0("figures/", analysisid, "_", an.desc, "/tsne/", analysisid, "_", an.descs, '_tsne')
#dir.figs$CD4 <- paste0("figures/", analysisid, "_", an.desc, "/CD4/", analysisid, "_", an.descs)
#dir.figs$CD8 <- paste0("figures/", analysisid, "_", an.desc, "/CD8/", analysisid, "_", an.descs)
#dir.figss <- paste0("figures/", analysisid, "_", an.desc, "/seurat/", analysisid, "_", an.descs)
#dir.figs.dpt <- paste0("figures/", analysisid, "_", an.desc, "/dpt/", analysisid, "_", an.descs)
#dir.figs.dpt.ni <- paste0("figures/", analysisid, "_", an.desc, "/dpt_not_integrated/", analysisid, "_", an.descs)
#dir.figs.bd <- paste0("figures/", analysisid, "_", an.desc, "/blood/", analysisid, "_", an.descs, "_blood")
#dir.figs.sk <- paste0("figures/", analysisid, "_", an.desc, "/skin/", analysisid, "_", an.descs, "_skin")

#=================================================================================================================#
# Analysis of processed data
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Load processed data
#-----------------------------------------------------------------------------------------------------------------#

# Load seurat object
bd <- readRDS(paste0('../201911/figures/GA_AN0296_sc_bd_blood_skin_tcell/GA_AN0296_sc_tcell_seurat.rds'))
#bd.markers <- read.csv(file = paste0('../201911/figures/GA_AN0296_sc_bd_blood_skin_tcell/GA_AN0296_sc_tcell_cluster_markers.csv'), row.names = 1)

# # Add sample name to cluster names
# #grid.arrange(DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'RNA_snn_res.1.9', split.by = 'orig.ident', label = T, label.size = 4, pt.size = 0.5, cols = NULL) + theme_bw() + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)), ts.ct, nrow = 1)
# #dupaxis.labels <- df.sm[!duplicated(df.sm$cluster),]$sample
# sample.cl <- c('blood_stim','blood_unstim','skin_unstim','skin_stim','blood_unstim','skin','stim','skin_stim',
#                'blood_stim','blood_stim','skin_unstim','unstim','skin_unstim','skin_stim','skin_unstim','skin_unstim',
#                'skin_unstim','skin_unstim','blood_unstim','blood_stim')
# bd.markers$cluster.name <- plyr::mapvalues(bd.markers$cluster, from = levels(factor(bd.markers$cluster)),
#                                            to = paste0(sample.cl,'_',levels(factor(bd.markers$cluster))))

# Load barcodes for each gate
path <- paste0('../../procdata/EX0010/SM_PD0010_flowjow_gating/Rhapsody_blood_skin_I/DBEC/FlowJo Gating/barcodes/')
file <- list.files(path, pattern = '.csv')[-c(5,6,11,12)]
fp <- file.path(path,file)
names(fp) <- gsub(' ','_', stringr::str_split(file, "\\.", n = 2, simplify = T)[, 1])
bc <- lapply(fp, function(x) read.csv(x, stringsAsFactors = F, col.names = 'barcode'))
bc <- bind_rows(bc, .id = 'subset')
bc$subset[duplicated(bc$barcode)]
table(bc$subset)

# Transfer new labels to seurat object
# #map <- data.frame(old = sort(as.numeric(names(Idents(bd)))), stringsAsFactors = F)
# map <- data.frame(old = (as.numeric(names(Idents(bd)))), stringsAsFactors = F)
# map$old <- as.numeric(map$old)
# map$new <- as.numeric(1:length(map$old))
# #map <- map[order(map$bc),]
# bd$bc <- as.numeric(names(Idents(bd)))
# bd$bc <- map$new[match(bd$bc, map$old)]
map <- read.csv('../../procdata/EX0010/SM_PD0010_flowjow_gating/Rhapsody_blood_skin_I/DBEC/FlowJo Gating/FLOWJO_mapping.csv', row.names = 1)
bd$bc <- as.numeric(names(Idents(bd)))
bd$bc <- map$CELLBARCODE_FLOWJO[match(bd$bc, map$Cell_Index)]

# Annotate cells according to flowjow gating
bd$subset <- bc$subset[match(bd$bc, bc$barcode)]
bd$subset[is.na(bd$subset)] <- 'n.d.'
table(bd$subset)
#saveRDS(bd, file = paste0(dir.figs, '_seurat_annotated.rds'))
#bd <- readRDS(file = paste0(dir.figs, '_seurat_annotated.rds'))
bd$tissue <- stringr::str_split(bd$Sample_Name, '_', n = 2, simplify = T)[, 1]
bd$status <- stringr::str_split(bd$Sample_Name, '_', n = 2, simplify = T)[, 2]
bd$subset <- plyr::mapvalues(bd$subset, from = levels(factor(bd$subset)), 
                             to = c('CD4- Tcm','CD4- Tem','CD4- Tn','CD4- Temra','CD4+ Tcm', 'CD4+ Tem', 'CD4+ Tn', 'CD4+ Temra', 'n.d.'))

#-----------------------------------------------------------------------------------------------------------------#
# Define gene sets
#-----------------------------------------------------------------------------------------------------------------#

## Prepare gene set
#molsig <- getGmt("../../databases/msigdb/msigdb.v7.0.symbols.gmt")
#molsigf <- molsig[grep("METABOLI", names(molsig)), ]
#molsigf <- molsigf[grep("KEGG", names(molsigf)), ]
#length(unique(names(molsigf)))
#geneSets <- subsetGeneSets(molsigf, rownames(ct$tp.sk))
#cbind(nGenes(geneSets))
##Add set size to its name
#geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), "_", nGenes(geneSets) ,"g", sep=""))

## Load gene sets from MSigDB
#molsig <- GSEABase::getBroadSets("../../databases/msigdb/msigdb_v7.0.xml")
#molsigf <- setNames(molsig@.Data, nm=names(molsig))
#molsigf <- molsigf[sapply(molsigf, function(x) bcSubCategory(collectionType(x))) %in% c('CP:KEGG')]
##molsigf <- molsigf[grep('RESTING',sapply(molsigf, function(x) x@setName))]
#molsigf <- molsigf[grep('GLYCOLYSIS|TCA|OXIDATIVE|PENTOSE_PHOS', names(molsigf))]
##names(molsigf) <- c('TCR_Signaling')

# Load table
#tscm.sig <- read.csv(paste0('../../databases/celltype/tscm/suptab5_tscm_vs_tcm_diff_genes.csv'), row.names = NULL)[,1:7]
#tscm.sig[tscm.sig$Tscm_vs_Tn_pval < 0.05 & tscm.sig$Tscm_vs_Tcm_pval < 0.05 & tscm.sig$Tscm_vs_Tem_pval < 0.05, ]
#tscm.sig[tscm.sig$Tscm_vs_Tn_pval < 0.05 & tscm.sig$Tscm_vs_Tcm_pval < 0.05 & tscm.sig$Tscm_vs_Tem_pval < 0.05 & 
#           tscm.sig$Tscm_vs_Tn_fc > 0 & tscm.sig$Tscm_vs_Tcm_fc > 0 & tscm.sig$Tscm_vs_Tem_fc > 0, ]
#tscm.sig[tscm.sig$Tscm_vs_Tcm_pval < 0.05 & tscm.sig$Tscm_vs_Tem_pval < 0.05, ]
df <- read.csv(paste0('../../databases/celltype/tscm/suptab2_cd8t_subsets_diff_genes.csv'), row.names = NULL, stringsAsFactors = F)
df$gene[!(df$gene %in% rownames(bd))]
df <- df[!(df$gene %in% 'UNKNOWN'),]

# Change aliases to symbols
genemap <- read.csv(paste0('../../databases/hgnc/20200709/biomart_genes.csv'), stringsAsFactors = F)
genemap[grep("PRAGMIN",genemap$Alias.symbol),]
genemap <- genemap[!duplicated(genemap$Alias.symbol),]
#genemap <- read.csv(paste0('../../databases/biomart/genemap_10x_GA_AN0260.csv'))
#idx <- match(df$gene, genemap$external_gene_name)
idx <- match(df$gene, genemap$Alias.symbol)
df$gene[!(is.na(idx))] <- genemap$Approved.symbol[idx[!(is.na(idx))]]
df$gene[!(df$gene %in% rownames(bd))]

# Define gene sets
fc.cutoff <- 0
gs <- list()
gs$Tscm_vs_Tn <- df[df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc > fc.cutoff, ]$gene
gs$Tscm_vs_Tcm <- df[df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tcm_fc > fc.cutoff, ]$gene
gs$Tscm_vs_Tem <- df[df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tem_fc > fc.cutoff, ]$gene
gs$Tn_vs_Tscm <- df[df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc < -fc.cutoff, ]$gene
gs$Tcm_vs_Tscm <- df[df$Tscm_vs_Tcm_pval < 0.05 & df$Tscm_vs_Tcm_fc < -fc.cutoff, ]$gene
gs$Tem_vs_Tscm <- df[df$Tscm_vs_Tem_pval < 0.05 & df$Tscm_vs_Tem_fc < -fc.cutoff, ]$gene
#gs$Tscm <- df[(df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc > fc.cutoff) & 
#                 (df$Tscm_vs_Tcm_pval < 0.05 & df$Tscm_vs_Tcm_fc > fc.cutoff) & 
#                 (df$Tscm_vs_Tem_pval < 0.05 & df$Tscm_vs_Tem_fc > fc.cutoff), ]$gene
#gs$Tscm_dn <- df[(df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc < -fc.cutoff) & 
#                 (df$Tscm_vs_Tcm_pval < 0.05 & df$Tscm_vs_Tcm_fc < -fc.cutoff) & 
#                 (df$Tscm_vs_Tem_pval < 0.05 & df$Tscm_vs_Tem_fc < -fc.cutoff), ]$gene
#FeaturePlot(pd$tp.bd, reduction = 'tsne', features = c('TDRD6','HAVCR1'))
#FeaturePlot(pd$tp.bd, reduction = 'tsne', features = c('TDRD6','HAVCR1'))
#df[(df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc < 0) & (df$Tscm_vs_Tcm_pval < 0.05 & df$Tscm_vs_Tcm_fc < 0) & (df$Tscm_vs_Tem_pval < 0.05 & df$Tscm_vs_Tem_fc < 0), ]$gene
#df[(df$Tscm_vs_Tn_pval < 0.05 & df$Tscm_vs_Tn_fc > 0) & (df$Tscm_vs_Tcm_pval < 0.05 & df$Tscm_vs_Tcm_fc > 0) & (df$Tscm_vs_Tem_pval < 0.05 & df$Tscm_vs_Tem_fc > 0), ]$gene
sapply(gs, length)
sapply(names(gs), function(x) sum(gs[[x]] %in% rownames(bd) / length(gs[[x]]) ))
gs <- Map(x=names(gs), function(x) gs[[x]][!duplicated(gs[[x]])])
#gs <- Map(x=names(gs), function(x) gs[[x]][gs[[x]] %in% rownames(bd)])
gs <- Map(x=names(gs), function(x) gs[[x]][gs[[x]] %in% rownames(bd)[rowSums(as.matrix(GetAssayData(bd, slot = 'data', assay = 'RNA'))) > 0]])
sapply(gs, length)

# Save list of genes from gene sets
gss <- as.data.frame(sapply(gs, "[", i = seq_len(max(sapply(gs, length)))), stringsAsFactors=F )
gss[is.na(gss)] <- ''
write.csv(gss, file = paste0(dir.figs$all, '_geneset.csv'), row.names = T)

# Extract relevant gene set
#gs <- subsetGeneSets(molsigf, rownames(ct))
#gs <- Map(x=sort(names(molsigf)), function(x) molsigf[[x]]@geneIds[molsigf[[x]]@geneIds %in% rownames(bd)])

#=================================================================================================================#
# Comparing gene set scores between tissues
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Calculate module scores for gene sets using Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Create color palette
col.hm <- colorRampPalette(c('grey80', 'skyblue2', 'royalblue1', 'royalblue4'))(50)
col.cl <- c('white','white','grey60','grey60')#c(blues[c(1,2,3,4)],'firebrick3')

# Define cells to include in the analysis
bd[['include_tissue']] <- bd$orig.ident %in% levels(factor(bd$orig.ident))
#bd <- subset(bd, subset = include_tissue == T)
#bd$type <- ifelse(stringr::str_split(bd$subset, ' ', n = 2, simplify = T)[,1]=='CD4','CD4','CD8')

# Split data by tissue
#tp <- list(bd = subset(bd, subset = cluster_split %in% levels(bd$cluster_split)[c(6)]), 
#           sk = subset(bd, subset = cluster_split %in% levels(bd$cluster_split)[c(1:4)]) )
#tp <- list(bd = subset(bd, subset = orig.ident %in% levels(factor(bd$orig.ident))[c(1)]), 
#           sk = subset(bd, subset = orig.ident %in% levels(factor(bd$orig.ident))[c(2)]) )
tp <- list(tp = bd)

# Create function to rescale scores for al sets together (the absolute values from 0 to 1 and then keep the negative sign for negative values)
rescale <- function(scores,to=c(0,1)) { scales::rescale(abs(scores),to=to) * ifelse(scores>=0,1,-1) }
rescale.all <- function(data,to=c(0,1)) { 
  df <- reshape2::melt(as.matrix(data))
  df$value <- scales::rescale(abs(df$value),to=to) * ifelse(df$value>=0,1,-1)
  df <- reshape2::dcast(df, Var1 ~ Var2)
  transform(df, row.names = Var1, Var1 = NULL)
}

# Calculate module scores (remove genes not detected in any cell for binning to work)
tp <- Map(x=names(tp), function(x) { #lapply(setNames(names(gs), nm = names(gs)) , function(y) {
  all <- AddModuleScore(object = subset(tp[[x]], features = rownames(tp[[x]])[
    rowSums(as.matrix(GetAssayData(tp[[x]], slot = 'data', assay = 'RNA'))) > 0]), features = gs, ctrl = 10, nbin = 24, seed = 1, 
    name = 'scores', assay = 'RNA')
  AddMetaData(object = tp[[x]], metadata = all@meta.data[,grep('scores',colnames(all@meta.data))], col.name = names(gs)) })
tp$tp@meta.data$Tscm_vs_Tn_all <- tp$tp@meta.data$Tscm_vs_Tn - tp$tp@meta.data$Tn_vs_Tscm
tp$tp@meta.data$Tscm_vs_Tcm_all <- tp$tp@meta.data$Tscm_vs_Tcm - tp$tp@meta.data$Tcm_vs_Tscm
tp$tp@meta.data$Tscm_vs_Tem_all <- tp$tp@meta.data$Tscm_vs_Tem - tp$tp@meta.data$Tem_vs_Tscm
tp$tp@meta.data$Tscm_all <- tp$tp@meta.data$Tscm_vs_Tcm_all + tp$tp@meta.data$Tscm_vs_Tem_all + tp$tp@meta.data$Tscm_vs_Tn_all

# Rescale all scores
#tp$bd@meta.data[,15:18] <- apply(tp$bd@meta.data[,15:18], 2, rescale, to=c(0,1))
#tp$sk@meta.data[,15:18] <- apply(tp$sk@meta.data[,15:18], 2, rescale, to=c(0,1))
#scaled_scores <- rescale.all(rbind(tp$bd@meta.data[,17,drop=F],tp$sk@meta.data[,17,drop=F]), to=c(0,1))
#scaled_scores <- apply(tp$tp@meta.data[,grep('KEGG',colnames(tp$tp@meta.data)),drop=F], 2, rescale.all, to=c(0,1))
scaled_scores <- lapply(setNames(colnames(tp$tp@meta.data)[grep('Tscm',colnames(tp$tp@meta.data))], 
                                 nm = colnames(tp$tp@meta.data)[grep('Tscm',colnames(tp$tp@meta.data))]), 
                        function(x) rescale.all(tp$tp@meta.data[,x,drop=F], to=c(0,1)) )
#tp$bd@meta.data[,17] <- scaled_scores[grep('bd', rownames(scaled_scores)),]
#tp$sk@meta.data[,17] <- scaled_scores[grep('sk', rownames(scaled_scores)),]
tp$tp@meta.data[,grep('Tscm',colnames(tp$tp@meta.data))] <- scaled_scores

#-----------------------------------------------------------------------------------------------------------------#
# Plot results in tSNE
#-----------------------------------------------------------------------------------------------------------------#

# Plot scores for lymphocytes from blood
gs.names <- setNames(colnames(tp$tp@meta.data)[grep('Tscm',colnames(tp$tp@meta.data))], 
                     nm = colnames(tp$tp@meta.data)[grep('Tscm',colnames(tp$tp@meta.data))])
df.gtm <- df.gta <- histplotta <- boxplotta <- list()
for(donor in sort(levels(bd$orig.ident))) {
  df.gta[[donor]] <- lapply(gs.names, function(y) {
    df.int2 <- Map(x=names(tp), function(x) { 
      df.int <- data.frame(tx = bd@reductions$tsne@cell.embeddings[,1], ty = bd@reductions$tsne@cell.embeddings[,2], 
                           #tx = c(as.vector(pd$tp.bd@reductions$tsne@cell.embeddings[,1]), as.vector(pd$tp.sk@reductions$tsne@cell.embeddings[,1])),
                           #ty = c(as.vector(pd$tp.bd@reductions$tsne@cell.embeddings[,2]), as.vector(pd$tp.sk@reductions$tsne@cell.embeddings[,2])),
                           #tx = bd@reductions$tsne@cell.embeddings[,1], ty = bd@reductions$tsne@cell.embeddings[,2],
                           #tx = pd[[ifelse(x=='blood','tp.bd','tp.sk')]]@reductions$tsne@cell.embeddings[,1], 
                           #ty = pd[[ifelse(x=='blood','tp.bd','tp.sk')]]@reductions$tsne@cell.embeddings[,2], 
                           subset = tp[[x]]@meta.data$subset, cell = rownames(tp[[x]]@meta.data), donor = tp[[x]]@meta.data$orig.ident, status = tp[[x]]@meta.data$status, 
                           score = tp[[x]]@meta.data[[y]][,1], tissue = tp[[x]]@meta.data$tissue, stringsAsFactors = T)[tp[[x]]@meta.data$include_tissue,]
      df.int <- df.int[df.int$donor == donor,]
      #df.int$score <- rescale(df.int$score,to=c(0,1))
      return(df.int) })
    #score = scales::rescale(tp[[x]]@meta.data[[geneSetName]],to=c(0,1)))[tp[[x]]$include_tissue,] })
    df.int2 <- bind_rows(df.int2)
    df.int2$tissue <- factor(df.int2$tissue, levels = levels(factor(df.int2$tissue))[c(1,2)])
    return(df.int2) })
  df.gtm[[donor]] <- bind_rows(df.gta[[donor]], .id='geneSetName')
}
df.gtm <- bind_rows(df.gtm, .id='group')
#df.gtm$geneSetName <- gsub('\\.',' ', gsub('_',' ', df.gtm$geneSetName))
#df.gtm$geneSetName <- factor(df.gtm$geneSetName, levels=levels(factor(df.gtm$geneSetName))[c(4,3,2,1,8,7,6,5)])

# Plot tSNE
col.hm2 <- grDevices::colorRampPalette(c(rev(RColorBrewer::brewer.pal(n = 9, name = 'YlGnBu')), (RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[2:9])))(50)
ts <- Map(x=setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(x) { 
    #df.ts <- as.data.frame(pd[[ifelse(x=='blood','tp.bd','tp.sk')]]@reductions$tsne@cell.embeddings)
    #df.ts <- df.ts[match(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue==x,]$cell, paste0(ifelse(x=='blood','bd_','sk_'),rownames(df.ts))),]
    #df.ts$score <- df.gtm[df.gtm$geneSetName==y & df.gtm$tissue == x,]$score
    #df.ts$cell <- df.gtm[df.gtm$geneSetName==y & df.gtm$tissue == x,]$cell
    #n <- 50*y
    #col.neg <- grDevices::colorRampPalette(c(rev(RColorBrewer::brewer.pal(n = 9, name = 'YlGnBu'))[1:8]))(n)
    #col.pos <- grDevices::colorRampPalette(c(RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')[3:9]))(50-n)
    ggplot(df.gtm[df.gtm$geneSetName == x,], aes_string(x='tx', y='ty', color='score')) + 
      geom_point(shape = 16, size = ifelse(x=='blood',0.8,1), alpha = 1) + 
      scale_color_gradientn('Score', colors = col.hm2, oob = scales::squish, limits = c(-1,1), breaks = seq(-1,1, length.out = 3), 
                            labels = scales::number_format(accuracy = 1)) + 
      scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(tissue~status, ncol = 2) + 
      facet_grid(status~tissue, scales='fixed', space='fixed') + 
      theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + theme(strip.text = element_text(size = 16, face = "plain", margin = margin(5,5,5,5))) + 
      ggtitle(paste0('CD3+ lymphocytes'), subtitle = paste0('', gsub('_',' ',x))) })

# Plot tSNE indicating cells with high scores for the comparison between Tscm and the other 3 compartments
tsp <- Map(x=names(tp), function(x) { 
  df.ts <- as.data.frame(bd@reductions$tsne@cell.embeddings)
  df.ts <- df.ts[match(df.gtm[df.gtm$geneSetName==unique(df.gtm$geneSetName)[7],]$cell, rownames(df.ts)),]
  df.ts[,unique(df.gtm$geneSetName)[7]] <- df.gtm[df.gtm$geneSetName == unique(df.gtm$geneSetName)[7],]$score
  df.ts[,unique(df.gtm$geneSetName)[8]] <- df.gtm[df.gtm$geneSetName == unique(df.gtm$geneSetName)[8],]$score
  df.ts[,unique(df.gtm$geneSetName)[9]] <- df.gtm[df.gtm$geneSetName == unique(df.gtm$geneSetName)[9],]$score
  df.ts$tissue <- df.gtm[df.gtm$geneSetName == unique(df.gtm$geneSetName)[7],]$tissue
  df.ts$status <- df.gtm[df.gtm$geneSetName == unique(df.gtm$geneSetName)[7],]$status
  df.ts$cell <- df.gtm[df.gtm$geneSetName==unique(df.gtm$geneSetName)[7],]$cell
  lp <- lapply(gs.names[grep('all',gs.names)][-4], function(y) {
    df.ts.int <- df.ts[order(df.ts[,y]),]
    df.ts.int$id <- as.numeric(factor(as.character(df.ts.int$cell), levels = (as.character(df.ts.int$cell))))
    #z <- setNames(c(0.55,0.05,0), gs.names[grep('all',gs.names)])
    z <- setNames(c(0,0,0), gs.names[grep('all',gs.names)][-4])
    ggplot(df.ts.int, aes_string(x='id', y=y)) + 
    geom_line() + geom_hline(yintercept = z[names(z) == y]) + theme_custom + xlab('Cell') + ylab(gsub('_',' ',y)) + 
    ggtitle(paste0('CD3+ lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) })
  df.ts$pos <- ifelse(df.ts$Tscm_vs_Tn_all > 0 & df.ts$Tscm_vs_Tcm_all > 0 & df.ts$Tscm_vs_Tem_all > 0, 'Yes', 'No')
  ts.int <- ggplot(df.ts %>% arrange(pos), aes_string(x='tSNE_1', y='tSNE_2', color='pos')) + 
    geom_point(shape = 16, size = ifelse(x=='blood',0.8,1), alpha = 0.6) + 
    scale_color_manual('Positive', values = c('grey80','royalblue3')) +
    scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~genotype, ncol = 4) + 
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) + facet_grid(status~tissue, scales='fixed', space='fixed') + 
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + theme(strip.text = element_text(size = 16, face = "plain", margin = margin(5,5,5,5))) + 
    ggtitle(paste0('CD3+ lymphocytes'), subtitle = 'Tscm vs Tn, Tcm and Tem') 
  return(list(ts = ts.int, lp = lp)) })

# #-----------------------------------------------------------------------------------------------------------------#
# # Plot results with statistical testing
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Calculate statistics for scores from each gene set
# ttest <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ttest <- as.matrix(pairwise.wilcox.test(df.gtm[df.gtm$geneSetName == y & df.gtm$donor == x,]$score, 
#                                             df.gtm[df.gtm$geneSetName == y & df.gtm$donor == x,]$tissue, p.adj = "bonferroni")$p.value)
#     ttest <- reshape2::melt(ttest)
#     #ttest[is.na(ttest)] <- 1
#     ttest <- ttest[!is.na(ttest$value),]
#     ttest$value <- as.numeric(formatC(ttest$value * length(unique(df.gtm$geneSetName)), format = 'g', digits = 2)) #bonferroni correction
#     ttest$value <- ifelse(ttest$value>1,1,ttest$value)
#     ttest$sig <- ifelse(ttest$value<0.05, '*', 'ns')
#     ttest$xs <- as.numeric(ttest$Var2)
#     ttest$xe <- as.numeric(ttest$Var1)
#     ttest <- ttest[order(ttest$xe, decreasing = T),]
#     ttest <- ttest[order(ttest$xs),]
#     ttest$yp <- 1:nrow(ttest)
#     return(ttest) }) })
# statistics <- list()
# statistics$t <- setNames(bind_rows(ttest$HSCT, .id='geneset')[,c(1:5)], nm=c('geneset','target','ref','adj.pval','sig'))
# #statistics$t <- statistics$t[order(statistics$t$sig),]
# statistics$t <- transform(statistics$t[order(statistics$t$adj.pval, decreasing = F),], row.names = 1:nrow(statistics$t))
# write.csv(statistics$t, file = paste0(dir.figs$t, '_statistics.csv'), row.names = T)
# 
# # Box plot for scores for subsets of both tissues
# #ifelse(x %in% levels(df.gtm$geneSetName)[c(1,4)],setNames(c('CD4','gdT'), c('CD4','gdT')),setNames(c('CD8','gdT'), c('CD8','gdT')))
# boxplotta <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='tissue', y='score', fill='tissue')) + 
#       geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(width = 0.2, height = 0, seed=1), show.legend = F) + 
#       geom_boxplot(outlier.colour = NA, color='grey40', lwd = 0.3, width=0.6, alpha = 0.6, show.legend = F) + 
#       geom_text(data = ttest[[x]][[y]], aes(y = 1.1 - 0.02*yp+(ya), x=(xs+xe-1)/2, 
#                                             label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 1, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttest[[x]][[y]], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = xs-0.1, xend = xe+1.1, 
#                        y = (1.03-0.05*yp+(ya)), yend = (1.03-0.05*yp+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm$score)<(-0.5),min(df.gtm$score)-0.01,-0.51),1.1+ya)) + 
#       #limits = c(ifelse(min(df.gtm[df.gtm$cluster==x,]$score)<0,min(df.gtm[df.gtm$cluster==x,]$score)-0.01,-0.01),1.01+ya)) + 
#       #coord_cartesian(ylim = c(-1,1), clip = "off") +
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,3,0.2,1),'cm'))#, axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Violin plot for scores for subsets of both tissues
# violinplotta <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='tissue', y='score', fill = 'tissue')) + 
#       #geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(0.2, seed=1), show.legend = F) + 
#       geom_violin(scale = "width", color = "black", size = 0.2, show.legend = F, alpha = 0.5) + 
#       geom_boxplot(outlier.colour = NA, color='black', lwd = 0.3, width=0.15, alpha = 0.8, show.legend = F) + 
#       geom_text(data = ttest[[x]][[y]], aes(y = 1.1 - 0.02*yp+(ya), x=(xs+xe-1)/2, 
#                                             label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 1, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttest[[x]][[y]], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = xs-0.1, xend = xe+1.1, 
#                        y = (1.03-0.05*yp+(ya)), yend = (1.03-0.05*yp+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)<(-0.5),min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)-0.01,-0.51),1.1+ya)) + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,3,0.2,1),'cm'))#, axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Plot frequency of scores for subsets of both tissues
# histplotta <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     #subsets <- levels(droplevels(df.gtm$cluster))[grep(y, levels(droplevels(df.gtm$cluster)))]
#     maxha <- ceiling(max(sapply(levels(df.gtm[df.gtm$geneSetName==y,]$tissue), function(w) sapply(x, function(z)
#       max(density(df.gtm[df.gtm$geneSetName==y & df.gtm$donor==z & df.gtm$tissue==w ,'score'])$y) ))))
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='score', color='tissue', fill='tissue')) + 
#       #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_vline(xintercept = 0, linetype = 'dashed') + 
#       geom_text(data = ttest[[x]][[y]], aes(y = 0.97*maxha, x=0.75, 
#                                             label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #annotate('text', x = 0.85, y = 1, size = 4, color = 'black', 
#       #         label = paste0('p=',formatC(t.test(df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='blood', 'score'], 
#       #                                            df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='skin', 'score'])$p.value,format = 'g', digits = 2))) +
#       #geom_vline(aes(xintercept = median(value)), linetype = 'dashed') +
#       #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#       #scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
#       scale_x_continuous('Score', limits = c(-1,1), breaks = c(seq(-1,1,0.5))) + 
#       scale_y_continuous('Density', limits = c(0,maxha+0.05*maxha), breaks = seq(0,maxha,maxha), expand = c(0,0)) + 
#       scale_fill_manual('Tissue', values = c(NA,col.cl[3])) + #c(blues[c(1,2,3,4)],NA)) + 
#       scale_color_manual('Tissue', values = rep('black', length(col.cl))) + 
#       #ggtitle('HSCT patient', subtitle = paste0('Gene set expression\nin host blood CD4 T cells')) + theme_custom + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) #+ facet_wrap(~sample, ncol = 1, scales = 'fixed', strip.position = 'right')
#     #facet_grid(cluster~., scales = 'fixed')
#   }) })
# 
# # Plot scores as split violin plots
# GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
#                            draw_group = function(self, data, ..., draw_quantiles = NULL) {
#                              data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
#                              grp <- data[1, "group"]
#                              newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
#                              newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
#                              newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
#                              
#                              if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
#                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
#                                                                          1))
#                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#                                aesthetics$alpha <- rep(1, nrow(quantiles))
#                                both <- cbind(quantiles, aesthetics)
#                                quantile_grob <- GeomPath$draw_panel(both, ...)
#                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#                              }
#                              else {
#                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
#                              }
#                            })
# geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
#                               draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
#                               show.legend = NA, inherit.aes = TRUE) {
#   layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
#         position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
#         params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
# }
# vp.gs <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,],  aes(x = gsub('_',' ',y), y = score, fill = tissue)) + 
#       geom_split_violin(scale = "width", color = "black", size = 0.2, show.legend = T, alpha = 0.6) + 
#       geom_hline(yintercept = 0, linetype = 'dashed') + 
#       #geom_text(data = ttest2[ttest2$geneSetName == x,], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2,
#       geom_text(data = ttest[[x]][[y]], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2, 
#                                             label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 1, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #geom_segment(data = ttest2[[x]], colour = "black", show.legend = F, size=0.3, 
#       #             #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#       #             aes(x = xs-0.1, xend = xe+1.1, 
#       #                 y = (1.03-0.1*yp+(ya)), yend = (1.03-0.1*yp+(ya))), inherit.aes = F) + 
#       #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2, seed=1), color = "black", show.legend = F) + 
#       #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.05, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size = 3) + 
#       #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
#       scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1.01,1+ya)) + 
#       scale_fill_manual('Tissue', values = col.cl[c(1,3)]) + xlab('') + ylab('Score') + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0(gsub('_',' ',y))) + #facet_wrap(~variable, ncol = 4) + 
#       #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
#       theme_custom + theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) }) })#,axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) })
# 
# # Calculate difference between average scores by tissue
# ds <- list()
# ds$t<- Map(x=levels(factor(df.gtm$donor)), function(x) { 
#   lapply(setNames((unique(statistics$t[statistics$t$adj.pval<0.05,]$geneset)), nm=(unique(statistics$t[statistics$t$adj.pval<0.05,]$geneset))), function(y) {
#     plyr::ddply(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], c('tissue'), summarise, avg = mean(score)) })})
# ds$t <- bind_rows(Map(x=levels(factor(df.gtm$donor)), function(x) bind_rows(ds$t[[x]], .id='geneset')), .id = 'donor')
# #ds$t <- plyr::ddply(ds$t, c('donor','geneset','tissue'), summarise, dif = last-first)
# ds$t <- ds$t %>% group_by(geneset) %>% summarize(dif = dplyr::last(avg)-dplyr::first(avg))
# ds$t <- transform(statistics$t[match(ds$t$geneset,statistics$t$geneset),], dif = ds$t$dif)
# ds$t$geneset <- sub('KEGG_','',ds$t$geneset)
# ds$t$geneset <- gsub('_',' ',ds$t$geneset)
# ds$t$geneset <- factor(ds$t$geneset, levels = unique(ds$t[order(ds$t$dif, decreasing = T),]$geneset))
# ds$t$logp <- -log10(ds$t$adj.pval)
# ds$t$logp <- ifelse(ds$t$logp>100,100,ds$t$logp)
# 
# # Plot bar plot
# #col.hm <- colorRampPalette(c('white', 'firebrick3', 'black'))(50)
# col.hm <- colorRampPalette(c('white', 'firebrick1', 'firebrick4','black'))(50)
# barplott <- #ggplot(ds$t[ds$t$geneset %in% ds$t$geneset[1:20],], aes(x = geneset, y = dif, fill = logp)) +
#   ggplot(ds$t, aes(x = geneset, y = dif, fill = logp)) +
#   geom_bar(stat = "identity", width=0.6, position=position_dodge(width = 0.7), color = 'black', alpha = 1) + #, color = 'grey80', fill = 'grey80') +
#   #geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
#   #              position = position_dodge(width = 0.7), width = 0.2, size = 0.5, show.legend = F) +
#   #geom_text(aes(y = percent+5, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5) +
#   #geom_text(data = ttest.cs, aes(y = max + 3, x = subset, label=ifelse(p.adj < 0.05,'*','ns')), #position = position_dodge(width = 0.7),
#   #          nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7),
#   #scale_x_discrete("", limits = levels(freqc$subset)) +
#   scale_y_continuous('Score [Skin - Blood]', limits = c(-0.31,0.31), breaks = seq(-0.4,0.4,0.1), expand = c(0,0)) + scale_x_discrete('') +
#   #scale_fill_manual("Group", values = col.sev) + xlab('Subset') +
#   scale_fill_gradientn(expression('-'*log[10]*'(adj.pval)'), colours = col.hm, 
#                        oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) +
#   #scale_fill_gradient(expression('-'*log[10]*'(adj.pval)'), low = 'white', high = 'firebrick3', na.value = 'red', 
#   #                    oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) +
#   ggtitle(paste('Skin vs blood lymphocytes'), subtitle = 'KEGG pathways') + theme_custom + #coord_flip(clip = "off") +
#   #theme(axis.text.y = element_text(size = 11, color = "black", angle = 0, hjust = 1), aspect.ratio = NULL)
#   theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust=0.5), aspect.ratio = NULL)
# 
# #=================================================================================================================#
# # Comparing gene set scores between genotypes for each tissue
# #=================================================================================================================#
# #-----------------------------------------------------------------------------------------------------------------#
# # Calculate module scores for gene sets using Seurat
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Split data by tissue
# #tp <- list(bd = subset(bd, subset = cluster_split %in% levels(bd$cluster_split)[c(6)]), 
# #           sk = subset(bd, subset = cluster_split %in% levels(bd$cluster_split)[c(1:4)]) )
# tg <- list(bd = subset(tp$tp, subset = orig.ident %in% levels(factor(tp$tp$orig.ident))[c(1)]), 
#            sk = subset(tp$tp, subset = orig.ident %in% levels(factor(tp$tp$orig.ident))[c(2)]) )
# 
# #-----------------------------------------------------------------------------------------------------------------#
# # Plot results with statistical testing
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Calculate statistics for scores from each gene set
# ttestg <- Map(x=sort(unique(df.gtm$geneSetName)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   ttest <- lapply(setNames(levels(factor(df.gtm$tissue)), nm=levels(factor(df.gtm$tissue))), function(y) {
#     ttest <- as.matrix(pairwise.wilcox.test(df.gtm[df.gtm$geneSetName == x & df.gtm$tissue == y,]$score, 
#                                             df.gtm[df.gtm$geneSetName == x & df.gtm$tissue == y,]$genotype, p.adj = "bonferroni")$p.value)
#     ttest <- reshape2::melt(ttest) })
#   ttest <- bind_rows(ttest, .id = 'tissue')
#   #ttest[is.na(ttest)] <- 1
#   ttest <- ttest[!is.na(ttest$value),]
#   ttest$value <- as.numeric(formatC(ttest$value * length(unique(df.gtm$geneSetName)) * length(unique(df.gtm$tissue)), format = 'g', digits = 2)) #bonferroni correction
#   ttest$value <- ifelse(ttest$value>1,1,ttest$value)
#   ttest$sig <- ifelse(ttest$value<0.05, '*', 'ns')
#   ttest$xs <- as.numeric(factor(ttest$tissue))
#   ttest$xe <- as.numeric(factor(ttest$tissue))
#   ttest <- ttest[order(ttest$xe, decreasing = T),]
#   ttest <- ttest[order(ttest$xs),]
#   ttest$yp <- 1:nrow(ttest)
#   return(ttest) })
# #ttestg <- bind_rows(Map(x=names(ttestg), function(x) bind_rows(ttestg[[x]], .id = 'geneSetName')), .id = 'tissue')
# statistics$tg <- setNames(bind_rows(ttestg, .id='geneset')[,c(1:6)], nm=c('geneset','tissue','target','ref','adj.pval','sig'))
# statistics$tg <- statistics$tg[order(statistics$tg$adj.pval, decreasing = F),]
# statistics$tg <- transform(statistics$tg[order(statistics$tg$tissue, decreasing = F),], row.names = 1:nrow(statistics$tg))
# write.csv(statistics$tg, file = paste0(dir.figs$tg, '_statistics.csv'), row.names = T)
# 
# # Box plot for scores for subsets of both tissues
# #ifelse(x %in% levels(df.gtm$geneSetName)[c(1,4)],setNames(c('CD4','gdT'), c('CD4','gdT')),setNames(c('CD8','gdT'), c('CD8','gdT')))
# boxplottg <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='tissue', y='score', fill='genotype')) + 
#       geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.85, seed=1), show.legend = F) + 
#       geom_boxplot(outlier.colour = NA, color='grey40', position=position_dodge(width = 0.85), lwd = 0.3, width=0.75, alpha = 0.6, show.legend = T) + 
#       geom_text(data = ttestg[[y]], aes(y = 1.1 +(ya), x=tissue, #(xs+xe-1)/2, 
#                                         label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestg[[y]], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(tissue))-0.1, xend = as.numeric(factor(tissue))+0.1, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm$score)<(-0.5),min(df.gtm$score)-0.01,-0.51),1.1+ya)) + 
#       #limits = c(ifelse(min(df.gtm[df.gtm$cluster==x,]$score)<0,min(df.gtm[df.gtm$cluster==x,]$score)-0.01,-0.01),1.01+ya)) + 
#       #coord_cartesian(ylim = c(-1,1), clip = "off") +
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'))#, axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Violin plot for scores for subsets of both tissues
# violinplottg <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='tissue', y='score', fill = 'genotype')) + 
#       #geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(0.2, seed=1), show.legend = F) + 
#       geom_violin(scale = "width", color = "black", size = 0.2, position=position_dodge(width = 0.85), width = 0.75, show.legend = F, alpha = 0.5) + 
#       geom_boxplot(outlier.colour = NA, color='black', position=position_dodge(width = 0.85), lwd = 0.3, width = 0.2, alpha = 0.8, show.legend = T) + 
#       geom_text(data = ttestg[[y]], aes(y = 1.1 +(ya), x=tissue, #(xs+xe-1)/2, 
#                                         label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestg[[y]], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(tissue))-0.25, xend = as.numeric(factor(tissue))+0.25, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)<(-0.5),min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)-0.01,-0.51),1.1+ya)) + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'))#, axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Plot frequency of scores for subsets of both tissues
# histplottg <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     #subsets <- levels(droplevels(df.gtm$cluster))[grep(y, levels(droplevels(df.gtm$cluster)))]
#     maxha <- ceiling(max(sapply(levels(df.gtm[df.gtm$geneSetName==y,]$tissue), function(w) sapply(x, function(z) sapply(levels(df.gtm$genotype), function(t)
#       max(density(df.gtm[df.gtm$geneSetName==y & df.gtm$donor==z & df.gtm$tissue==w & df.gtm$genotype==t,'score'])$y) )))))
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,], aes_string(x='score', color='genotype', fill='genotype')) + 
#       #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_vline(xintercept = 0, linetype = 'dashed') + 
#       geom_text(data = ttestg[[y]], aes(y = 0.97*maxha, x=0.75, 
#                                         label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #annotate('text', x = 0.85, y = 1, size = 4, color = 'black', 
#       #         label = paste0('p=',formatC(t.test(df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='blood', 'score'], 
#       #                                            df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='skin', 'score'])$p.value,format = 'g', digits = 2))) +
#       #geom_vline(aes(xintercept = median(value)), linetype = 'dashed') +
#       #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#       #scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
#       scale_x_continuous('Score', limits = c(-1,1), breaks = c(seq(-1,1,0.5))) + 
#       scale_y_continuous('Density', limits = c(0,maxha+0.05*maxha), breaks = seq(0,maxha,maxha), expand = c(0,0)) + 
#       scale_fill_manual('Tissue', values = c(NA,col.cl[3])) + #c(blues[c(1,2,3,4)],NA)) + 
#       scale_color_manual('Tissue', values = rep('black', length(col.cl))) + 
#       #ggtitle('HSCT patient', subtitle = paste0('Gene set expression\nin host blood CD4 T cells')) + theme_custom + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) + facet_wrap(~tissue, ncol = 1, scales = 'fixed', strip.position = 'right')
#     #facet_grid(cluster~., scales = 'fixed')
#   }) })
# 
# # Plot scores as split violin plots
# vp.gg <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x,],  aes(x = tissue, y = score, fill = genotype)) + 
#       geom_split_violin(scale = "width", color = "black", size = 0.2, show.legend = T, alpha = 0.6) + 
#       geom_hline(yintercept = 0, linetype = 'dashed') + 
#       #geom_text(data = ttest2[ttest2$geneSetName == x,], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2,
#       geom_text(data = ttestg[[y]], aes(y = 1 +(ya), x=tissue, #(xs+xe-1)/2, 
#                                         label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #geom_segment(data = ttest2[[x]], colour = "black", show.legend = F, size=0.3, 
#       #             #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#       #             aes(x = xs-0.1, xend = xe+1.1, 
#       #                 y = (1.03-0.1*yp+(ya)), yend = (1.03-0.1*yp+(ya))), inherit.aes = F) + 
#       #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2, seed=1), color = "black", show.legend = F) + 
#       #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.05, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size = 3) + 
#       #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
#       scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1.01,1+ya)) + 
#       scale_fill_manual('Tissue', values = col.cl[c(1,3)]) + xlab('') + ylab('Score') + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0(gsub('_',' ',y))) + #facet_wrap(~variable, ncol = 4) + 
#       #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
#       theme_custom + theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) }) })#,axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) })
# 
# # Calculate difference between average scores by tissue
# ds$tg<- Map(x=levels(factor(df.gtm$tissue)), function(x) { 
#   lapply(setNames((unique(statistics$tg[statistics$tg$tissue == x & statistics$tg$adj.pval<0.05,]$geneset)), 
#                   nm=(unique(statistics$tg[statistics$tg$tissue == x & statistics$tg$adj.pval<0.05,]$geneset))), function(y) {
#                     plyr::ddply(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x,], c('genotype'), summarise, avg = mean(score)) })})
# ds$tg <- bind_rows(Map(x=levels(factor(df.gtm$tissue)), function(x) bind_rows(ds$tg[[x]], .id='geneset')), .id = 'tissue')
# #ds$tg <- plyr::ddply(ds$tg, c('donor','geneset','tissue'), summarise, dif = last-first)
# ds$tg <- ds$tg %>% group_by(paste0(tissue,'_',geneset)) %>% summarize(tissue = unique(tissue), geneset = unique(geneset), dif = dplyr::last(avg)-dplyr::first(avg))
# ds$tg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { 
#   int <- statistics$tg[statistics$tg$tissue == x,]
#   transform(int[match(ds$tg[ds$tg$tissue == x,]$geneset, int$geneset),], dif = ds$tg[ds$tg$tissue == x,]$dif) })
# ds$tg <- bind_rows(ds$tg)#, .id = 'tissue')
# ds$tg$geneset <- sub('KEGG_','',ds$tg$geneset)
# ds$tg$geneset <- gsub('_',' ',ds$tg$geneset)
# #ds$tg$geneset <- factor(ds$tg$geneset, levels = unique(ds$tg[order(ds$tg$dif, decreasing = T),]$geneset))
# ds$tg$logp <- -log10(ds$tg$adj.pval)
# ds$tg$logp <- ifelse(ds$tg$logp>100,100,ds$tg$logp)
# 
# # Plot bar plot
# barplottg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { 
#   ggplot(ds$tg[ds$tg$tissue == x, ], aes(x = factor(geneset, levels = unique(ds$tg[ds$tg$tissue == x, ][order(ds$tg[ds$tg$tissue == x, ]$dif, decreasing = T),]$geneset)), y = dif, fill = logp)) + 
#     geom_bar(stat = "identity", width=0.6, position=position_dodge(width = 0.7), color = 'black', alpha = 1) + 
#     scale_y_continuous('Score [Host - Donor]', limits = c(-0.21,0.21), breaks = seq(-0.4,0.4,0.1), expand = c(0,0)) + scale_x_discrete('') + 
#     scale_fill_gradientn(expression('-'*log[10]*'(adj.pval)'), colours = col.hm, 
#                          #oob = scales::squish, limits = c(0,20), breaks = seq(0,20, length.out = 5)) + 
#                          oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) + 
#     ggtitle(paste0('Host vs donor lymphocytes (', x, ')'), subtitle = 'KEGG pathways') + theme_custom + #coord_flip(clip = "off") + 
#     #theme(axis.text.y = element_text(size = 11, color = "black", angle = 0, hjust = 1), aspect.ratio = NULL) })
#     theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust=0.5), aspect.ratio = NULL) })
# 
# 
# 
# 
# 
# #=================================================================================================================#
# # Comparing gene set scores between subsets for each tissue
# #=================================================================================================================#
# #-----------------------------------------------------------------------------------------------------------------#
# # Plot results with statistical testing
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Calculate statistics for scores from each gene set
# df.gtm$celltype <- factor(df.gtm$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,6,5,1)])
# ttests <- Map(x=sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   ttest <- lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ttest <- as.matrix(pairwise.wilcox.test(df.gtm[df.gtm$geneSetName == y & df.gtm$celltype == x,]$score, 
#                                             df.gtm[df.gtm$geneSetName == y & df.gtm$celltype == x,]$tissue, p.adj = "bonferroni")$p.value) })
#   transform(reshape2::melt(ttest), geneSetName = L1, L1 = NULL) })
# ttests <- bind_rows(ttests, .id = 'celltype')
# #ttest[is.na(ttest)] <- 1
# ttests <- ttests[!is.na(ttests$value),]
# ttests$value <- as.numeric(formatC(ttests$value * length(unique(df.gtm$geneSetName)) * length(unique(ttests$celltype)), format = 'g', digits = 2)) #bonferroni correction
# ttests$value <- ifelse(ttests$value>1,1,ttests$value)
# ttests$sig <- ifelse(ttests$value<0.05, '*', 'ns')
# ttests$xs <- as.numeric(factor(ttests$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)]))
# ttests$xe <- as.numeric(factor(ttests$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)]))
# ttests <- ttests[order(ttests$xe, decreasing = T),]
# ttests <- ttests[order(ttests$xs),]
# ttests$yp <- 1:nrow(ttests)
# ttests$celltype <- factor(ttests$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,6,5,1)])
# #statistics$ts <- setNames(bind_rows(ttests, .id='geneset')[,c(1:6)], nm=c('geneset','tissue','target','ref','adj.pval','sig'))
# statistics$ts <- setNames(ttests[,c(1,5,2,3,4,6)], nm=c('celltype','geneset','target','ref','adj.pval','sig'))
# statistics$ts <- statistics$ts[order(statistics$ts$adj.pval, decreasing = F),]
# statistics$ts <- transform(statistics$ts[order(statistics$ts$celltype, decreasing = F),], row.names = 1:nrow(statistics$ts))
# write.csv(statistics$ts, file = paste0(dir.figs$ts, '_statistics.csv'), row.names = T)
# 
# # Box plot for scores for subsets of both tissues
# #ifelse(x %in% levels(df.gtm$geneSetName)[c(1,4)],setNames(c('CD4','gdT'), c('CD4','gdT')),setNames(c('CD8','gdT'), c('CD8','gdT')))
# boxplotts <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x & df.gtm$celltype %in% unique(ttests$celltype),], aes_string(x='celltype', y='score', fill='tissue')) + 
#       geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.85, seed=1), show.legend = F) + 
#       geom_boxplot(outlier.colour = NA, color='grey40', position=position_dodge(width = 0.85), lwd = 0.3, width=0.75, alpha = 0.6, show.legend = T) +
#       geom_text(data = ttests[ttests$geneSetName == y,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                              label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttests[ttests$geneSetName == y,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(celltype))-0.1, xend = as.numeric(factor(celltype))+0.1, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm$score)<(-0.5),min(df.gtm$score)-0.01,-0.51),1.1+ya)) + 
#       #limits = c(ifelse(min(df.gtm[df.gtm$cluster==x,]$score)<0,min(df.gtm[df.gtm$cluster==x,]$score)-0.01,-0.01),1.01+ya)) + 
#       #coord_cartesian(ylim = c(-1,1), clip = "off") +
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Violin plot for scores for subsets of both tissues
# violinplotts <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x & df.gtm$celltype %in% unique(ttests$celltype),], aes_string(x='celltype', y='score', fill = 'tissue')) + 
#       #geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(0.2, seed=1), show.legend = F) + 
#       geom_violin(scale = "width", color = "black", size = 0.2, position=position_dodge(width = 0.85), width = 0.75, show.legend = F, alpha = 0.5) + 
#       geom_boxplot(outlier.colour = NA, color='black', position=position_dodge(width = 0.85), lwd = 0.3, width = 0.2, alpha = 0.8, show.legend = T) + 
#       geom_text(data = ttests[ttests$geneSetName == y,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                              label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttests[ttests$geneSetName == y,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(celltype))-0.25, xend = as.numeric(factor(celltype))+0.25, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)<(-0.5),min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)-0.01,-0.51),1.1+ya)) + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Plot frequency of scores for subsets of both tissues
# histplotts <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     #subsets <- levels(droplevels(df.gtm$cluster))[grep(y, levels(droplevels(df.gtm$cluster)))]
#     maxha <- ceiling(max(sapply(levels(df.gtm$celltype), function(w) sapply(x, function(z) 
#       max(density(df.gtm[df.gtm$geneSetName==y & df.gtm$donor==z & df.gtm$celltype==w,'score'])$y) ))))
#     #df.gtm$celltype <- factor(df.gtm$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)])
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x & df.gtm$celltype %in% unique(ttests$celltype),], aes_string(x='score', color='tissue', fill='tissue')) + 
#       #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_vline(xintercept = 0, linetype = 'dashed') + 
#       geom_text(data = ttests[ttests$geneSetName == y,], aes(y = 0.97*maxha, x=0.75, 
#                                                              label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #annotate('text', x = 0.85, y = 1, size = 4, color = 'black', 
#       #         label = paste0('p=',formatC(t.test(df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='blood', 'score'], 
#       #                                            df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='skin', 'score'])$p.value,format = 'g', digits = 2))) +
#       #geom_vline(aes(xintercept = median(value)), linetype = 'dashed') +
#       #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#       #scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
#       scale_x_continuous('Score', limits = c(-1,1), breaks = c(seq(-1,1,0.5))) + 
#       scale_y_continuous('Density', limits = c(0,maxha+0.05*maxha), breaks = seq(0,maxha,maxha), expand = c(0,0)) + 
#       scale_fill_manual('Tissue', values = c(NA,col.cl[3])) + #c(blues[c(1,2,3,4)],NA)) + 
#       scale_color_manual('Tissue', values = rep('black', length(col.cl))) + 
#       #ggtitle('HSCT patient', subtitle = paste0('Gene set expression\nin host blood CD4 T cells')) + theme_custom + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) + facet_wrap(~celltype, ncol = 1, scales = 'fixed', strip.position = 'right')
#     #facet_grid(cluster~., scales = 'fixed')
#   }) })
# 
# # Plot scores as split violin plots
# vp.gss <- Map(x=levels(factor(df.gtm$donor)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$donor %in% x & df.gtm$celltype %in% unique(ttests$celltype),],  aes(x = celltype, y = score, fill = tissue)) + 
#       geom_split_violin(scale = "width", color = "black", size = 0.2, show.legend = T, alpha = 0.6) + 
#       geom_hline(yintercept = 0, linetype = 'dashed') + 
#       #geom_text(data = ttest2[ttest2$geneSetName == x,], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2,
#       geom_text(data = ttests[ttests$geneSetName == y,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                              label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #geom_segment(data = ttest2[[x]], colour = "black", show.legend = F, size=0.3, 
#       #             #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#       #             aes(x = xs-0.1, xend = xe+1.1, 
#       #                 y = (1.03-0.1*yp+(ya)), yend = (1.03-0.1*yp+(ya))), inherit.aes = F) + 
#       #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2, seed=1), color = "black", show.legend = F) + 
#       #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.05, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size = 3) + 
#       #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
#       scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1.01,1.1+ya)) + 
#       scale_fill_manual('Tissue', values = col.cl[c(1,3)]) + xlab('') + ylab('Score') + 
#       ggtitle('HSCT patient lymphocytes', subtitle = paste0(gsub('_',' ',y))) + #facet_wrap(~variable, ncol = 4) + 
#       #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
#       theme_custom + theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) }) })
# 
# # Calculate difference between average scores by tissue
# ds$ts <- lapply(setNames(as.character(unique(ttests$celltype)), nm = as.character(unique(ttests$celltype))), function(y) { 
#   Map(x=setNames((unique(statistics$ts[statistics$ts$celltype == y & statistics$ts$adj.pval<0.05,]$geneset)), 
#                  nm=(unique(statistics$ts[statistics$ts$celltype == y & statistics$ts$adj.pval<0.05,]$geneset))), function(x) {
#                    plyr::ddply(df.gtm[df.gtm$geneSetName==x & df.gtm$celltype %in% y,], c('tissue'), summarise, avg = mean(score)) })})
# ds$ts <- bind_rows(Map(x=levels(factor(df.gtm$celltype)), function(x) bind_rows(ds$ts[[x]], .id='geneset')), .id = 'celltype')
# #ds$ts <- plyr::ddply(ds$ts, c('donor','geneset','tissue'), summarise, dif = last-first)
# ds$ts <- ds$ts %>% group_by(paste0(celltype,'_',geneset)) %>% summarize(celltype = unique(celltype), geneset = unique(geneset), dif = dplyr::last(avg)-dplyr::first(avg))
# ds$ts <- Map(x=as.character(unique(ttests$celltype)), function(x) { 
#   int <- statistics$ts[statistics$ts$celltype == x,]
#   transform(int[match(ds$ts[ds$ts$celltype == x,]$geneset, int$geneset),], dif = ds$ts[ds$ts$celltype == x,]$dif) })
# ds$ts <- bind_rows(ds$ts)#, .id = 'tissue')
# ds$ts$geneset <- sub('KEGG_','',ds$ts$geneset)
# ds$ts$geneset <- gsub('_',' ',ds$ts$geneset)
# #ds$ts$geneset <- factor(ds$ts$geneset, levels = unique(ds$ts[order(ds$ts$dif, decreasing = T),]$geneset))
# ds$ts$logp <- -log10(ds$ts$adj.pval)
# ds$ts$logp <- ifelse(ds$ts$logp>100,100,ds$ts$logp)
# 
# # Plot bar plot
# barplotts <- Map(x=as.character(unique(ds$ts$celltype)), function(x) { 
#   ggplot(ds$ts[ds$ts$celltype == x, ], aes(x = factor(geneset, levels = unique(ds$ts[ds$ts$celltype == x, ][order(ds$ts[ds$ts$celltype == x, ]$dif, decreasing = T),]$geneset)), y = dif, fill = logp)) + 
#     geom_bar(stat = "identity", width=0.6, position=position_dodge(width = 0.7), color = 'black', alpha = 1) + 
#     scale_y_continuous('Score [Skin - Blood]', breaks = seq(-0.4,0.4,0.1), expand = c(0,0.01)) + scale_x_discrete('') + 
#     scale_fill_gradientn(expression('-'*log[10]*'(adj.pval)'), colours = col.hm, 
#                          oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) + 
#     ggtitle(paste0('Skin vs blood ', x, ''), subtitle = 'KEGG pathways') + theme_custom + #coord_flip(clip = "off") + 
#     #theme(axis.text.y = element_text(size = 11, color = "black", angle = 0, hjust = 1), aspect.ratio = NULL) })
#     theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust=0.5), aspect.ratio = NULL) })
# 
# 
# 
# 
# #=================================================================================================================#
# # Comparing gene set scores between genotypes for each subsets of each tissue
# #=================================================================================================================#
# #-----------------------------------------------------------------------------------------------------------------#
# # Plot results with statistical testing
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Calculate statistics for scores from each gene set
# table(bd$celltype, paste(bd$orig.ident, bd$genotype))
# ttestsg <- Map(x=sort(as.character(unique(bd$celltype)))[c(2,3,4,7)], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   ttest <- lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ttest <- lapply(setNames(sort(levels(df.gtm$tissue)), nm=sort(levels(df.gtm$tissue))), function(z) { 
#       ttest <- as.matrix(pairwise.wilcox.test(df.gtm[df.gtm$geneSetName == y & df.gtm$celltype == x & df.gtm$tissue == z,]$score, 
#                                               df.gtm[df.gtm$geneSetName == y & df.gtm$celltype == x & df.gtm$tissue == z,]$genotype, p.adj = "bonferroni")$p.value) }) })
#   transform(reshape2::melt(ttest), geneSetName = L1, tissue = L2, L1 = NULL, L2 = NULL) })
# ttestsg <- bind_rows(ttestsg, .id = 'celltype')
# #ttest[is.na(ttest)] <- 1
# ttestsg <- ttestsg[!is.na(ttestsg$value),]
# ttestsg$value <- as.numeric(formatC(ttestsg$value * length(unique(df.gtm$geneSetName)) * length(unique(ttestsg$celltype)) * length(unique(ttestsg$tissue)), format = 'g', digits = 2)) #bonferroni correction
# ttestsg$value <- ifelse(ttestsg$value>1,1,ttestsg$value)
# ttestsg$sig <- ifelse(ttestsg$value<0.05, '*', 'ns')
# ttestsg$xs <- as.numeric(factor(ttestsg$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)]))
# ttestsg$xe <- as.numeric(factor(ttestsg$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)]))
# ttestsg <- ttestsg[order(ttestsg$xe, decreasing = T),]
# ttestsg <- ttestsg[order(ttestsg$xs),]
# ttestsg$yp <- 1:nrow(ttestsg)
# ttestsg$celltype <- factor(ttestsg$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,6,5,1)])
# #statistics$tsg <- setNames(bind_rows(ttestsg, .id='geneset')[,c(1:6)], nm=c('geneset','tissue','target','ref','adj.pval','sig'))
# statistics$tsg <- setNames(ttestsg[,c(6,1,5,2,3,4,7)], nm=c('tissue','celltype','geneset','target','ref','adj.pval','sig'))
# statistics$tsg <- statistics$tsg[order(statistics$tsg$adj.pval, decreasing = F),]
# statistics$tsg <- transform(statistics$tsg[order(statistics$tsg$celltype, decreasing = F),], row.names = 1:nrow(statistics$tsg))
# statistics$tsg <- transform(statistics$tsg[order(statistics$tsg$tissue, decreasing = F),], row.names = 1:nrow(statistics$tsg))
# write.csv(statistics$tsg, file = paste0(dir.figs$tsg, '_statistics.csv'), row.names = T)
# 
# # Box plot for scores for subsets of both tissues
# #ifelse(x %in% levels(df.gtm$geneSetName)[c(1,4)],setNames(c('CD4','gdT'), c('CD4','gdT')),setNames(c('CD8','gdT'), c('CD8','gdT')))
# boxplottsg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$celltype %in% unique(ttestsg$celltype),], aes_string(x='celltype', y='score', fill='genotype')) + 
#       geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.85, seed=1), show.legend = F) + 
#       geom_boxplot(outlier.colour = NA, color='grey40', position=position_dodge(width = 0.85), lwd = 0.3, width=0.75, alpha = 0.6, show.legend = T) +
#       geom_text(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(celltype))-0.1, xend = as.numeric(factor(celltype))+0.1, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm$score)<(-0.5),min(df.gtm$score)-0.01,-0.51),1.1+ya)) + 
#       #limits = c(ifelse(min(df.gtm[df.gtm$cluster==x,]$score)<0,min(df.gtm[df.gtm$cluster==x,]$score)-0.01,-0.01),1.01+ya)) + 
#       #coord_cartesian(ylim = c(-1,1), clip = "off") +
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Violin plot for scores for subsets of both tissues
# violinplottsg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$celltype %in% unique(ttestsg$celltype),], aes_string(x='celltype', y='score', fill = 'genotype')) + 
#       #geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(0.2, seed=1), show.legend = F) + 
#       geom_violin(scale = "width", color = "black", size = 0.2, position=position_dodge(width = 0.85), width = 0.75, show.legend = F, alpha = 0.5) + 
#       geom_boxplot(outlier.colour = NA, color='black', position=position_dodge(width = 0.85), lwd = 0.3, width = 0.2, alpha = 0.8, show.legend = T) + 
#       geom_text(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(celltype))-0.25, xend = as.numeric(factor(celltype))+0.25, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)<(-0.5),min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)-0.01,-0.51),1.1+ya)) + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       #facet_wrap(~tissue, ncol = 1, scales = 'free', strip.position = 'right') + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Plot frequency of scores for subsets of both tissues
# histplottsg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     #subsets <- levels(droplevels(df.gtm$cluster))[grep(y, levels(droplevels(df.gtm$cluster)))]
#     maxha <- ceiling(max(sapply(unique(ttestsg$celltype), function(w) sapply(x, function(z) sapply(x, function(z) sapply(levels(df.gtm$genotype), function(t)
#       max(density(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue==z & df.gtm$celltype==w & df.gtm$genotype==t,'score'])$y) ))))))
#     #df.gtm$celltype <- factor(df.gtm$celltype, levels = sort(as.character(unique(bd$celltype)))[c(2,3,4,7,5,1)])
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$t %in% x & df.gtm$celltype %in% unique(ttestsg$celltype),], aes_string(x='score', color='genotype', fill='genotype')) + 
#       #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_vline(xintercept = 0, linetype = 'dashed') + 
#       geom_text(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], aes(y = 0.97*maxha, x=0.75, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #annotate('text', x = 0.85, y = 1, size = 4, color = 'black', 
#       #         label = paste0('p=',formatC(t.test(df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='blood', 'score'], 
#       #                                            df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='skin', 'score'])$p.value,format = 'g', digits = 2))) +
#       #geom_vline(aes(xintercept = median(value)), linetype = 'dashed') +
#       #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#       #scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
#       scale_x_continuous('Score', limits = c(-1,1), breaks = c(seq(-1,1,0.5))) + 
#       scale_y_continuous('Density', limits = c(0,maxha+0.05*maxha), breaks = seq(0,maxha,maxha), expand = c(0,0)) + 
#       scale_fill_manual('Tissue', values = c(NA,col.cl[3])) + #c(blues[c(1,2,3,4)],NA)) + 
#       scale_color_manual('Tissue', values = rep('black', length(col.cl))) + 
#       #ggtitle('HSCT patient', subtitle = paste0('Gene set expression\nin host blood CD4 T cells')) + theme_custom + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) + facet_wrap(~celltype, ncol = 1, scales = 'fixed', strip.position = 'right')
#     #facet_grid(cluster~., scales = 'fixed')
#   }) })
# 
# # Plot scores as split violin plots
# vp.gssg <- Map(x=levels(factor(df.gtm$tissue)), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$celltype %in% unique(ttestsg$celltype),],  aes(x = celltype, y = score, fill = genotype)) + 
#       geom_split_violin(scale = "width", color = "black", size = 0.2, show.legend = T, alpha = 0.6) + 
#       geom_hline(yintercept = 0, linetype = 'dashed') + 
#       #geom_text(data = ttest2[ttest2$geneSetName == x,], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2,
#       geom_text(data = ttestsg[ttestsg$geneSetName == y & ttestsg$tissue %in% x,], aes(y = 1.1 +(ya), x=celltype, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #geom_segment(data = ttest2[[x]], colour = "black", show.legend = F, size=0.3, 
#       #             #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#       #             aes(x = xs-0.1, xend = xe+1.1, 
#       #                 y = (1.03-0.1*yp+(ya)), yend = (1.03-0.1*yp+(ya))), inherit.aes = F) + 
#       #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2, seed=1), color = "black", show.legend = F) + 
#       #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.05, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size = 3) + 
#       #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
#       scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1.01,1.1+ya)) + 
#       scale_fill_manual('Tissue', values = col.cl[c(1,3)]) + xlab('') + ylab('Score') + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
#       theme_custom + theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) }) })
# 
# # Calculate difference between average scores by tissue
# ds$tsg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue))), function(z) { 
#   lapply(setNames(as.character(unique(ttestsg$celltype)), nm = as.character(unique(ttestsg$celltype))), function(y) { 
#     Map(x=setNames((unique(statistics$tsg[statistics$tsg$celltype == y & statistics$tsg$tissue == z & statistics$tsg$adj.pval<0.05,]$geneset)), 
#                    nm=(unique(statistics$tsg[statistics$tsg$celltype == y & statistics$tsg$tissue == z & statistics$tsg$adj.pval<0.05,]$geneset))), function(x) {
#                      plyr::ddply(df.gtm[df.gtm$geneSetName==x & df.gtm$tissue==z & df.gtm$celltype %in% y,], c('genotype'), summarise, avg = mean(score)) })})})
# ds$tsg <- bind_rows(Map(x=levels(factor(df.gtm$tissue)), function(x) {
#   bind_rows(lapply(setNames(as.character(unique(ttestsg$celltype)), nm = as.character(unique(ttestsg$celltype))), function(y) { 
#     bind_rows(ds$tsg[[x]][[y]], .id='geneset')}), .id = 'celltype') }), .id = 'tissue')
# #ds$tsg <- plyr::ddply(ds$tsg, c('donor','geneset','tissue'), summarise, dif = last-first)
# ds$tsg <- ds$tsg %>% group_by(paste0(tissue,'_',celltype,'_',geneset)) %>% summarize(tissue = unique(tissue), celltype = unique(celltype), geneset = unique(geneset), dif = dplyr::last(avg)-dplyr::first(avg))
# ds$tsg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue))), function(z) { 
#   Map(x=as.character(unique(ttestsg$celltype)), function(x) { 
#     int <- statistics$tsg[statistics$tsg$celltype == x & statistics$tsg$tissue == z,]
#     transform(int[match(ds$tsg[ds$tsg$celltype == x & ds$tsg$tissue == z,]$geneset, int$geneset),], dif = ds$tsg[ds$tsg$celltype == x & ds$tsg$tissue == z,]$dif) }) })
# ds$tsg <- bind_rows(Map(x=levels(factor(df.gtm$tissue)), function(x) bind_rows(ds$tsg[[x]]) ))
# ds$tsg$geneset <- sub('KEGG_','',ds$tsg$geneset)
# ds$tsg$geneset <- gsub('_',' ',ds$tsg$geneset)
# #ds$tsg$geneset <- factor(ds$tsg$geneset, levels = unique(ds$tsg[order(ds$tsg$dif, decreasing = T),]$geneset))
# ds$tsg$logp <- -log10(ds$tsg$adj.pval)
# ds$tsg$logp <- ifelse(ds$tsg$logp>100,100,ds$tsg$logp)
# 
# # Plot bar plot
# barplottsg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue))), function(z) { 
#   Map(x=as.character(unique(ds$tsg[ds$tsg$tissue == z,]$celltype)), function(x) { 
#     ggplot(ds$tsg[ds$tsg$celltype == x & ds$tsg$tissue == z, ], aes(x = factor(geneset, levels = unique(ds$tsg[ds$tsg$celltype == x, ][order(ds$tsg[ds$tsg$celltype == x, ]$dif, decreasing = T),]$geneset)), y = dif, fill = logp)) + 
#       geom_bar(stat = "identity", width=0.6, position=position_dodge(width = 0.7), color = 'black', alpha = 1) + 
#       scale_y_continuous('Score [Host - Donor]', expand = c(0,0.01)) + scale_x_discrete('') + #, breaks = seq(-0.4,0.4,0.1)
#       scale_fill_gradientn(expression('-'*log[10]*'(adj.pval)'), colours = col.hm, 
#                            oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) + 
#       ggtitle(paste0('Host vs donor ', x, ''), subtitle = 'KEGG pathways') + theme_custom + #coord_flip(clip = "off") + 
#       #theme(axis.text.y = element_text(size = 11, color = "black", angle = 0, hjust = 1), aspect.ratio = NULL) })})
#       theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust=0.5), aspect.ratio = NULL) })})
# 
# 
# 
# 
# #=================================================================================================================#
# # Comparing gene set scores between genotypes for each subsets of each tissue
# #=================================================================================================================#
# #-----------------------------------------------------------------------------------------------------------------#
# # Plot results with statistical testing
# #-----------------------------------------------------------------------------------------------------------------#
# 
# # Define cells to include as control for comparison (only memory cells for blood, where host cells fall within the 
# # donor cells, and all cells for skin)
# bd[['include_bs']] <- bd$subset %in% levels(factor(bd$subset))[c(3,6)]#[c(2,3,5,6,8,11)]
# table(bd$include_bs, bd$subset)
# table(bd$genotype[bd$include_bs], bd$subset[bd$include_bs])
# tg$bd[['include_bs']] <- tg$bd$subset %in% levels(factor(tg$bd$subset))[c(2,4)]
# table(tg$bd$include_bs, tg$bd$subset)
# table(tg$bd$genotype[tg$bd$include_bs], tg$bd$subset[tg$bd$include_bs])
# #tb <- subset(tg$bd, subset = include_bs == T)
# 
# # Calculate statistics for scores from each gene set
# table(tg$bd$subset, paste(tg$bd$orig.ident, tg$bd$genotype))
# ttestsb <- Map(x=sort(as.character(unique(subset(tg$bd, subset = include_bs == T)$subset))), function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   ttest <- lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     lapply(setNames(sort(levels(df.gtm$tissue)), nm=sort(levels(df.gtm$tissue)))[1], function(z) { 
#       as.matrix(pairwise.wilcox.test(df.gtm[df.gtm$geneSetName == y & df.gtm$subset == x & df.gtm$tissue == z,]$score, 
#                                               df.gtm[df.gtm$geneSetName == y & df.gtm$subset == x & df.gtm$tissue == z,]$genotype, p.adj = "bonferroni")$p.value) }) })
#   transform(reshape2::melt(ttest), geneSetName = L1, tissue = L2, L1 = NULL, L2 = NULL) })
# ttestsb <- bind_rows(ttestsb, .id = 'subset')
# #ttest[is.na(ttest)] <- 1
# ttestsb <- ttestsb[!is.na(ttestsb$value),]
# ttestsb$value <- as.numeric(formatC(ttestsb$value * length(unique(ttestsb$geneSetName)) * length(unique(ttestsb$subset)) * length(unique(ttestsb$tissue)), format = 'g', digits = 2)) #bonferroni correction
# ttestsb$value <- ifelse(ttestsb$value>1,1,ttestsb$value)
# ttestsb$sig <- ifelse(ttestsb$value<0.05, '*', 'ns')
# ttestsb$xs <- as.numeric(factor(ttestsb$subset, levels = sort(as.character(unique(bd$subset)))[c(3,6)]))
# ttestsb$xe <- as.numeric(factor(ttestsb$subset, levels = sort(as.character(unique(bd$subset)))[c(3,6)]))
# ttestsb <- ttestsb[order(ttestsb$xe, decreasing = T),]
# ttestsb <- ttestsb[order(ttestsb$xs),]
# ttestsb$yp <- 1:nrow(ttestsb)
# ttestsb$subset <- factor(ttestsb$subset, levels = sort(as.character(unique(bd$subset)))[c(3,6)])
# #statistics$bmg <- setNames(bind_rows(ttestsb, .id='geneset')[,c(1:6)], nm=c('geneset','tissue','target','ref','adj.pval','sig'))
# statistics$bmg <- setNames(ttestsb[,c(6,1,5,2,3,4,7)], nm=c('tissue','subset','geneset','target','ref','adj.pval','sig'))
# statistics$bmg <- statistics$bmg[order(statistics$bmg$adj.pval, decreasing = F),]
# statistics$bmg <- transform(statistics$bmg[order(statistics$bmg$subset, decreasing = F),], row.names = 1:nrow(statistics$bmg))
# statistics$bmg <- transform(statistics$bmg[order(statistics$bmg$tissue, decreasing = F),], row.names = 1:nrow(statistics$bmg))
# write.csv(statistics$bmg, file = paste0(dir.figs$bmg, '_statistics.csv'), row.names = T)
# 
# # Box plot for scores for subsets of both tissues
# #ifelse(x %in% levels(df.gtm$geneSetName)[c(1,4)],setNames(c('CD4','gdT'), c('CD4','gdT')),setNames(c('CD8','gdT'), c('CD8','gdT')))
# boxplottb <- Map(x=levels(factor(df.gtm$tissue))[1], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$subset %in% unique(ttestsb$subset),], aes_string(x='subset', y='score', fill='genotype')) + 
#       geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.85, seed=1), show.legend = F) + 
#       geom_boxplot(outlier.colour = NA, color='grey40', position=position_dodge(width = 0.85), lwd = 0.3, width=0.75, alpha = 0.6, show.legend = T) +
#       geom_text(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], aes(y = 1.1 +(ya), x=subset, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(subset))-0.1, xend = as.numeric(factor(subset))+0.1, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm$score)<(-0.5),min(df.gtm$score)-0.01,-0.51),1.1+ya)) + 
#       #limits = c(ifelse(min(df.gtm[df.gtm$cluster==x,]$score)<0,min(df.gtm[df.gtm$cluster==x,]$score)-0.01,-0.01),1.01+ya)) + 
#       #coord_cartesian(ylim = c(-1,1), clip = "off") +
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Violin plot for scores for subsets of both tissues
# violinplottb <- Map(x=levels(factor(df.gtm$tissue))[1], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1#1-max(df.gtm[df.gtm$cluster==x,]$score)
#     #subsets <- levels(df.gtm$cluster)[grep(y, levels(df.gtm$cluster))]
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$subset %in% unique(ttestsb$subset),], aes_string(x='subset', y='score', fill = 'genotype')) + 
#       #geom_jitter(size = 1, shape = 21, alpha = 1, stroke = 0.2, position=position_jitter(0.2, seed=1), show.legend = F) + 
#       geom_violin(scale = "width", color = "black", size = 0.2, position=position_dodge(width = 0.85), width = 0.75, show.legend = F, alpha = 0.5) + 
#       geom_boxplot(outlier.colour = NA, color='black', position=position_dodge(width = 0.85), lwd = 0.3, width = 0.2, alpha = 0.8, show.legend = T) + 
#       geom_text(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], aes(y = 1.1 +(ya), x=subset, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       geom_segment(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], colour = "black", show.legend = F, size=0.3, 
#                    #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#                    aes(x = as.numeric(factor(subset))-0.25, xend = as.numeric(factor(subset))+0.25, #x = xs-0.1, xend = xe+1.1, 
#                        y = (0.98+(ya)), yend = (0.98+(ya))), inherit.aes = F) + 
#       scale_fill_manual('', values = col.cl[c(1,3)]) + #scale_shape_manual('', values = c(21,24)) + 
#       #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1,0.5), limits = c(-1.1,1.1)) + 
#       #scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.2), limits = c(-0.6,1.01)) +
#       scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(-1,1+ya,0.5), limits = c(-1.01,1.1+ya)) + #c(ifelse(min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)<(-0.5),min(df.gtm[df.gtm$geneSetName==x & df.gtm$cluster %in% subsets,]$score)-0.01,-0.51),1.1+ya)) + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       #facet_wrap(~tissue, ncol = 1, scales = 'free', strip.position = 'right') + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) #+ facet_wrap(~variable, ncol = 5, scales = 'free')
#   }) })
# 
# # Plot frequency of scores for subsets of both tissues
# histplottb <- Map(x=levels(factor(df.gtm$tissue))[1], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     #subsets <- levels(droplevels(df.gtm$cluster))[grep(y, levels(droplevels(df.gtm$cluster)))]
#     maxha <- ceiling(max(sapply(levels(ttestsb$subset), function(w) sapply(x, function(z) sapply(x, function(z) sapply(levels(df.gtm$genotype), function(t)
#       max(density(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue==z & df.gtm$subset==w & df.gtm$genotype==t,'score'])$y) ))))))
#     #df.gtm$subset <- factor(df.gtm$subset, levels = sort(as.character(unique(bd$subset)))[c(2,3,4,7,5,1)])
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$t %in% x & df.gtm$subset %in% unique(ttestsb$subset),], aes_string(x='score', color='genotype', fill='genotype')) + 
#       #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
#       geom_vline(xintercept = 0, linetype = 'dashed') + 
#       geom_text(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], aes(y = 0.97*maxha, x=0.75, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3.5, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #annotate('text', x = 0.85, y = 1, size = 4, color = 'black', 
#       #         label = paste0('p=',formatC(t.test(df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='blood', 'score'], 
#       #                                            df.gta[[geneSetName]][df.gta[[geneSetName]]$sample=='skin', 'score'])$p.value,format = 'g', digits = 2))) +
#       #geom_vline(aes(xintercept = median(value)), linetype = 'dashed') +
#       #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#       #scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
#       scale_x_continuous('Score', limits = c(-1,1), breaks = c(seq(-1,1,0.5))) + 
#       scale_y_continuous('Density', limits = c(0,maxha+0.05*maxha), breaks = seq(0,maxha,maxha), expand = c(0,0)) + 
#       scale_fill_manual('Tissue', values = c(NA,col.cl[3])) + #c(blues[c(1,2,3,4)],NA)) + 
#       scale_color_manual('Tissue', values = rep('black', length(col.cl))) + 
#       #ggtitle('HSCT patient', subtitle = paste0('Gene set expression\nin host blood CD4 T cells')) + theme_custom + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm')) + facet_wrap(~subset, ncol = 1, scales = 'fixed', strip.position = 'right')
#     #facet_grid(cluster~., scales = 'fixed')
#   }) })
# 
# # Plot scores as split violin plots
# vp.gsb <- Map(x=levels(factor(df.gtm$tissue))[1], function(x) { #lapply(setNames(c('CD4','CD8','gdT'), c('CD4','CD8','gdT')), function(y) {
#   lapply(setNames(sort(unique(df.gtm$geneSetName)), nm=sort(unique(df.gtm$geneSetName))), function(y) {
#     ya <- 0.1
#     ggplot(df.gtm[df.gtm$geneSetName==y & df.gtm$tissue %in% x & df.gtm$subset %in% unique(ttestsb$subset),],  aes(x = subset, y = score, fill = genotype)) + 
#       geom_split_violin(scale = "width", color = "black", size = 0.2, show.legend = T, alpha = 0.6) + 
#       geom_hline(yintercept = 0, linetype = 'dashed') + 
#       #geom_text(data = ttest2[ttest2$geneSetName == x,], aes(y = 1.1 - 0.1*yp+(ya), x=(xs+xe-2)/2,
#       geom_text(data = ttestsb[ttestsb$geneSetName == y & ttestsb$tissue %in% x,], aes(y = 1.1 +(ya), x=subset, #(xs+xe-1)/2, 
#                                                                                        label=paste0('p=',value)), #position = position_dodge(width = 0.7), 
#                 nudge_x = 0, size = 3, color = "black", inherit.aes = F) + #position = position_dodge(width = 0.7), 
#       #geom_segment(data = ttest2[[x]], colour = "black", show.legend = F, size=0.3, 
#       #             #aes(x = ifelse(value<0.05, xs-0.1, NA), xend = ifelse(value<0.05, xe+1.1, NA), 
#       #             aes(x = xs-0.1, xend = xe+1.1, 
#       #                 y = (1.03-0.1*yp+(ya)), yend = (1.03-0.1*yp+(ya))), inherit.aes = F) + 
#       #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2, seed=1), color = "black", show.legend = F) + 
#       #geom_text(data=de$tp.mg[de$tp.mg$gene %in% genelist,], aes(x=gene,y=5.05, label=paste0('p=',formatC(p_val_adj,format = 'g', digits = 2))), inherit.aes = F, size = 3) + 
#       #scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + 
#       scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1.01,1.1+ya)) + 
#       scale_fill_manual('Tissue', values = col.cl[c(1,3)]) + xlab('') + ylab('Score') + 
#       ggtitle(paste0('HSCT ',x,' lymphocytes'), subtitle = paste0('', gsub('_',' ',y))) + theme_custom + 
#       #ggtitle('HSCT skin', subtitle = paste0("Y-linked gene: ", x)) + #facet_wrap(~variable, ncol = 4) + 
#       theme_custom + theme(aspect.ratio = NULL, plot.margin = unit(c(0.2,1,0.2,2),'cm'), axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1)) }) })
# 
# # Calculate difference between average scores
# ds$bmg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue)))[1], function(z) { 
#   lapply(setNames(as.character(unique(ttestsb$subset)), nm = as.character(unique(ttestsb$subset))), function(y) { 
#     Map(x=setNames((unique(statistics$bmg[statistics$bmg$subset == y & statistics$bmg$tissue == z & statistics$bmg$adj.pval<0.05,]$geneset)), 
#                    nm=(unique(statistics$bmg[statistics$bmg$subset == y & statistics$bmg$tissue == z & statistics$bmg$adj.pval<0.05,]$geneset))), function(x) {
#                      plyr::ddply(df.gtm[df.gtm$geneSetName==x & df.gtm$tissue==z & df.gtm$subset %in% y,], c('genotype'), summarise, avg = mean(score)) })})})
# ds$bmg <- bind_rows(Map(x=levels(factor(df.gtm$tissue)), function(x) {
#   bind_rows(lapply(setNames(as.character(unique(ttestsb$subset)), nm = as.character(unique(ttestsb$subset))), function(y) { 
#     bind_rows(ds$bmg[[x]][[y]], .id='geneset')}), .id = 'subset') }), .id = 'tissue')
# #ds$bmg <- plyr::ddply(ds$bmg, c('donor','geneset','tissue'), summarise, dif = last-first)
# ds$bmg <- ds$bmg %>% group_by(paste0(tissue,'_',subset,'_',geneset)) %>% summarize(tissue = unique(tissue), subset = unique(subset), geneset = unique(geneset), dif = dplyr::last(avg)-dplyr::first(avg))
# ds$bmg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue))), function(z) { 
#   Map(x=as.character(unique(ttestsb$subset)), function(x) { 
#     int <- statistics$bmg[statistics$bmg$subset == x & statistics$bmg$tissue == z,]
#     transform(int[match(ds$bmg[ds$bmg$subset == x & ds$bmg$tissue == z,]$geneset, int$geneset),], dif = ds$bmg[ds$bmg$subset == x & ds$bmg$tissue == z,]$dif) }) })
# ds$bmg <- bind_rows(Map(x=levels(factor(df.gtm$tissue)), function(x) bind_rows(ds$bmg[[x]]) ))
# ds$bmg$geneset <- sub('KEGG_','',ds$bmg$geneset)
# ds$bmg$geneset <- gsub('_',' ',ds$bmg$geneset)
# #ds$bmg$geneset <- factor(ds$bmg$geneset, levels = unique(ds$bmg[order(ds$bmg$dif, decreasing = T),]$geneset))
# ds$bmg$logp <- -log10(ds$bmg$adj.pval)
# ds$bmg$logp <- ifelse(ds$bmg$logp>100,100,ds$bmg$logp)
# 
# # Plot bar plot
# barplotbmg <- lapply(setNames(levels(factor(df.gtm$tissue)), nm = levels(factor(df.gtm$tissue)))[1], function(z) { 
#   Map(x=as.character(unique(ds$bmg[ds$bmg$tissue == z,]$subset)), function(x) { 
#     ggplot(ds$bmg[ds$bmg$subset == x & ds$bmg$tissue == z, ], aes(x = factor(geneset, levels = unique(ds$bmg[ds$bmg$subset == x, ][order(ds$bmg[ds$bmg$subset == x, ]$dif, decreasing = T),]$geneset)), y = dif, fill = logp)) + 
#       geom_bar(stat = "identity", width=0.6, position=position_dodge(width = 0.7), color = 'black', alpha = 1) + 
#       scale_y_continuous('Score [Host - Donor]', expand = c(0,0.01)) + scale_x_discrete('') + #, breaks = seq(-0.4,0.4,0.1)
#       scale_fill_gradientn(expression('-'*log[10]*'(adj.pval)'), colours = col.hm, 
#                            oob = scales::squish, limits = c(0,100), breaks = seq(0,100, length.out = 5)) + 
#       ggtitle(paste0('Host vs donor ', x, ''), subtitle = 'KEGG pathways') + theme_custom + #coord_flip(clip = "off") + 
#       #theme(axis.text.y = element_text(size = 11, color = "black", angle = 0, hjust = 1), aspect.ratio = NULL) })})
#       theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust=0.5), aspect.ratio = NULL) })})





#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Save tSNE plots
lapply(names(ts), function(x) ggsave(paste0(dir.figs$tsne,'_',x,".pdf"), ts[[x]], height = 6, width = 6.5) )
ggsave(paste0(dir.figs$tsne,'_triple_positive.pdf'), tsp$tp$ts, height = 6, width = 6.5)
lapply(names(tsp$tp$lp), function(x) ggsave(paste0(dir.figs$tsne,'_distribution_',x,".pdf"), tsp$tp$lp[[x]], height = 4, width = 4) )

# # Save boxplots for comparison between tissues
# #lapply(names(boxplotm), function(x) { ggsave(paste0(dir.figs$all, '_boxplot_',gsub(' ','_',x),'.pdf'), boxplotm[[x]], height = 3.8, width = 5) })
# #lapply(names(histplotm), function(x) { ggsave(paste0(dir.figs$all, '_hist_',gsub(' ','_',x),'.pdf'), histplotm[[x]], height = 5, width = 4.5) })
# #lapply(names(vp.gs), function(x) { ggsave(paste0(dir.figs$all, '_violin_',gsub(' ','_',x),'.pdf'), vp.gs[[x]], height = 3.8, width = 4.5) })
# #lapply(names(boxplotta), function(x) { lapply(names(boxplotta[[x]]), function(y) { ggsave(paste0(dir.figs[['all']], '_boxplot_',gsub(' ','_',x),'','.pdf'), boxplotta[[x]][[y]], height = 3.5, width = 3.5) })})
# #lapply(names(histplotta), function(x) { lapply(names(histplotta[[x]]), function(y) { ggsave(paste0(dir.figs[['all']], '_hist_',gsub(' ','_',x),'','.pdf'), histplotta[[x]][[y]], height = 3.4, width = 4.5) })})
# #lapply(names(violinplotta), function(x) { lapply(names(violinplotta[[x]]), function(y) { ggsave(paste0(dir.figs[['all']], '_violin_',sub(' ','_',x),'','.pdf'), violinplotta[[x]][[y]], height = 3.5, width = 3.5) })})
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_boxplot_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = boxplotta[[x]], nrow = 3, ncol = 3), height = 11, width = 12) })
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_hist_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = histplotta[[x]], nrow = 3, ncol = 3), height = 11, width = 14) })
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_violin_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplotta[[x]], nrow = 3, ncol = 3), height = 11, width = 12) })
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_boxplot_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = boxplotta[[x]], nrow = 2, ncol = 2), height = 8, width = 8) })
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_hist_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = histplotta[[x]], nrow = 2, ncol = 2), height = 8, width = 10) })
# #lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs[['all']], '_violin_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplotta[[x]], nrow = 2, ncol = 2), height = 8, width = 8) })
# lapply(names(boxplotta), function(x) { ggsave(paste0(dir.figs$t, '_boxplot','','.pdf'), gridExtra::marrangeGrob(grobs = boxplotta[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 4.5) })
# lapply(names(histplotta), function(x) { ggsave(paste0(dir.figs$t, '_hist','','.pdf'), gridExtra::marrangeGrob(grobs = histplotta[[x]], nrow = 1, ncol = 1, top = NULL), height = 3, width = 5.5) })
# lapply(names(violinplotta), function(x) { ggsave(paste0(dir.figs$t, '_violin_split','','.pdf'), gridExtra::marrangeGrob(grobs = violinplotta[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 4.5) })
# lapply(names(vp.gs), function(x) { ggsave(paste0(dir.figs$t, '_violin_merged','','.pdf'), gridExtra::marrangeGrob(grobs = vp.gs[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 4.5) })
# ggsave(paste0(dir.figs$t, '_barplot','.pdf'), barplott, height = 8.5, width = 25)
# 
# # Save boxplots for comparison between genotype for each tissue
# lapply(names(boxplottg), function(x) { ggsave(paste0(dir.figs$tg, '_boxplot','','.pdf'), gridExtra::marrangeGrob(grobs = boxplottg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 5.5) })
# lapply(names(histplottg), function(x) { ggsave(paste0(dir.figs$tg, '_hist','','.pdf'), gridExtra::marrangeGrob(grobs = histplottg[[x]], nrow = 1, ncol = 1, top = NULL), height = 4.5, width = 5.8) })
# lapply(names(violinplottg), function(x) { ggsave(paste0(dir.figs$tg, '_violin_split','','.pdf'), gridExtra::marrangeGrob(grobs = violinplottg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 5.5) })
# lapply(names(vp.gg), function(x) { ggsave(paste0(dir.figs$tg, '_violin_merged','','.pdf'), gridExtra::marrangeGrob(grobs = vp.gg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 5.5) })
# lapply(names(barplottg)[1], function(x) { ggsave(paste0(dir.figs$tg, '_barplot_',x,'.pdf'), barplottg[[x]], height = 6.8, width = 4) })
# lapply(names(barplottg)[2], function(x) { ggsave(paste0(dir.figs$tg, '_barplot_',x,'.pdf'), barplottg[[x]], height = 8.3, width = 9.2) })
# 
# # Save boxplots for comparison between tissue for each subset
# lapply(names(boxplotts), function(x) { ggsave(paste0(dir.figs$ts, '_boxplot','','.pdf'), gridExtra::marrangeGrob(grobs = boxplotts[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.8, width = 8) })
# lapply(names(histplotts), function(x) { ggsave(paste0(dir.figs$ts, '_hist','','.pdf'), gridExtra::marrangeGrob(grobs = histplotts[[x]], nrow = 1, ncol = 1, top = NULL), height = 7, width = 5.8) })
# lapply(names(violinplotts), function(x) { ggsave(paste0(dir.figs$ts, '_violin_split','','.pdf'), gridExtra::marrangeGrob(grobs = violinplotts[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.8, width = 8) })
# lapply(names(vp.gss), function(x) { ggsave(paste0(dir.figs$ts, '_violin_merged','','.pdf'), gridExtra::marrangeGrob(grobs = vp.gss[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.8, width = 8) })
# #lapply(names(barplotts), function(x) { ggsave(paste0(dir.figs$ts, '_barplot_',x,'.pdf'), barplotts[[x]], height = 8.5, width = 25) })
# ggsave(paste0(dir.figs$ts, '_barplot','.pdf'), gridExtra::marrangeGrob(grobs = barplotts, nrow = 1, ncol = 1, top = NULL), height = 8.5, width = 25)
# 
# # Save boxplots for comparison between genotypes for each subset of each tissue
# lapply(names(boxplottsg), function(x) { ggsave(paste0(dir.figs$tsg, '_boxplot_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = boxplottsg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(histplottsg), function(x) { ggsave(paste0(dir.figs$tsg, '_hist_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = histplottsg[[x]], nrow = 1, ncol = 1, top = NULL), height = 5.2, width = 5.8) })
# lapply(names(violinplottsg), function(x) { ggsave(paste0(dir.figs$tsg, '_violin_split_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplottsg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(vp.gssg), function(x) { ggsave(paste0(dir.figs$tsg, '_violin_merged_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = vp.gssg[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(barplottsg), function(x) { ggsave(paste0(dir.figs$tsg, '_barplot_',x,'.pdf'), gridExtra::marrangeGrob(grobs = barplottsg[[x]], nrow = 1, ncol = 1, top = NULL), height = 7.5, width = 5) })
# 
# # Save boxplots for comparison between genotypes for each memory subset of blood
# lapply(names(boxplottb), function(x) { ggsave(paste0(dir.figs$bmg, '_boxplot_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = boxplottb[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(histplottb), function(x) { ggsave(paste0(dir.figs$bmg, '_hist_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = histplottb[[x]], nrow = 1, ncol = 1, top = NULL), height = 5.2, width = 5.8) })
# lapply(names(violinplottb), function(x) { ggsave(paste0(dir.figs$bmg, '_violin_split_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplottb[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(vp.gsb), function(x) { ggsave(paste0(dir.figs$bmg, '_violin_merged_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = vp.gsb[[x]], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(barplotbmg), function(x) { ggsave(paste0(dir.figs$bmg, '_barplot_',x,'.pdf'), gridExtra::marrangeGrob(grobs = barplotbmg[[x]], nrow = 1, ncol = 1, top = NULL), height = 4, width = 3.2) })
# 
# # Save violin plots only for significant gene sets
# lapply(names(violinplotta), function(x) { ggsave(paste0(dir.figs$t, '_violin_split_sig','','.pdf'), gridExtra::marrangeGrob(grobs = violinplotta[[x]][statistics$t[statistics$t$adj.pval<0.05,]$geneset], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 4.5) })
# lapply(names(violinplottg), function(x) { ggsave(paste0(dir.figs$tg, '_violin_split_sig','','.pdf'), gridExtra::marrangeGrob(grobs = violinplottg[[x]][statistics$tg[statistics$tg$adj.pval<0.05,]$geneset], nrow = 1, ncol = 1, top = NULL), height = 3.5, width = 5.5) })
# lapply(names(violinplotts), function(x) { ggsave(paste0(dir.figs$ts, '_violin_split_sig','','.pdf'), gridExtra::marrangeGrob(grobs = violinplotts[[x]][statistics$ts[statistics$ts$adj.pval<0.05,]$geneset], nrow = 1, ncol = 1, top = NULL), height = 3.8, width = 8) })
# lapply(names(violinplottsg), function(x) { ggsave(paste0(dir.figs$tsg, '_violin_split_sig_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplottsg[[x]][statistics$tsg[statistics$tsg$adj.pval<0.05 & statistics$tsg$tissue == x,]$geneset], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })
# lapply(names(violinplottb), function(x) { ggsave(paste0(dir.figs$bmg, '_violin_split_sig_',gsub(' ','_',x),'','.pdf'), gridExtra::marrangeGrob(grobs = violinplottb[[x]][statistics$bmg[statistics$bmg$adj.pval<0.05,]$geneset], nrow = 1, ncol = 1, top = NULL), height = 3.75, width = 6.5) })

#-----------------------------------------------------------------------------------------------------------------#
# End of script
#-----------------------------------------------------------------------------------------------------------------#

#=================================================================================================================#
# R session info
#=================================================================================================================#

#> devtools::session_info()
#─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#setting  value                       
#version  R version 3.5.1 (2018-07-02)
#os       macOS  10.14.6              
#system   x86_64, darwin15.6.0        
#ui       RStudio                     
#language (EN)                        
#collate  en_US.UTF-8                 
#ctype    en_US.UTF-8                 
#tz       Europe/Berlin               
#date     2020-06-24                  
#
#─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#package              * version    date       lib source                           
#abind                  1.4-5      2016-07-21 [1] CRAN (R 3.5.0)                   
#annotate             * 1.60.0     2018-10-30 [1] Bioconductor                     
#AnnotationDbi        * 1.44.0     2018-10-30 [1] Bioconductor                     
#assertthat             0.2.0      2017-04-11 [1] CRAN (R 3.5.0)                   
#backports              1.1.2      2017-12-13 [1] CRAN (R 3.5.0)                   
#bibtex                 0.4.2      2017-06-30 [1] CRAN (R 3.5.0)                   
#bindr                  0.1.1      2018-03-13 [1] CRAN (R 3.5.0)                   
#bindrcpp               0.2.2      2018-03-29 [1] CRAN (R 3.5.0)                   
#Biobase              * 2.42.0     2018-10-30 [1] Bioconductor                     
#BiocGenerics         * 0.28.0     2018-10-30 [1] Bioconductor                     
#BiocParallel           1.16.2     2018-11-28 [1] Bioconductor                     
#bit                    1.1-14     2018-05-29 [1] CRAN (R 3.5.0)                   
#bit64                  0.9-7      2017-05-08 [1] CRAN (R 3.5.0)                   
#bitops                 1.0-6      2013-08-17 [1] CRAN (R 3.5.0)                   
#blob                   1.1.1      2018-03-25 [1] CRAN (R 3.5.0)                   
#boot                   1.3-20     2017-08-06 [1] CRAN (R 3.5.1)                   
#callr                  3.1.0      2018-12-10 [1] CRAN (R 3.5.0)                   
#car                    3.0-2      2018-08-23 [1] CRAN (R 3.5.0)                   
#carData                3.0-2      2018-09-30 [1] CRAN (R 3.5.0)                   
#caTools                1.17.1.1   2018-07-20 [1] CRAN (R 3.5.0)                   
#cellranger             1.1.0      2016-07-27 [1] CRAN (R 3.5.0)                   
#circlize               0.4.5      2018-11-21 [1] CRAN (R 3.5.0)                   
#class                  7.3-14     2015-08-30 [1] CRAN (R 3.5.1)                   
#cli                    1.0.1      2018-09-25 [1] CRAN (R 3.5.0)                   
#clue                   0.3-57     2019-02-25 [1] CRAN (R 3.5.2)                   
#cluster                2.0.7-1    2018-04-13 [1] CRAN (R 3.5.1)                   
#codetools              0.2-15     2016-10-05 [1] CRAN (R 3.5.1)                   
#colorspace             1.3-2      2016-12-14 [1] CRAN (R 3.5.0)                   
#ComplexHeatmap       * 2.0.0      2019-05-02 [1] Bioconductor                     
#cowplot                0.9.3      2018-07-15 [1] CRAN (R 3.5.0)                   
#crayon                 1.3.4      2017-09-16 [1] CRAN (R 3.5.0)                   
#curl                   3.2        2018-03-28 [1] CRAN (R 3.5.0)                   
#data.table             1.11.8     2018-09-30 [1] CRAN (R 3.5.0)                   
#DBI                    1.0.0      2018-05-02 [1] CRAN (R 3.5.0)                   
#DelayedArray           0.8.0      2018-10-30 [1] Bioconductor                     
#DEoptimR               1.0-8      2016-11-19 [1] CRAN (R 3.5.0)                   
#desc                   1.2.0      2018-05-01 [1] CRAN (R 3.5.0)                   
#destiny              * 2.12.0     2018-10-30 [1] Bioconductor                     
#devtools               2.0.1      2018-10-26 [1] CRAN (R 3.5.1)                   
#digest                 0.6.18     2018-10-10 [1] CRAN (R 3.5.0)                   
#dplyr                * 0.7.8      2018-11-10 [1] CRAN (R 3.5.0)                   
#e1071                  1.7-0      2018-07-28 [1] CRAN (R 3.5.0)                   
#fitdistrplus           1.0-11     2018-09-10 [1] CRAN (R 3.5.0)                   
#forcats                0.3.0      2018-02-19 [1] CRAN (R 3.5.0)                   
#foreign                0.8-71     2018-07-20 [1] CRAN (R 3.5.0)                   
#fs                     1.2.6      2018-08-23 [1] CRAN (R 3.5.0)                   
#future                 1.10.0     2018-10-17 [1] CRAN (R 3.5.0)                   
#future.apply           1.0.1      2018-08-26 [1] CRAN (R 3.5.0)                   
#gbRd                   0.4-11     2012-10-01 [1] CRAN (R 3.5.0)                   
#gdata                  2.18.0     2017-06-06 [1] CRAN (R 3.5.0)                   
#geneplotter            1.60.0     2018-10-30 [1] Bioconductor                     
#GenomeInfoDb           1.18.1     2018-11-12 [1] Bioconductor                     
#GenomeInfoDbData       1.2.0      2018-11-02 [1] Bioconductor                     
#GenomicRanges          1.34.0     2018-10-30 [1] Bioconductor                     
#GetoptLong             0.1.7      2018-06-10 [1] CRAN (R 3.5.0)                   
#ggplot2              * 3.1.0      2018-10-25 [1] CRAN (R 3.5.0)                   
#ggrepel                0.8.0      2018-05-09 [1] CRAN (R 3.5.0)                   
#ggridges               0.5.1      2018-09-27 [1] CRAN (R 3.5.0)                   
#ggthemes               4.0.1      2018-08-24 [1] CRAN (R 3.5.0)                   
#GlobalOptions          0.1.0      2018-06-09 [1] CRAN (R 3.5.0)                   
#globals                0.12.4     2018-10-11 [1] CRAN (R 3.5.0)                   
#glue                   1.3.0      2018-07-17 [1] CRAN (R 3.5.0)                   
#gplots                 3.0.1      2016-03-30 [1] CRAN (R 3.5.0)                   
#graph                * 1.60.0     2018-10-30 [1] Bioconductor                     
#gridExtra            * 2.3        2017-09-09 [1] CRAN (R 3.5.0)                   
#GSEABase             * 1.44.0     2018-10-30 [1] Bioconductor                     
#GSVA                   1.30.0     2018-10-30 [1] Bioconductor                     
#gtable                 0.2.0      2016-02-26 [1] CRAN (R 3.5.0)                   
#gtools                 3.8.1      2018-06-26 [1] CRAN (R 3.5.0)                   
#haven                  2.0.0      2018-11-22 [1] CRAN (R 3.5.0)                   
#hms                    0.4.2      2018-03-10 [1] CRAN (R 3.5.0)                   
#htmltools              0.3.6      2017-04-28 [1] CRAN (R 3.5.0)                   
#htmlwidgets            1.3        2018-09-30 [1] CRAN (R 3.5.0)                   
#httpuv                 1.4.5      2018-07-19 [1] CRAN (R 3.5.0)                   
#httr                   1.4.0      2018-12-11 [1] CRAN (R 3.5.1)                   
#ica                    1.0-2      2018-05-24 [1] CRAN (R 3.5.0)                   
#igraph                 1.2.2      2018-07-27 [1] CRAN (R 3.5.0)                   
#IRanges              * 2.16.0     2018-10-30 [1] Bioconductor                     
#irlba                  2.3.2      2018-01-11 [1] CRAN (R 3.5.0)                   
#jsonlite               1.6        2018-12-07 [1] CRAN (R 3.5.0)                   
#KernSmooth             2.23-15    2015-06-29 [1] CRAN (R 3.5.1)                   
#laeken                 0.4.6      2014-08-19 [1] CRAN (R 3.5.0)                   
#later                  0.7.5      2018-09-18 [1] CRAN (R 3.5.0)                   
#lattice                0.20-38    2018-11-04 [1] CRAN (R 3.5.0)                   
#lazyeval               0.2.1      2017-10-29 [1] CRAN (R 3.5.0)                   
#listenv                0.7.0      2018-01-21 [1] CRAN (R 3.5.0)                   
#lmtest                 0.9-36     2018-04-04 [1] CRAN (R 3.5.0)                   
#lsei                   1.2-0      2017-10-23 [1] CRAN (R 3.5.0)                   
#magrittr               1.5        2014-11-22 [1] CRAN (R 3.5.0)                   
#MASS                   7.3-51.1   2018-11-01 [1] CRAN (R 3.5.1)                   
#Matrix                 1.2-15     2018-11-01 [1] CRAN (R 3.5.1)                   
#matrixStats            0.54.0     2018-07-23 [1] CRAN (R 3.5.0)                   
#memoise                1.1.0      2017-04-21 [1] CRAN (R 3.5.0)                   
#metap                  1.0        2018-07-25 [1] CRAN (R 3.5.0)                   
#mime                   0.6        2018-10-05 [1] CRAN (R 3.5.0)                   
#munsell                0.5.0      2018-06-12 [1] CRAN (R 3.5.0)                   
#nnet                   7.3-12     2016-02-02 [1] CRAN (R 3.5.1)                   
#npsurv                 0.4-0      2017-10-14 [1] CRAN (R 3.5.0)                   
#openxlsx             * 4.1.0      2018-05-26 [1] CRAN (R 3.5.0)                   
#outliers               0.14       2011-01-24 [1] CRAN (R 3.5.0)                   
#pbapply                1.3-4      2018-01-10 [1] CRAN (R 3.5.0)                   
#pbmcapply              1.3.0      2018-11-05 [1] CRAN (R 3.5.0)                   
#pheatmap               1.0.12     2019-01-04 [1] CRAN (R 3.5.2)                   
#pillar                 1.3.0      2018-07-14 [1] CRAN (R 3.5.0)                   
#pkgbuild               1.0.2      2018-10-16 [1] CRAN (R 3.5.0)                   
#pkgconfig              2.0.2      2018-08-16 [1] CRAN (R 3.5.0)                   
#pkgload                1.0.2      2018-10-29 [1] CRAN (R 3.5.0)                   
#plotly                 4.8.0      2018-07-20 [1] CRAN (R 3.5.1)                   
#plyr                   1.8.4      2016-06-08 [1] CRAN (R 3.5.0)                   
#png                    0.1-7      2013-12-03 [1] CRAN (R 3.5.0)                   
#prettyunits            1.0.2      2015-07-13 [1] CRAN (R 3.5.0)                   
#processx               3.2.1      2018-12-05 [1] CRAN (R 3.5.1)                   
#promises               1.0.1      2018-04-13 [1] CRAN (R 3.5.0)                   
#proxy                  0.4-22     2018-04-08 [1] CRAN (R 3.5.0)                   
#ps                     1.2.1      2018-11-06 [1] CRAN (R 3.5.1)                   
#purrr                  0.2.5      2018-05-29 [1] CRAN (R 3.5.0)                   
#R.methodsS3            1.7.1      2016-02-16 [1] CRAN (R 3.5.0)                   
#R.oo                   1.22.0     2018-04-22 [1] CRAN (R 3.5.0)                   
#R.utils                2.7.0      2018-08-27 [1] CRAN (R 3.5.0)                   
#R6                     2.3.0      2018-10-04 [1] CRAN (R 3.5.0)                   
#RANN                   2.6        2018-07-16 [1] CRAN (R 3.5.0)                   
#RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 3.5.0)                   
#Rcpp                   1.0.0      2018-11-07 [1] CRAN (R 3.5.0)                   
#RcppEigen              0.3.3.5.0  2018-11-24 [1] CRAN (R 3.5.1)                   
#RCurl                  1.95-4.11  2018-07-15 [1] CRAN (R 3.5.0)                   
#Rdpack                 0.10-1     2018-10-04 [1] CRAN (R 3.5.0)                   
#readxl                 1.1.0      2018-04-20 [1] CRAN (R 3.5.0)                   
#remotes                2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
#reshape2               1.4.3      2017-12-11 [1] CRAN (R 3.5.0)                   
#reticulate             1.10       2018-08-05 [1] CRAN (R 3.5.0)                   
#rio                    0.5.16     2018-11-26 [1] CRAN (R 3.5.0)                   
#rjson                  0.2.20     2018-06-08 [1] CRAN (R 3.5.0)                   
#rlang                  0.4.0      2019-06-25 [1] CRAN (R 3.5.2)                   
#robustbase             0.93-3     2018-09-21 [1] CRAN (R 3.5.0)                   
#ROCR                   1.0-7      2015-03-26 [1] CRAN (R 3.5.0)                   
#rprojroot              1.3-2      2018-01-03 [1] CRAN (R 3.5.0)                   
#RSQLite                2.1.1      2018-05-06 [1] CRAN (R 3.5.0)                   
#rstudioapi             0.8        2018-10-02 [1] CRAN (R 3.5.0)                   
#rsvd                   1.0.0      2018-11-06 [1] CRAN (R 3.5.0)                   
#Rtsne                  0.15       2018-11-10 [1] CRAN (R 3.5.0)                   
#S4Vectors            * 0.20.1     2018-11-09 [1] Bioconductor                     
#scales                 1.0.0      2018-08-09 [1] CRAN (R 3.5.0)                   
#scatterplot3d          0.3-41     2018-03-14 [1] CRAN (R 3.5.0)                   
#SDMTools               1.1-221    2014-08-05 [1] CRAN (R 3.5.0)                   
#sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 3.5.1)                   
#Seurat               * 3.0.0.9000 2018-12-12 [1] Github (satijalab/seurat@aa9ffda)
#shape                  1.4.4      2018-02-07 [1] CRAN (R 3.5.0)                   
#shiny                  1.2.0      2018-11-02 [1] CRAN (R 3.5.0)                   
#shinythemes            1.1.2      2018-11-06 [1] CRAN (R 3.5.0)                   
#SingleR              * 0.2.1      2018-12-12 [1] local                            
#smoother               1.1        2015-04-16 [1] CRAN (R 3.5.0)                   
#sp                     1.3-1      2018-06-05 [1] CRAN (R 3.5.0)                   
#stringi                1.2.4      2018-07-20 [1] CRAN (R 3.5.0)                   
#stringr                1.3.1      2018-05-10 [1] CRAN (R 3.5.0)                   
#SummarizedExperiment   1.12.0     2018-10-30 [1] Bioconductor                     
#survival               2.43-3     2018-11-26 [1] CRAN (R 3.5.0)                   
#tibble                 1.4.2      2018-01-22 [1] CRAN (R 3.5.0)                   
#tidyr                  0.8.2      2018-10-28 [1] CRAN (R 3.5.0)                   
#tidyselect             0.2.5      2018-10-11 [1] CRAN (R 3.5.0)                   
#tsne                   0.1-3      2016-07-15 [1] CRAN (R 3.5.0)                   
#TTR                    0.23-4     2018-09-20 [1] CRAN (R 3.5.0)                   
#usethis                1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
#vcd                    1.4-4      2017-12-06 [1] CRAN (R 3.5.0)                   
#VIM                    4.7.0      2017-04-11 [1] CRAN (R 3.5.0)                   
#viridisLite            0.3.0      2018-02-01 [1] CRAN (R 3.5.0)                   
#withr                  2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
#XML                  * 3.98-1.16  2018-08-19 [1] CRAN (R 3.5.0)                   
#xtable                 1.8-3      2018-08-29 [1] CRAN (R 3.5.0)                   
#xts                    0.11-2     2018-11-05 [1] CRAN (R 3.5.0)                   
#XVector                0.22.0     2018-10-30 [1] Bioconductor                     
#yaml                   2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
#zip                    1.0.0      2017-04-25 [1] CRAN (R 3.5.0)                   
#zlibbioc               1.28.0     2018-10-30 [1] Bioconductor                     
#zoo                    1.8-4      2018-09-19 [1] CRAN (R 3.5.0)                   
#
#[1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library