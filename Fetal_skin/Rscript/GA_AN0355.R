#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a healthy donor
# Date: 2020.03.13
# Adapted from https://github.com/10XGenomics/single-cell-3prime-snp-clustering
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger with genotyping through SNV using snpclust
# - genotype cells from the skin and merged blood and skin (after removing overlaping cell barcodes from skin) 
#   of the HSCT patient (host/donor)
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(dplyr)
#library(SingleR)
#library(grid)
#library(pheatmap)
#library(RColorBrewer)
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
                      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = "black"), 
                      axis.ticks.y = element_line(size = 0.4, colour = "black"),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour="black"),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = "black", size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 14, face = "bold", margin = margin(0,0,5,0)), 
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
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])

# Define samples analyzed (all)
#sa <- c("cd")

# Define reference sample name
ct <- c("blood")
tr <- c("skin")

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c("sc_skin_blood_genotyping_new")
an.descs <- c("sc_genotyping")

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx('../../ngs_sample_list.xlsm', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
#nsl <- nsl[nsl$Experiment_ID == experimentid1 | nsl$Experiment_ID == experimentid2, ]
nsl <- nsl[nsl$Experiment_ID == experimentid, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], plyr::summarise, sum = NA)

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
dir.create(paste0("figures/", analysisid, "_", an.desc, "/barcodes"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/tsne"))
dir.figs <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
dir.figs.bc <- paste0("figures/", analysisid, "_", an.desc, "/barcodes/", analysisid, "_", an.descs)
#dir.figs.ts <- paste0("figures/", analysisid, "_", an.desc, "/tsne/", analysisid, "_", an.descs)
#dir.figs.um <- paste0("figures/", analysisid, "_", an.desc, "/umap/", analysisid, "_", an.descs)
#dir.figs.hm <- paste0("figures/", analysisid, "_", an.desc, "/heatmap/", analysisid, "_", an.descs)
#dir.figs.vp <- paste0("figures/", analysisid, "_", an.desc, "/violin/", analysisid, "_", an.descs)

#-----------------------------------------------------------------------------------------------------------------#
# Identify host and donor cells
#-----------------------------------------------------------------------------------------------------------------#

# Read in summary json file
json.sk<-jsonlite::fromJSON("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/summary.json")
#json.bd<-jsonlite::fromJSON("../../procdata/EX0008/GA_PD0008_snpclust/tp_blood_1/summary.json")
json.bd<-jsonlite::fromJSON("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/summary.json") # parameters 80_2_1_1
#json.mg<-jsonlite::fromJSON("../../procdata/EX0008/tpMerge/outs/summary.json")
json.mg<-jsonlite::fromJSON("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/summary.json")

# Load barcodes for skin and blood
bc.sk <- readLines("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/raw_allele_bc_matrices_mex/ref/barcodes.tsv")
#bc.bd <- readLines("../../procdata/EX0008/GA_PD0008_snpclust/tp_blood_1/raw_allele_bc_matrices_mex/ref/barcodes.tsv")
bc.bd <- readLines("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/raw_allele_bc_matrices_mex/ref/barcodes.tsv")
#bc.mg <- readLines("../../procdata/EX0008/tpMerge/outs/raw_allele_bc_matrices_mex/ref/barcodes.tsv")
bc.mg <- readLines("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/raw_allele_bc_matrices_mex/ref/barcodes.tsv")

# Mix should have all barcodes of skin (17 barcodes overlapped with blood and were removed from blood for snpclust analysis)
#length(intersect(bc.mg,bc.sk))==length(bc.sk)
intersect(bc.bd,bc.sk)
#length(bc.sk)+length(bc.bd)==length(bc.mg)
length(bc.sk)+length(bc.bd)

# Extract genotype calls
#table(bc.sk %in% read.delim('../201902/figures/GA_AN0162_sc_skin_blood_genotyping/GA_AN0162_sc_skin_genotyping_barcodes_host.tsv',as.is = T,header = F)[,1])
gt.sk <- setNames(ifelse(json.sk$model2_call==0, 'host', 'donor'), nm = bc.sk)
gt.bd <- setNames(ifelse(json.bd$model2_call==0, 'donor', 'host'), nm = bc.bd) #INVERT!!!
gt.mg <- setNames(ifelse(json.mg$model2_call==0, 'host', 'donor'), nm = bc.mg) #INVERT!!!
#tb.gt <- data.frame(cbind(skin=rev(table(gt.sk)), blood=rev(table(gt.bd)), merged=rev(table(gt.mg))))
tb.gt <- data.frame(cbind(skin=rev(table(gt.sk)), merged=rev(table(gt.mg))))
#gt.mg <- setNames(ifelse(json.mg$model2_call==0, 'host', 'donor'), nm = bc.mg)
#gt.mg[names(gt.mg) %in% intersect(bc.bd,bc.sk)] #10 cells from skin were wrongly called as host in the merged data
table(gt.mg[names(gt.mg) %in% bc.bd]); table(gt.mg[names(gt.mg) %in% bc.sk])
#intersect(names(gt.mg2[gt.mg2 %in% 'host']),names(gt.mg[gt.mg %in% 'donor'])) %in% bc.bd '2 for the old analysis'

## Filter skin cells from merged data with blood
#json.mg <- lapply(setNames(names(json.mg), nm = names(json.mg)), function(x) {
#  json.mg[[x]][json.mg[[x]]] })

# Define barcodes
bc.gt.sk <- data.frame(genotype=gt.sk, barcode=names(gt.sk), row.names = 1:length(gt.sk))[,2:1]
bc.gt.bd <- data.frame(genotype=gt.bd, barcode=names(gt.bd), row.names = 1:length(gt.bd))[,2:1]
bc.gt.mg <- data.frame(genotype=gt.mg, barcode=names(gt.mg), row.names = 1:length(gt.mg))[,2:1]

#-----------------------------------------------------------------------------------------------------------------#
# Identify probability of sample having 2 genotypes
#-----------------------------------------------------------------------------------------------------------------#

# Decide which model to use (use model 2, for 2 genotypes, if at least 90% of cells fit)
mp.sk <- sum(json.sk$model2_probability>0.70) / sum(json.sk$model1_probability>0.70)
mp.bd <- sum(json.bd$model2_probability>0.70) / sum(json.bd$model1_probability>0.70)
mp.mg <- sum(json.mg$model2_probability>0.70) / sum(json.mg$model1_probability>0.70)
#mp.mg <- sum(json.mg$model2_probability>0.75) / sum(json.mg$model1_probability>0.75)
tb.mp <- data.frame(skin = paste0(round(mp.sk, digits = 3)*100, "%"), 
                    #blood = paste0(round(mp.bd, digits = 3)*100, "%"), 
                    merged = paste0(round(mp.mg, digits = 3)*100, "%"), row.names = "Probability of 2 genotypes")

#-----------------------------------------------------------------------------------------------------------------#
# Identify amount of cells for each genotype
#-----------------------------------------------------------------------------------------------------------------#

# Number of cells for each genotype
nc.sk <- table(json.sk$model2_call)
table(json.sk$model2_g1_call)/length(json.sk$model2_g1_call)
table(json.sk$model2_g0_call)/length(json.sk$model2_g0_call)
nc.bd <- table(json.bd$model2_call)
table(json.bd$model2_g1_call)/length(json.bd$model2_g1_call)
table(json.bd$model2_g0_call)/length(json.bd$model2_g0_call)
nc.mg <- table(json.mg$model2_call)
table(json.mg$model2_g1_call)/length(json.mg$model2_g1_call)
table(json.mg$model2_g0_call)/length(json.mg$model2_g0_call)
#nc <- data.frame(skin = as.vector(nc.sk), blood = as.vector(nc.bd), row.names = c("genotype 1", "genotype 2"))
#tb.nc.sk <- data.frame(Frequency = rev(as.vector(nc.sk)), row.names = c("genotype S1", "genotype S2"))
#tb.nc.bd <- data.frame(Frequency = rev(as.vector(nc.bd)), row.names = c("genotype B1", "genotype B2"))
#tb.nc <- rbind(tb.nc.sk, tb.nc.bd)
#tb.nc <- data.frame(cbind(skin=nc.sk, blood=rev(nc.bd), merged=nc.mg), row.names = paste('genotype',1:2)) #INVERT!!!
#tb.nc <- data.frame(cbind(skin=nc.sk, blood=rev(nc.bd), merged=nc.mg), row.names = c('host','donor')) #INVERT!!!
tb.nc <- data.frame(cbind(skin=nc.sk, merged=nc.mg), row.names = paste('genotype',1:2))

#-----------------------------------------------------------------------------------------------------------------#
# Frequency of SNV detection
#-----------------------------------------------------------------------------------------------------------------#

# Read in reference and alt allele counts at each snv (row) in each cell (column)
ref.sk<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/raw_allele_bc_matrices_mex/ref/matrix.mtx"))
alt.sk<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/raw_allele_bc_matrices_mex/alt/matrix.mtx"))
#ref_rs<-rowSums(ref)
#alt_rs<-rowSums(alt)
snps_count.sk<-rowSums(ref.sk+alt.sk)
snps_freq.sk<-alt.sk/(alt.sk+ref.sk+1)
median(snps_count.sk)
median(colSums(ref.sk>0))
median(colSums(alt.sk>0))
hist(colSums(alt.sk>0))
median(colSums(alt.sk+ref.sk>0))
ref.bd<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_blood_1/raw_allele_bc_matrices_mex/ref/matrix.mtx"))
alt.bd<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_blood_1/raw_allele_bc_matrices_mex/alt/matrix.mtx"))
#ref_rs<-rowSums(ref)
#alt_rs<-rowSums(alt)
snps_count.bd<-rowSums(ref.bd+alt.bd)
snps_freq.bd<-alt.bd/(alt.bd+ref.bd+1)
median(snps_count.bd)
median(colSums(ref.bd>0))
median(colSums(alt.bd>0))
hist(colSums(alt.bd>0))
median(colSums(alt.bd+ref.bd>0))
ref.mg<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/raw_allele_bc_matrices_mex/ref/matrix.mtx"))
alt.mg<-as.matrix(Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/raw_allele_bc_matrices_mex/alt/matrix.mtx"))
#ref_rs<-rowSums(ref)
#alt_rs<-rowSums(alt)
snps_count.mg<-rowSums(ref.mg+alt.mg)
snps_freq.mg<-alt.mg/(alt.mg+ref.mg+1)
median(snps_count.mg)
median(colSums(ref.mg>0))
median(colSums(alt.mg>0))
hist(colSums(alt.mg>0))
median(colSums(alt.mg+ref.mg>0))

# Plot frequency of cells where distinct SNV were detected
hist.cs.sk <- ggplot(data.frame(value = rowSums(alt.sk+ref.sk>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,25,1), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed") +
  #scale_x_continuous("Number of cells containing a SNV", sec.axis = dup_axis(name = "", breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of cells per SNV", breaks = c(seq(0,25,5), median(rowSums(alt.sk+ref.sk>0)))) + 
  scale_y_continuous("Number of SNV") + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT skin', subtitle = "Frequency of cells\ncontaining a distinct SNV") + theme_custom + theme(aspect.ratio = 0.7)
hist.cs.sk
hist.cs.bd <- ggplot(data.frame(value = rowSums(alt.bd+ref.bd>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,25,1), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed") +
  #scale_x_continuous("Number of cells containing a SNV", sec.axis = dup_axis(name = "", breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of cells per SNV", breaks = c(seq(0,25,5), median(rowSums(alt.bd+ref.bd>0)))) + 
  scale_y_continuous("Number of SNV") + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT blood', subtitle = "Frequency of cells\ncontaining a distinct SNV") + theme_custom + theme(aspect.ratio = 0.7)
hist.cs.bd
hist.cs.mg <- ggplot(data.frame(value = rowSums(alt.mg+ref.mg>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,25,1), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed") +
  #scale_x_continuous("Number of cells containing a SNV", sec.axis = dup_axis(name = "", breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of cells per SNV", breaks = c(seq(0,25,5), median(rowSums(alt.mg+ref.mg>0)))) + 
  scale_y_continuous("Number of SNV") + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT merged', subtitle = "Frequency of cells\ncontaining a distinct SNV") + theme_custom + theme(aspect.ratio = 0.7)
hist.cs.mg

# Plot frequency of SNV detected per cell
hist.sc.sk <- ggplot(data.frame(value = colSums(alt.sk+ref.sk>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,240,10), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed")+
  #scale_x_continuous("Number of SNV per cell", sec.axis = dup_axis(name = "", breaks = median(colSums(alt.sk>0)), labels = round(median(colSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of SNV per cell", breaks = c(seq(0,240,60), median(colSums(alt.sk+ref.sk>0)))) + 
  scale_y_continuous("Number of cells", limits=c(0,300)) + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT skin', "Frequency of SNV\ndetected per cell") + theme_custom + theme(aspect.ratio = 0.7)
hist.sc.sk
hist.sc.bd <- ggplot(data.frame(value = colSums(alt.bd+ref.bd>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,240,10), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed")+
  #scale_x_continuous("Number of SNV per cell", sec.axis = dup_axis(name = "", breaks = median(colSums(alt.sk>0)), labels = round(median(colSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of SNV per cell", breaks = c(seq(0,240,60)[-3], median(colSums(alt.bd+ref.bd>0)))) + 
  scale_y_continuous("Number of cells", limits=c(0,1100), breaks = seq(0,1100,250)) + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT blood', "Frequency of SNV\ndetected per cell") + theme_custom + theme(aspect.ratio = 0.7)
hist.sc.bd
hist.sc.mg <- ggplot(data.frame(value = colSums(alt.mg+ref.mg>0)), aes(value)) + 
  geom_histogram(color = "black", breaks = seq(0,240,10), position = "identity", alpha = 0.3, size = 0.2) + 
  geom_vline(aes(xintercept = median(value)), linetype = "dashed")+
  #scale_x_continuous("Number of SNV per cell", sec.axis = dup_axis(name = "", breaks = median(colSums(alt.sk>0)), labels = round(median(colSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
  scale_x_continuous("Number of SNV per cell", breaks = c(seq(0,250,60), median(colSums(alt.mg+ref.mg>0)))) + 
  scale_y_continuous("Number of cells", limits=c(0,1100), breaks = seq(0,1100,250)) + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
  ggtitle('HSCT merged', "Frequency of SNV\ndetected per cell") + theme_custom + theme(aspect.ratio = 0.7)
hist.sc.mg

#hist.cs.sk <- ggplot(data.frame(value = rowSums(alt.sk+ref.sk>0)), aes(value)) + 
#  geom_histogram(color = "black", breaks = seq(0,100,2), position = "identity", alpha = 0.3, size = 0.2) + 
#  geom_vline(aes(xintercept = median(value)), linetype = "dashed") +
#  #scale_x_continuous("Number of cells containing a SNV", sec.axis = dup_axis(name = "", breaks = median(rowSums(alt.sk+ref.sk>0)), labels = round(median(rowSums(alt.sk+ref.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#  scale_x_continuous("Number of cells containing the SNV", breaks = c(seq(0,100,20), median(rowSums(alt.sk+ref.sk>0))), labels = c(seq(0,100,20), round(median(rowSums(alt.sk+ref.sk>0))))) + 
#  scale_y_continuous("Number of distinct SNV") + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
#  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
#  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
#  ggtitle("Frequency of cells\ncontaining a distinct SNV") + theme_custom + theme(aspect.ratio = NULL)
#hist.cs.sk
#
## Plot frequency of SNV detected per cell
#hist.sc.sk <- ggplot(data.frame(value = colSums(alt.sk+ref.sk>0)), aes(value)) + 
#  geom_histogram(color = "black", breaks = seq(0,250,10), position = "identity", alpha = 0.3, size = 0.2) + 
#  geom_vline(aes(xintercept = median(value)), linetype = "dashed") +
#  #scale_x_continuous("Number of SNV per cell", sec.axis = dup_axis(name = "", breaks = median(colSums(alt.sk+ref.sk>0)), labels = round(median(colSums(alt.sk+ref.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
#  scale_x_continuous("Number of SNV per cell", breaks = c(seq(0,250,50), median(colSums(alt.sk+ref.sk>0))), labels = c(seq(0,250,50), round(median(colSums(alt.sk+ref.sk>0))))) + 
#  scale_y_continuous("Number of cells") + #, limits = c(0,4000), breaks = seq(0,4000,1000), expand = c(0,0)) + 
#  #scale_fill_manual("p-value:", values = c("black", "firebrick3")) + 
#  #scale_color_manual("p-value:", values = c("black", "firebrick3")) + 
#  ggtitle("Frequency of SNV\ndetected per cell") + theme_custom + theme(aspect.ratio = NULL)
#hist.sc.sk

#-----------------------------------------------------------------------------------------------------------------#
# Number of allele counts
#-----------------------------------------------------------------------------------------------------------------#

# Read in genotype likelihood at each snv in each cell
rr.sk<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/likelihood_allele_bc_matrices_mex/0|0/matrix.mtx")
ra.sk<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/likelihood_allele_bc_matrices_mex/1|0/matrix.mtx")
aa.sk<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/likelihood_allele_bc_matrices_mex/1|1/matrix.mtx")
rr.bd<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/likelihood_allele_bc_matrices_mex/0|0/matrix.mtx")
ra.bd<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/likelihood_allele_bc_matrices_mex/1|0/matrix.mtx")
aa.bd<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/likelihood_allele_bc_matrices_mex/1|1/matrix.mtx")
#rr.mg<-Matrix::readMM("../../procdata/EX0008/tpMerge/outs/likelihood_allele_bc_matrices_mex/0|0/matrix.mtx")
#ra.mg<-Matrix::readMM("../../procdata/EX0008/tpMerge/outs/likelihood_allele_bc_matrices_mex/1|0/matrix.mtx")
#aa.mg<-Matrix::readMM("../../procdata/EX0008/tpMerge/outs/likelihood_allele_bc_matrices_mex/1|1/matrix.mtx")
rr.mg<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/likelihood_allele_bc_matrices_mex/0|0/matrix.mtx")
ra.mg<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/likelihood_allele_bc_matrices_mex/1|0/matrix.mtx")
aa.mg<-Matrix::readMM("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/likelihood_allele_bc_matrices_mex/1|1/matrix.mtx")

# Calculate the inferred genotype at each snv in each cell
rr_rs.sk<-rowSums(as.matrix(rr.sk))
ra_rs.sk<-rowSums(as.matrix(ra.sk))
aa_rs.sk<-rowSums(as.matrix(aa.sk))
ga.sk<-apply(cbind(rr_rs.sk,ra_rs.sk,aa_rs.sk),1,which.max)
table(ga.sk)
rr_rs.bd<-rowSums(as.matrix(rr.bd))
ra_rs.bd<-rowSums(as.matrix(ra.bd))
aa_rs.bd<-rowSums(as.matrix(aa.bd))
ga.bd<-apply(cbind(rr_rs.bd,ra_rs.bd,aa_rs.bd),1,which.max)
table(ga.bd)
rr_rs.mg<-rowSums(as.matrix(rr.mg))
ra_rs.mg<-rowSums(as.matrix(ra.mg))
aa_rs.mg<-rowSums(as.matrix(aa.mg))
ga.mg<-apply(cbind(rr_rs.mg,ra_rs.mg,aa_rs.mg),1,which.max)
table(ga.mg)
#tb.gs <- data.frame(Frequency = as.vector(table(gt)), row.names = c("ref/ref", "ref/alt", "alt/alt"))
#tb.gs <- data.frame(cbind(skin=table(ga.sk), blood=table(ga.bd), merged=table(ga.mg)), row.names = c("ref/ref", "ref/alt", "alt/alt"))
tb.gs <- data.frame(cbind(skin=table(ga.sk), merged=table(ga.mg)), row.names = c("ref/ref", "ref/alt", "alt/alt"))

#-----------------------------------------------------------------------------------------------------------------#
# Compare genotypes obtained to reference blood sample (should be equal to genotype from donor skin cells)
#-----------------------------------------------------------------------------------------------------------------#

# Load identified SNV
snp.sk <- read.delim("../../procdata/EX0008/GA_PD0008_snpclust/tp_skin_1/raw_allele_bc_matrices_mex/alt/genes.tsv",as.is = T,header = F)[,1]
#snp.bd <- read.delim("../../procdata/EX0008/GA_PD0008_snpclust/tp_blood_1/raw_allele_bc_matrices_mex/alt/genes.tsv",as.is = T,header = F)[,1]
snp.bd <- read.delim("../../procdata/EX0008/GA_PD0008_snpclust_80_2_1_1_1e-3/raw_allele_bc_matrices_mex/alt/genes.tsv",as.is = T,header = F)[,1]
#snp.mg <- read.delim("../../procdata/EX0008/tpMerge/outs/raw_allele_bc_matrices_mex/alt/genes.tsv",as.is = T,header = F)[,1]
snp.mg <- read.delim("../../procdata/EX0008/GA_PD0008_snpclust_merged/tp_merged_1/raw_allele_bc_matrices_mex/alt/genes.tsv",as.is = T,header = F)[,1]
length(intersect(snp.sk,snp.bd))
length(intersect(snp.sk,snp.mg))

#gen.bd=data.frame(json.bd$model1_g0_call,row.names = snp.bd)
#gen.sk1=data.frame(json.sk$model2_g0_call,row.names = snp.sk)
#gen.sk2=data.frame(json.sk$model2_g1_call,row.names = snp.sk)
#gen=merge(gen.bd,gen.sk1,by="row.names")
#gen=merge(gen,gen.sk2,by="row.names")

# Merge data and filter for common SNV
#gen.bd <- data.frame(snv=snp.bd, blood=json.bd$model1_g0_call)
gen.sk1 <- data.frame(snv=snp.sk, skin1=json.sk$model2_g0_call)
gen.sk2 <- data.frame(snv=snp.sk, skin2=json.sk$model2_g1_call)
gen.bd1 <- data.frame(snv=snp.bd, blood1=json.bd$model2_g1_call)
gen.bd2 <- data.frame(snv=snp.bd, blood2=json.bd$model2_g0_call)
gen.mg1 <- data.frame(snv=snp.mg, merged1=json.mg$model2_g0_call)
gen.mg2 <- data.frame(snv=snp.mg, merged2=json.mg$model2_g1_call)
#common.snv <- merge(data.frame(snv=snp.sk), data.frame(snv=snp.bd), by="snv")$snv
#gen <- merge(gen.bd, gen.sk1, by="snv")
#gen <- merge(gen, gen.sk2, by="snv")
gen.sk <- merge(gen.sk1, gen.sk2, by="snv")
gen.bd <- merge(gen.bd1, gen.bd2, by="snv")
gen <- merge(gen.sk, gen.bd, by="snv")
gen.mg <- merge(gen.mg1, gen.mg2, by="snv")
gen2 <- merge(gen.sk, gen.mg, by="snv")
gen3 <- merge(gen.bd, gen.mg, by="snv")
#gen3 <- merge(gen3, gen.sk, by="snv")

# Check similarity between SNV identified in blood sample and in each genotype from skin sample
table(gen.sk1$skin1 == gen.sk2$skin2)
table(gen['blood1'] == gen['blood2'])
sum(gen.bd1$blood1 == gen.bd2$blood2)/dim(gen.bd1)[1]
sim.s1s2 <- sum(gen['skin1'] == gen['skin2'])/dim(gen)[1]
sim.s1b1 <- sum(gen['skin1'] == gen['blood1'])/dim(gen)[1]
sim.s1b2 <- sum(gen['skin1'] == gen['blood2'])/dim(gen)[1]
#sim.s2s2 <- sum(gen['skin2'] == gen['skin1'])/dim(gen)[1]
sim.s2b1 <- sum(gen['skin2'] == gen['blood1'])/dim(gen)[1]
sim.s2b2 <- sum(gen['skin2'] == gen['blood2'])/dim(gen)[1]
sim.b1b2 <- sum(gen['blood1'] == gen['blood2'])/dim(gen)[1]
sim.s1m1 <- sum(gen2['skin1'] == gen2['merged1'])/dim(gen2)[1]
sim.s1m2 <- sum(gen2['skin1'] == gen2['merged2'])/dim(gen2)[1]
sim.s2m1 <- sum(gen2['skin2'] == gen2['merged1'])/dim(gen2)[1]
sim.s2m2 <- sum(gen2['skin2'] == gen2['merged2'])/dim(gen2)[1]
sim.s1s22 <- sum(gen2['skin1'] == gen2['skin2'])/dim(gen2)[1]
sim.m1m2 <- sum(gen2['merged1'] == gen2['merged2'])/dim(gen2)[1]
sim.b1m1 <- sum(gen3['blood1'] == gen3['merged1'])/dim(gen3)[1]
sim.b1m2 <- sum(gen3['blood1'] == gen3['merged2'])/dim(gen3)[1]
sim.b2m1 <- sum(gen3['blood2'] == gen3['merged1'])/dim(gen3)[1]
sim.b2m2 <- sum(gen3['blood2'] == gen3['merged2'])/dim(gen3)[1]
sim.b1b22 <- sum(gen3['blood1'] == gen3['blood2'])/dim(gen3)[1]
sim.m1m22 <- sum(gen3['merged1'] == gen3['merged2'])/dim(gen3)[1]

#ss <- data.frame("blood vs skin genotype 1" = sbd.g1, "blood vs skin genotype 2" = sbd.g2, 
#                 "skin genotypes 1 vs 2" = ssk, row.names = "Similarity")
#tb.ss <-data.frame(Similarity = paste0(round(c(sim.s1s2, sim.s1b1, sim.s1b2, sim.s2b1, sim.s2b2), digits = 2)*100, "%"), 
#                   row.names = c('skin g1 vs skin g2','skin g1 vs blood g1','skin g1 vs blood g2',
#                                 'skin g2 vs blood g1','skin g2 vs blood g2'))
#tb.ss <-data.frame(Similarity = paste0(round(c(sim.s1b2,sim.s2b2), digits = 2)*100, "%"), 
#                   row.names = c('skin g1 vs blood donor','skin g2 vs blood donor'))
#tb.ss2 <-data.frame(Similarity = paste0(round(c(sim.s1s22, sim.m1m2, sim.s1m1, sim.s1m2, sim.s2m1, sim.s2m2), digits = 2)*100, "%"), 
#                   row.names = c('skin g1 vs skin g2','merged g1 vs merged g2','skin g1 vs merged g1','skin g1 vs merged g2',
#                                 'skin g2 vs merged g1','skin g2 vs merged g2'))
#tb.ss3 <-data.frame(Similarity = paste0(round(c(sim.b1b22, sim.m1m22, sim.b1m1, sim.b1m2, sim.b2m1, sim.b2m2), digits = 2)*100, "%"), 
#                    row.names = c('blood g1 vs blood g2','merged g1 vs merged g2','blood g1 vs merged g1','blood g1 vs merged g2',
#                                  'blood g2 vs merged g1','blood g2 vs merged g2'))
tb.ss <-data.frame(Similarity = paste0(round(c(sim.s1b2,sim.s2b2), digits = 2)*100, "%"), 
                   row.names = c('skin host vs blood donor','skin donor vs blood donor'))
tb.ss2 <-data.frame(Similarity = paste0(round(c(sim.s1s22, sim.m1m2, sim.s1m1, sim.s1m2, sim.s2m1, sim.s2m2), digits = 2)*100, "%"), 
                    row.names = c('skin host vs skin donor','merged host vs merged donor','skin host vs merged host','skin host vs merged donor',
                                  'skin donor vs merged host','skin donor vs merged donor'))
tb.ss3 <-data.frame(Similarity = paste0(round(c(sim.b1b22, sim.m1m22, sim.b1m1, sim.b1m2, sim.b2m1, sim.b2m2), digits = 2)*100, "%"), 
                    row.names = c('blood host vs blood donor','merged host vs merged donor','blood host vs merged host','blood host vs merged donor',
                                  'blood donor vs merged host','blood donor vs merged donor'))

#-----------------------------------------------------------------------------------------------------------------#
# Skin genotypes similarity to blood for each single cell
#-----------------------------------------------------------------------------------------------------------------#

# Create data table
countLimit <- 1
cellGenotypeMatrix=Matrix::Matrix(0,nrow = nrow(ref.sk), ncol = ncol(ref.sk), sparse = TRUE)   #duplicate(ref.sk)
cellGenotypeMatrix[which(ref.sk>=countLimit)]=1
cellGenotypeMatrix[which(alt.sk>=countLimit)]=3
cellGenotypeMatrix[which((ref.sk>=countLimit) & (alt.sk>=countLimit))]=2
#sanity chekcs
sum(cellGenotypeMatrix==2 | cellGenotypeMatrix==3)==sum(alt.sk>=countLimit)
sum(cellGenotypeMatrix==2 | cellGenotypeMatrix==1)==sum(ref.sk>=countLimit)

# Create functions to compare SNV per cell
#compareSNVToBloodFromMerged <- function(cellGenotype){
#  #genotypeBlood=data.frame(snpNames = snp.bd, blood_Geno = json.bd$model1_g0_call)
#  genotypeBlood=data.frame(snpNames = snp.mg, blood_Geno = json.mg$model2_g1_call)#[bc.mg %in% bc.bd,]#[!(bc.mg %in% bc.sk),]  #Donor
#  genotypeSkin=data.frame(snpNames =snp.sk, skinCellGeno = cellGenotype-1) #-1 for different allele-codes...
#  genoMerge=merge(genotypeBlood, genotypeSkin,by="snpNames")
#  cSame=sum(genoMerge["blood_Geno"]==genoMerge["skinCellGeno"])
#  cAll=sum(genoMerge["skinCellGeno"]>=0)
#  return(cSame/cAll)
#}
table(bc.mg %in% bc.sk)
compareSNVToBloodGenotypesFromMerged <- function(cellGenotype){
  #genotypeBlood=data.frame(snpNames = snp.bd, blood_Geno = json.bd$model1_g0_call)
  #genotypeBlood=data.frame(snpNames = snp.mg, blood_Geno = json.mg$model2_g0_call)[bc.mg %in% bc.bd,]#[!(bc.mg %in% bc.sk),]
  genotypeBlood1=data.frame(snpNames = snp.mg, blood_Geno = json.mg$model2_g0_call)#[json.mg$model2_call==0,]
  genotypeBlood2=data.frame(snpNames = snp.mg, blood_Geno = json.mg$model2_g1_call)#[json.mg$model2_call==1,]
  genotypeSkin=data.frame(snpNames =snp.sk, skinCellGeno = cellGenotype-1) #-1 for different allele-codes...
  genoMerge1=merge(genotypeBlood1, genotypeSkin,by="snpNames")
  genoMerge2=merge(genotypeBlood2, genotypeSkin,by="snpNames")
  cSame1=sum(genoMerge1["blood_Geno"]==genoMerge1["skinCellGeno"])
  cSame2=sum(genoMerge2["blood_Geno"]==genoMerge2["skinCellGeno"])
  cAll1=sum(genoMerge1["skinCellGeno"]>=0)
  cAll2=sum(genoMerge2["skinCellGeno"]>=0)
  return(c(cSame1/cAll1,cSame2/cAll2))
}
compareSNVToSkin <- function(cellGenotype){
  cSame1=sum(cellGenotype==(json.sk$model2_g0_call+1)) #+1 for the same reason as -1 above
  cSame2=sum(cellGenotype==(json.sk$model2_g1_call+1))
  cAll=sum(cellGenotype>=1)
  return(c(cSame1/cAll,cSame2/cAll))
}

# Compare SNV between samples
skinComparison=apply(cellGenotypeMatrix,2,compareSNVToSkin)
#bloodComparison_old=apply(cellGenotypeMatrix,2,compareSNVToBloodFromMerged)
bloodComparison=apply(cellGenotypeMatrix,2,compareSNVToBloodGenotypesFromMerged) #lapply is element-wise...

# Compate SNV per cell
skinButNotBlood=setdiff(snp.sk,snp.bd)
skinNblood=snp.sk[! snp.sk %in% skinButNotBlood]
bloodWskinSNV=data.frame(bloodSNV=json.bd$model1_g0_call,row.names = snp.bd)[skinNblood,]

#validCells <- colSums((ref.sk+alt.sk)>=countLimit)>25
validCells <- colSums((ref.sk+alt.sk)>=countLimit)>0
hist(skinComparison[1,validCells & json.sk$model2_call==0],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno1-assigned cells compared to Geno1(Skin)")
hist(skinComparison[1,validCells] ,xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="All cells compared to Geno1")
hist(skinComparison[2,validCells & json.sk$model2_call==1],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno2-assigned cells compared to Geno2(Skin)")
hist(skinComparison[2,validCells] ,xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="All cells compared to Geno2")
hist(bloodComparison[1,validCells & json.sk$model2_call==0],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno1-assigned cells (skin) compared to Geno1 merged")
hist(bloodComparison[2,validCells & json.sk$model2_call==0],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno1-assigned cells (skin) compared to Geno2 merged")
hist(bloodComparison[validCells & json.sk$model2_call==1],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno2-assigned cells (skin) compared to merged")
hist(bloodComparison[validCells ],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="all cells (skin) compared to merged")
hist(bloodComparison,xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="all cells (skin, no minimal SNV-count per cell) compared to merged")
hist(bloodComparison[json.sk$model2_call==1],xlim=c(0,1),ylim=c(0,450),breaks = seq(0,1,0.05),
     main="Geno2-assigned cells (skin,no minimal SNV-count per cell) compared to merged")
hist(skinComparison[1,validCells] ,xlim=c(0.1,0.7),ylim=c(0,60),breaks = seq(0,1,0.005),
     main="All cells compared to Geno2")
hist(skinComparison[2,validCells] ,xlim=c(0.1,0.7),ylim=c(0,60),breaks = seq(0,1,0.005),
     main="All cells compared to Geno2")
hist(bloodComparison[validCells ] ,xlim=c(0.1,0.7),ylim=c(0,60),breaks = seq(0,1,0.005),
     main="all cells (skin) compared to merged")

ifSkinGeno1=json.sk$model2_call==0
ifSkinGeno2=json.sk$model2_call==1
identSkinGeno1=skinComparison[1,]
identSkinGeno2=skinComparison[2,]
identBloodGeno1=bloodComparison[1,]
identBloodGeno2=bloodComparison[2,]

plot(identSkinGeno1[validCells],identSkinGeno2[ validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "skin donor similarity", 
     main = "donor vs. host identity:\nall cells")
plot(identSkinGeno1[ifSkinGeno1 & validCells],identBloodGeno1[ifSkinGeno1 & validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\ncells genotyped as host")
plot(identSkinGeno1[ifSkinGeno1 & validCells],identBloodGeno2[ifSkinGeno1 & validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\ncells genotyped as host")
plot(identSkinGeno1[validCells],identBloodGeno1[validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\nall cells")
plot(identSkinGeno1[validCells],identBloodGeno2[validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\nall cells")

plot(identSkinGeno2[ifSkinGeno2 & validCells],identBloodGeno1[ifSkinGeno2 & validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin donor similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\ncells genotyped as donor")
plot(identSkinGeno2[ifSkinGeno2 & validCells],identBloodGeno2[ifSkinGeno2 & validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin donor similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\ncells genotyped as donor")
plot(identSkinGeno2[validCells],identBloodGeno1[validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\nall cells")
plot(identSkinGeno2[validCells],identBloodGeno2[validCells],
     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
     main = "donor vs. blood identity:\nall cells")

# Merge data to draw plots
dot.sim <- df <- list()
#df$Skin_genotype_1_vs_skin_genotype_2 <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                                    Similarity_to_skin_genotype_2=identSkinGeno2, Genotype=gt.sk)
#df$Blood_genotype_1_vs_blood_genotype_2 <- data.frame(Similarity_to_blood_genotype_1=identBloodGeno1, 
#                                                      Similarity_to_blood_genotype_2=identBloodGeno2, Genotype=gt.sk)
#df$Skin_genotype_1_vs_blood_genotype_1 <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                                     Similarity_to_blood_genotype_1=identBloodGeno1, Genotype=gt.sk)
#df$Skin_genotype_1_vs_blood_genotype_2 <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                                     Similarity_to_blood_genotype_2=identBloodGeno2, Genotype=gt.sk)
#df$Skin_genotype_2_vs_blood_genotype_1 <- data.frame(Similarity_to_skin_genotype_2=identSkinGeno2, 
#                                                     Similarity_to_blood_genotype_1=identBloodGeno1, Genotype=gt.sk)
#df$Skin_genotype_2_vs_blood_genotype_2 <- data.frame(Similarity_to_skin_genotype_2=identSkinGeno2, 
#                                                     Similarity_to_blood_genotype_2=identBloodGeno2, Genotype=gt.sk)
#df$Skin_genotype_1_vs_skin_genotype_2 <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                                    Similarity_to_skin_genotype_2=identSkinGeno2, 
#                                                    Genotype=factor(json.sk$model2_call+1))#paste0('genotype ',json.sk$model2_call+1))#gt.sk)
#df$Blood_host_vs_blood_donor <- data.frame(Similarity_to_blood_host=identBloodGeno1, 
#                                           Similarity_to_blood_donor=identBloodGeno2, Genotype=factor(json.sk$model2_call+1))
#df$Skin_genotype_1_vs_blood_host <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                               Similarity_to_blood_host=identBloodGeno1, Genotype=factor(json.sk$model2_call+1))
#df$Skin_genotype_1_vs_blood_donor <- data.frame(Similarity_to_skin_genotype_1=identSkinGeno1, 
#                                                Similarity_to_blood_donor=identBloodGeno2, Genotype=factor(json.sk$model2_call+1))
#df$Skin_genotype_2_vs_blood_host <- data.frame(Similarity_to_skin_genotype_2=identSkinGeno2, 
#                                               Similarity_to_blood_host=identBloodGeno1, Genotype=factor(json.sk$model2_call+1))
#df$Skin_genotype_2_vs_blood_donor <- data.frame(Similarity_to_skin_genotype_2=identSkinGeno2, 
#                                                Similarity_to_blood_donor=identBloodGeno2, Genotype=factor(json.sk$model2_call+1))
df$skin_host_vs_skin_donor <- data.frame(Similarity_to_skin_host=identSkinGeno1, 
                                                    Similarity_to_skin_donor=identSkinGeno2, Genotype=gt.sk)
df$blood_host_vs_blood_donor <- data.frame(Similarity_to_blood_host=identBloodGeno1, 
                                           Similarity_to_blood_donor=identBloodGeno2, Genotype=gt.sk)
df$skin_host_vs_blood_host <- data.frame(Similarity_to_skin_host=identSkinGeno1, 
                                               Similarity_to_blood_host=identBloodGeno1, Genotype=gt.sk)
df$skin_host_vs_blood_donor <- data.frame(Similarity_to_skin_host=identSkinGeno1, 
                                                Similarity_to_blood_donor=identBloodGeno2, Genotype=gt.sk)
df$skin_donor_vs_blood_host <- data.frame(Similarity_to_skin_donor=identSkinGeno2, 
                                               Similarity_to_blood_host=identBloodGeno1, Genotype=gt.sk)
df$skin_donor_vs_blood_donor <- data.frame(Similarity_to_skin_donor=identSkinGeno2, 
                                                Similarity_to_blood_donor=identBloodGeno2, Genotype=gt.sk)

# Draw plots
#t <- ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.8, alpha = 0.5) +
#  scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) +
#  scale_x_continuous(breaks = c(-20,0,20)) + facet_wrap(~variable, ncol = 4) + theme_custom + xlab('') + ylab('')
#ggplot(data.frame(x=identSkinGeno1, y=identSkinGeno2)) + geom_point(aes(x=x,y=y)) + theme_custom
#ggplot(data.frame(x=identSkinGeno1[ifSkinGeno1], y=identSkinGeno2[ifSkinGeno1])) + geom_point(aes(x=x,y=y)) + theme_custom
#ggplot(df$Skin_genotype_1_vs_skin_genotype_2, 
#                         aes_string(x=colnames(df$Skin_genotype_1_vs_skin_genotype_2)[1],y=colnames(df$Skin_genotype_1_vs_skin_genotype_2)[2], color='Genotype')) + 
#  geom_point(shape=16, size=1, alpha=0.7) + theme_custom + 
#  scale_x_continuous(gsub('_',' ',colnames(df$Skin_genotype_1_vs_skin_genotype_2)[1]), limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
#  scale_y_continuous(gsub('_',' ',colnames(df$Skin_genotype_1_vs_skin_genotype_2)[2]), limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
#  scale_color_manual('Genotype', values = c('black', 'royalblue4'), 
#                     guide = guide_legend(override.aes = list(shape = 16,size = 4, alpha = 0.7))) + 
#  ggtitle(gsub('_',' ','Skin_genotype_1_vs_skin_genotype_2'))
dot.sim <- Map(x=names(df), function(x) 
  ggplot(df[[x]], aes_string(x=colnames(df[[x]])[1],y=colnames(df[[x]])[2], color='Genotype', fill='Genotype')) + 
    geom_point(shape=21, size=1.5, alpha=0.7) + theme_custom + 
    scale_x_continuous(gsub('_',' ',colnames(df[[x]])[1]), limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
    scale_y_continuous(gsub('_',' ',colnames(df[[x]])[2]), limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
    scale_fill_manual('Skin', values = c('white', 'black')) +#c('black', 'white')) +
    scale_color_manual('Skin', values = c('grey40', 'grey40'), 
                       guide = guide_legend(override.aes = list(shape = 21,size = 4, alpha = 0.7))) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + ggtitle('HSCT skin', subtitle = paste(gsub('_',' ',x), 'SNV similarity')) )

#-----------------------------------------------------------------------------------------------------------------#
# Skin genotypes similarity to blood for each single cell
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
#                                       to = c('CD4 T','CD4 T','CD4 T','CD8 T','NK cell','CD8 T','NKT',
#                                              'B cell','gdT','CD4 T','B cell','Treg','n.d.'))
#pd$hd.bd$celltype <- factor(pd$hd.bd$celltype, levels = levels(pd$hd.bd$celltype)[c(1,2,6,7,4,3,5,8)])
pd$tp.bd$subset <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), 
                                   to = c('CD4 Tn','CD8 Tn','B cell','CD4 Tm','CD8 Tm','NK cell','gdT',
                                          'CD4 Tn','gdT','CD8 Tn','B cell','B cell','Treg'))
pd$tp.bd$subset <- factor(pd$tp.bd$subset, levels = levels(pd$tp.bd$subset)[c(1,4,2,5,8,7,6,3)])
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
#bc <- list(tp.sk = gsub('-1', '', read.delim('../202003/figures/GA_AN0355_sc_skin_blood_genotyping/barcodes/GA_AN0355_sc_genotyping_barcodes_skin_host.tsv',as.is = T,header = F)[,1]), 
#           tp.bd = gsub('-1', '', read.delim('../202003/figures/GA_AN0355_sc_skin_blood_genotyping/barcodes/GA_AN0355_sc_genotyping_barcodes_merged_blood_only_host.tsv',as.is = T,header = F)[,1]))
bc <- list(tp.sk = gsub('-1', '', bc.gt.sk[bc.gt.sk$genotype == "host",]$barcode), 
           tp.bd = gsub('-1', '', bc.gt.mg[bc.gt.mg$barcode %in% bc.bd & bc.gt.mg$genotype == 'host',]$barcode))

# Keep only the barcodes present in seurat object after fitlering steps
bc$tp.sk <- bc$tp.sk[bc$tp.sk %in% names(Idents(pd$tp.sk))]
bc$tp.bd <- bc$tp.bd[bc$tp.bd %in% names(Idents(pd$tp.bd))]

# Classify idents according to barcodes from host and donor genotypes
pd$tp.sk[['genotype']] <- names(Idents(pd$tp.sk))
pd$tp.sk[['genotype']] <- ifelse(pd$tp.sk$genotype %in% bc$tp.sk, 'host', 'donor')
pd$tp.bd[['genotype']] <- names(Idents(pd$tp.bd))
pd$tp.bd[['genotype']] <- ifelse(pd$tp.bd$genotype %in% bc$tp.bd, 'host', 'donor')

# Number of filtered cells
length(gt.sk) #2051
length(pd$tp.sk[['genotype']]$genotype) #2048
length(gt.bd) #7508
length(pd$tp.bd[['genotype']]$genotype) #7473

# Rearrange cluster numbers to match subsets in all four samples
for(x in names(pd)) pd[[x]]$cluster_old <- pd[[x]]$cluster
#pd$hd.bd$cluster_new <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), to = c(2,3,4,5,9,6,8,12,10,1,11,7,13))
#pd$hd.bd$cluster <- plyr::mapvalues(pd$hd.bd$cluster, from = levels(factor(pd$hd.bd$cluster)), to = c(2,3,4,5,10,6,9,12,8,1,11,7,13))
#pd$hd.bd$cluster <- factor(pd$hd.bd$cluster, levels = sort(as.numeric(levels(factor(pd$hd.bd$cluster)))))
#pd$hd.sk$cluster <- plyr::mapvalues(pd$hd.sk$cluster, from = levels(factor(pd$hd.sk$cluster)), to = c(2,1,3,4,7,5,12,11,6,10,13,9,8))
#pd$hd.sk$cluster <- factor(pd$hd.sk$cluster, levels = sort(as.numeric(levels(factor(pd$hd.sk$cluster)))))
pd$tp.bd$cluster <- plyr::mapvalues(pd$tp.bd$cluster, from = levels(factor(pd$tp.bd$cluster)), to = c(2,4,11,3,6,10,8,1,9,5,12,13,7))
pd$tp.bd$cluster <- factor(pd$tp.bd$cluster, levels = sort(as.numeric(levels(factor(pd$tp.bd$cluster)))))
pd$tp.sk$cluster <- plyr::mapvalues(pd$tp.sk$cluster, from = levels(factor(pd$tp.sk$cluster)), to = c(9,1,4,2,10,3,5,6,7,12,8,11,13))
pd$tp.sk$cluster <- factor(pd$tp.sk$cluster, levels = sort(as.numeric(levels(factor(pd$tp.sk$cluster)))))
#for(x in names(pd)) pd[[x]]$subset_cluster <- paste0(pd[[x]]@meta.data$cluster_new,' - ',pd[[x]]@meta.data$celltype)
#pd$hd.bd$subset_cluster <- factor(pd$hd.bd$subset_cluster, levels = levels(factor(pd$hd.bd$subset_cluster))[c(1,6:13,2:5)])
#pd$hd.sk$subset_cluster <- factor(pd$hd.sk$subset_cluster, levels = levels(factor(pd$hd.sk$subset_cluster))[c(1,6:13,2:5)])
#pd$tp.bd$subset_cluster <- factor(pd$tp.bd$subset_cluster, levels = levels(factor(pd$tp.bd$subset_cluster))[c(1,6:13,2:5)])
#pd$tp.sk$subset_cluster <- factor(pd$tp.sk$subset_cluster, levels = levels(factor(pd$tp.sk$subset_cluster))[c(1,6:13,2:5)])
#for(x in names(pd)) pd[[x]]$subset_cluster <- paste0(ifelse(as.numeric(pd[[x]]@meta.data$cluster)<10, paste0('  ',pd[[x]]@meta.data$cluster), paste0(pd[[x]]@meta.data$cluster)),' - ',pd[[x]]@meta.data$celltype)

# Annotate subsets
pd$tp.bd$subset_host <- as.character(pd$tp.bd$subset)
pd$tp.bd$subset_host[names(pd$tp.bd$subset) %in% bc$tp.bd] <- c('host')
pd$tp.bd$subset_host <- factor(pd$tp.bd$subset_host, levels = levels(factor(pd$tp.bd$subset_host))[c(3,2,5,4,6,9,8,1,7)])
pd$tp.sk$subset_host <- as.character(pd$tp.sk$subset)
pd$tp.sk$subset_host[names(pd$tp.sk$subset) %in% bc$tp.sk] <- c('host')
pd$tp.sk$subset_host <- factor(pd$tp.sk$subset_host, levels = levels(factor(pd$tp.sk$subset_host))[c(2,3,4,8,7,6,1,5)])

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of genotypes between subsets
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentages
absHost <- t(cbind(skin=table(pd$tp.sk$genotype), blood=table(pd$tp.bd$genotype)))
#absHost <- rbind(absHost, total=colSums(absHost))
#pctHost <- round(100*(pctHost / rowSums(pctHost)), digits = 1)
#pctHost <- as.matrix(100*(absHost[1:2,] / rowSums(absHost[1:2,])))
pctHost <- as.matrix(100*(absHost / rowSums(absHost)))
pctHost <- reshape2::melt(pctHost)

# Plot percentages
barplot.pct <- ggplot(pctHost, aes(x = Var1, y = value)) + 
  geom_bar(aes(fill = Var2), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = plyr::ddply(reshape2::melt(as.matrix(absHost)), "Var1", plyr::summarise, total = sum(value)), 
            aes(y = 95, label = total), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = pctHost[pctHost$Var2=='host',], 
            aes(y = value+5, label =  paste0(round(value, digits = 1),' %') ), size = 3, angle = 0, hjust = 0.5) + 
  scale_x_discrete('') + 
  scale_y_continuous("% of total cells", breaks = seq(0,100,20),
                     #ceiling(pctHost[pctHost$Var1=='skin' & pctHost$Var2=='host',]$value), 
                     #ceiling(pctHost[pctHost$Var1=='blood' & pctHost$Var2=='host',]$value)), 
                     expand = c(0,0)) + #labels = scales::percent_format(), 
  scale_fill_manual("Genotype", values = c('white', 'black')) +
  ggtitle("HSCT patient", subtitle = "Genotype distribution") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 0, hjust = 0.5), aspect.ratio = NULL)
barplot.pct

## Cummulative percentage for hsct skin subsets split between host and donor
#freq.tdc <- data.frame('subset' = rownames(table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)), 
#                       'skin_donor' = table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)[,1])
#freq.thc <- data.frame('subset' = rownames(table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)), 
#                       'skin_host' = table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)[,2])
#freqc <- merge(freq.tdc, freq.thc, by = "subset", all = T)
#freqc <- reshape2::melt(freqc, id.vars = "subset", value.name = "Freq")
#freqc$variable <- sub("\\.", ")", sub("\\.", "(", gsub("_", " ", freqc$variable)))
#freqc$variable <- factor(plyr::mapvalues(freqc$variable, from = levels(factor(freqc$variable)), to = c('donor', 'host')), 
#                         levels = c('donor', 'host'))
#
## Calculate percentage from total cells
#cn.totalc <- plyr::ddply(freqc, "subset", plyr::summarise, sum = sum(Freq))
#for(i in 1:nrow(freqc)) freqc[i, "percent"] <- 100*freqc[i, "Freq"] / cn.totalc[cn.totalc$subset == freqc$subset[i], "sum"]
#freqc$subset <- factor(freqc$subset, levels = levels(factor(freqc$subset))[c(2,3,4,7,6,5,1)])
#
## Plot distribution of genotypes per subset
#freqc.dist <- ggplot(freqc, aes(x = subset, y = percent)) + 
#  geom_bar(aes(fill = variable), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
#  geom_text(data = cn.totalc[cn.totalc$subset %in% levels(cn.totalc$subset)[-6],], aes(y = 95, label = sum), size = 3, angle = 90, hjust = 1) + 
#  geom_text(data = freqc[freqc$variable %in% levels(freqc$variable)[c(2)] & freqc$subset %in% levels(freqc$subset)[6:7],], 
#            aes(y = percent+5, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5) + 
#  geom_text(data = freqc[freqc$variable %in% levels(freqc$variable)[c(2)] & freqc$subset %in% levels(freqc$subset)[1:4],], 
#            aes(y = percent-5, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
#  geom_text(data = cn.totalc[cn.totalc$subset %in% levels(cn.totalc$subset)[6],], aes(y = 95, label = sum), 
#            size = 3, angle = 90, hjust = 1, color = 'white') + 
#  geom_text(data = freqc[freqc$variable %in% levels(freqc$variable)[c(2)] & freqc$subset %in% levels(freqc$subset)[5],], 
#            aes(y = percent-15, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
#  scale_x_discrete("", limits = levels(freqc$subset)) + 
#  scale_y_continuous("Percentage of total cells", breaks = seq(0,100,20), expand = c(0,0)) + 
#  scale_fill_manual("Genotype", values = c('white', 'black')) +
#  ggtitle("scRNA-Seq", subtitle = "Genotype distribution per subset") + theme_custom + #theme(aspect.ratio = NULL)
#  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
#freqc.dist

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of genotypes between subsets within T cells
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentages
absHostT <- t(cbind(skin=table(pd$tp.sk$genotype[pd$tp.sk$subset %in% levels(pd$tp.sk$subset)[-c(6,7)]]), 
                    blood=table(pd$tp.bd$genotype[pd$tp.bd$subset %in% levels(pd$tp.bd$subset)[-c(7,8)]]) ))
#absHost <- rbind(absHost, total=colSums(absHost))
#pctHost <- round(100*(pctHost / rowSums(pctHost)), digits = 1)
#pctHost <- as.matrix(100*(absHost[1:2,] / rowSums(absHost[1:2,])))
pctHostT <- as.matrix(100*(absHostT / rowSums(absHostT)))
pctHostT <- reshape2::melt(pctHostT)

# Plot percentages
barplot.pctT <- ggplot(pctHostT, aes(x = Var1, y = value)) + 
  geom_bar(aes(fill = Var2), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = plyr::ddply(reshape2::melt(as.matrix(absHostT)), "Var1", plyr::summarise, total = sum(value)), 
            aes(y = 95, label = total), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = pctHostT[pctHostT$Var2=='host',], 
            aes(y = value+5, label =  paste0(round(value, digits = 1),' %') ), size = 3, angle = 0, hjust = 0.5) + 
  scale_x_discrete('') + 
  scale_y_continuous("% of total cells", breaks = seq(0,100,20),
                     #ceiling(pctHost[pctHost$Var1=='skin' & pctHost$Var2=='host',]$value), 
                     #ceiling(pctHost[pctHost$Var1=='blood' & pctHost$Var2=='host',]$value)), 
                     expand = c(0,0)) + #labels = scales::percent_format(), 
  scale_fill_manual("Genotype", values = c('white', 'black')) +
  ggtitle("HSCT patient", subtitle = "Genotype distribution within T cells") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 0, hjust = 0.5), aspect.ratio = NULL)
barplot.pctT

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of genotypes between subsets split in T cells and non-T cells
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentages
absHostTnT <- t(cbind(skin_T_cell=table(pd$tp.sk$genotype[pd$tp.sk$subset %in% levels(pd$tp.sk$subset)[-c(6,7)]]), 
                      blood_T_cell=table(pd$tp.bd$genotype[pd$tp.bd$subset %in% levels(pd$tp.bd$subset)[-c(7,8)]]), 
                      `skin_non-T_cell`=table(pd$tp.sk$genotype[pd$tp.sk$subset %in% levels(pd$tp.sk$subset)[c(6,7)]]), 
                      `blood_non-T_cell`=table(pd$tp.bd$genotype[pd$tp.bd$subset %in% levels(pd$tp.bd$subset)[c(7,8)]]) ))
#absHost <- rbind(absHost, total=colSums(absHost))
#pctHost <- round(100*(pctHost / rowSums(pctHost)), digits = 1)
#pctHost <- as.matrix(100*(absHost[1:2,] / rowSums(absHost[1:2,])))
pctHostTnT <- as.matrix(100*(absHostTnT / rowSums(absHostTnT)))
pctHostTnT <- reshape2::melt(pctHostTnT)
pctHostTnT$Var1 <- factor(gsub('_',' ',pctHostTnT$Var1), levels = levels(factor(gsub('_',' ',pctHostTnT$Var1)))[c(4,3,2,1)])

pctHostTnTn <- plyr::ddply(reshape2::melt(as.matrix(absHostTnT)), "Var1", plyr::summarise, total = sum(value))
pctHostTnTn$Var1 <- gsub('_',' ',pctHostTnTn$Var1)

# Plot percentages
barplot.pctTnT <- ggplot(pctHostTnT, aes(x = Var1, y = value)) + 
  geom_bar(aes(fill = Var2), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = pctHostTnTn, aes(y = 95, label = total), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = pctHostTnT[pctHostTnT$Var2=='host',], 
            aes(y = value+5, label =  paste0(round(value, digits = 1),' %') ), size = 3, angle = 0, hjust = 0.5) + 
  scale_x_discrete('') + 
  scale_y_continuous("% of total cells", breaks = seq(0,100,20),
                     #ceiling(pctHost[pctHost$Var1=='skin' & pctHost$Var2=='host',]$value), 
                     #ceiling(pctHost[pctHost$Var1=='blood' & pctHost$Var2=='host',]$value)), 
                     expand = c(0,0)) + #labels = scales::percent_format(), 
  scale_fill_manual("Genotype", values = c('white', 'black')) +
  ggtitle("HSCT patient", subtitle = "Genotype distribution") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
barplot.pctTnT

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of subsets between hsct skin genotypes
#-----------------------------------------------------------------------------------------------------------------#

# Cummulative percentage for hsct skin subsets split between host and donor
freq.tdc.sk <- data.frame('subset' = rownames(table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)), 
                       'skin_donor' = table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)[,1])
freq.thc.sk <- data.frame('subset' = rownames(table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)), 
                       'skin_host' = table(pd$tp.sk@meta.data$subset, pd$tp.sk@meta.data$genotype)[,2])
freqc.sk <- merge(freq.tdc.sk, freq.thc.sk, by = "subset", all = T)
freqc.sk <- reshape2::melt(freqc.sk, id.vars = "subset", value.name = "Freq")
freqc.sk$variable <- sub("\\.", ")", sub("\\.", "(", gsub("_", " ", freqc.sk$variable)))
freqc.sk$variable <- factor(plyr::mapvalues(freqc.sk$variable, from = levels(factor(freqc.sk$variable)), to = c('donor', 'host')), 
                         levels = c('donor', 'host'))

# Calculate percentage from total cells
cn.totalc.sk <- plyr::ddply(freqc.sk, "subset", plyr::summarise, sum = sum(Freq))
for(i in 1:nrow(freqc.sk)) freqc.sk[i, "percent"] <- 100*freqc.sk[i, "Freq"] / cn.totalc.sk[cn.totalc.sk$subset == freqc.sk$subset[i], "sum"]
freqc.sk$subset <- factor(freqc.sk$subset, levels = levels(factor(freqc.sk$subset))[c(2,3,4,7,6,5,1)])

# Plot distribution of genotypes per subset
freqc.dist.sk <- ggplot(freqc.sk, aes(x = subset, y = percent)) + 
  geom_bar(aes(fill = variable), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = cn.totalc.sk[cn.totalc.sk$subset %in% levels(cn.totalc.sk$subset)[-6],], aes(y = 95, label = sum), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = freqc.sk[freqc.sk$variable %in% levels(freqc.sk$variable)[c(2)] & freqc.sk$subset %in% levels(freqc.sk$subset)[6:7],], 
            aes(y = percent+5, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5) + 
  geom_text(data = freqc.sk[freqc.sk$variable %in% levels(freqc.sk$variable)[c(2)] & freqc.sk$subset %in% levels(freqc.sk$subset)[1:4],], 
            aes(y = percent-5, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
  geom_text(data = cn.totalc.sk[cn.totalc.sk$subset %in% levels(cn.totalc.sk$subset)[6],], aes(y = 95, label = sum), 
            size = 3, angle = 90, hjust = 1, color = 'white') + 
  geom_text(data = freqc.sk[freqc.sk$variable %in% levels(freqc.sk$variable)[c(2)] & freqc.sk$subset %in% levels(freqc.sk$subset)[5],], 
            aes(y = percent-15, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
  scale_x_discrete("", limits = levels(freqc.sk$subset)) + 
  scale_y_continuous("Percentage of total cells", breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_manual("Genotype", values = c('white', 'black')) +
  ggtitle(unique(pd$tp.sk$sample), subtitle = "Genotype distribution per subset") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
freqc.dist.sk

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of subsets between hsct blood genotypes
#-----------------------------------------------------------------------------------------------------------------#

# Cummulative percentage for hsct skin subsets split between host and donor
freq.tdc.bd <- data.frame('subset' = rownames(table(pd$tp.bd@meta.data$subset, pd$tp.bd@meta.data$genotype)), 
                       'blood_donor' = table(pd$tp.bd@meta.data$subset, pd$tp.bd@meta.data$genotype)[,1])
freq.thc.bd <- data.frame('subset' = rownames(table(pd$tp.bd@meta.data$subset, pd$tp.bd@meta.data$genotype)), 
                       'blood_host' = table(pd$tp.bd@meta.data$subset, pd$tp.bd@meta.data$genotype)[,2])
freqc.bd <- merge(freq.tdc.bd, freq.thc.bd, by = "subset", all = T)
freqc.bd <- reshape2::melt(freqc.bd, id.vars = "subset", value.name = "Freq")
freqc.bd$variable <- sub("\\.", ")", sub("\\.", "(", gsub("_", " ", freqc.bd$variable)))
freqc.bd$variable <- factor(plyr::mapvalues(freqc.bd$variable, from = levels(factor(freqc.bd$variable)), to = c('donor', 'host')), 
                         levels = c('donor', 'host'))

# Calculate percentage from total cells
cn.totalc.bd <- plyr::ddply(freqc.bd, "subset", plyr::summarise, sum = sum(Freq))
for(i in 1:nrow(freqc.bd)) freqc.bd[i, "percent"] <- 100*freqc.bd[i, "Freq"] / cn.totalc.bd[cn.totalc.bd$subset == freqc.bd$subset[i], "sum"]
freqc.bd$subset <- factor(freqc.bd$subset, levels = levels(factor(freqc.bd$subset))[c(3,2,5,4,6,8,7,1)])

# Plot distribution of genotypes per subset
freqc.dist.bd <- ggplot(freqc.bd, aes(x = subset, y = percent)) + 
  geom_bar(aes(fill = variable), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = cn.totalc.bd[cn.totalc.bd$subset %in% levels(cn.totalc.bd$subset),], aes(y = 95, label = sum), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = freqc.bd[freqc.bd$variable %in% levels(freqc.bd$variable)[c(2)] & freqc.bd$subset %in% levels(freqc.bd$subset),], 
            aes(y = percent+5, label =  paste0(round(percent, digits = 2),'%') ), size = 3, angle = 0, hjust = 0.5) + 
  #geom_text(data = cn.totalc.bd[cn.totalc.bd$subset %in% levels(cn.totalc.bd$subset)[6],], aes(y = 95, label = sum), 
  #          size = 3, angle = 90, hjust = 1, color = 'white') + 
  #geom_text(data = freqc.bd[freqc.bd$variable %in% levels(freqc.bd$variable)[c(2)] & freqc.bd$subset %in% levels(freqc.bd$subset)[5],], 
  #          aes(y = percent-15, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
  scale_x_discrete("", limits = levels(freqc.bd$subset)) + 
  scale_y_continuous("Percentage of total cells", breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_manual("Genotype", values = c('white', 'black')) +
  ggtitle(unique(pd$tp.bd$sample), subtitle = "Genotype distribution per subset") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
freqc.dist.bd

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of genotypes between CD4 and CD8 subsets
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentages
#freqs <- rbind(merge(transform(freq.tdc.sk[1:2,], subset = paste('skin',subset)), transform(freq.thc.sk[1:2,], subset = paste('skin',subset)), by = 'subset'), 
#               merge(transform(freq.tdc.bd[4:5,], subset = paste('blood',sub('m','',subset))), transform(freq.thc.bd[4:5,], subset = paste('blood',sub('m','',subset))), by = 'subset'))
freqs <- merge(merge(freq.tdc.sk[1:2,], freq.thc.sk[1:2,], by = 'subset'), 
               transform(merge(freq.tdc.bd[c(2,4),], freq.thc.bd[c(2,4),], by = 'subset'), subset = sub('m','',subset)), by = 'subset')
freqs <- reshape2::melt(freqs, id.vars = "subset", value.name = "Freq")
#freqs$tissue <- stringr::str_split(freqs$variable, "_", simplify = T)[,1]
#freqs$genotype <- stringr::str_split(freqs$variable, "_", simplify = T)[,2]
freqs$variable <- sub('_',' ',freqs$variable)
freqs$subset <- factor(freqs$subset, levels=levels(factor(as.character(freqs$subset)))[2:1])

# Calculate percentage from total cells
#cn.totalcs <- plyr::ddply(freqs, c('subset','tissue'), plyr::summarise, sum = sum(Freq))
#cn.totalcs <- plyr::ddply(freqs, c('subset','genotype'), plyr::summarise, sum = sum(Freq))
cn.totalcs <- plyr::ddply(freqs, c('variable'), plyr::summarise, sum = sum(Freq))
#for(i in 1:nrow(freqs)) freqs[i, "percent"] <- 100*freqs[i, "Freq"] / cn.totalcs[cn.totalcs$subset == freqs$subset[i] & cn.totalcs$tissue == freqs$tissue[i], "sum"]
#for(i in 1:nrow(freqs)) freqs[i, "percent"] <- 100*freqs[i, "Freq"] / cn.totalcs[cn.totalcs$subset == freqs$subset[i] & cn.totalcs$genotype == freqs$genotype[i], "sum"]
for(i in 1:nrow(freqs)) freqs[i, "percent"] <- 100*freqs[i, "Freq"] / cn.totalcs[cn.totalcs$variable == freqs$variable[i], "sum"]
freqs$subset <- factor(freqs$subset, levels = levels(factor(freqs$subset)))

# Plot percentages
barplot.pcts <-  ggplot(freqs, aes(x = variable, y = percent)) + 
  geom_bar(aes(fill = subset), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = cn.totalcs[cn.totalcs$variable %in% levels(factor(cn.totalcs$variable)),], aes(y = 95, label = sum), size = 3, angle = 90, hjust = 1) + 
  geom_text(data = freqs[freqs$subset %in% levels(freqs$subset)[2],], color = 'white', 
            aes(y = percent-5, label =  paste0(round(percent, digits = 2),'%') ), size = 3, angle = 0, hjust = 0.5) + 
  #geom_text(data = cn.totalc.bd[cn.totalc.bd$subset %in% levels(cn.totalc.bd$subset)[6],], aes(y = 95, label = sum), 
  #          size = 3, angle = 90, hjust = 1, color = 'white') + 
  #geom_text(data = freqc.bd[freqc.bd$variable %in% levels(freqc.bd$variable)[c(2)] & freqc.bd$subset %in% levels(freqc.bd$subset)[5],], 
  #          aes(y = percent-15, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
  scale_x_discrete("", limits = levels(freqs$variable)) + 
  scale_y_continuous("Percentage of total cells", breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_manual("Subset", values = c('white', 'black'), breaks = levels(freqs$subset)[2:1]) +
  ggtitle('HSCT patient', subtitle = "CD4 / CD8 distribution per genotype") + theme_custom + #theme(aspect.ratio = NULL)
  theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
barplot.pcts

#-----------------------------------------------------------------------------------------------------------------#
# Frequencies of naive and memory between host CD4 and CD8 subsets
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentages
freqsh <- freq.thc.bd[c(1:4),]
freqsh <- reshape2::melt(freqsh, id.vars = "subset", value.name = "Freq")
freqsh$variable <- sub('_',' ',freqsh$variable)
freqsh$subset <- factor(freqsh$subset, levels=levels(factor(as.character(freqsh$subset)))[c(3,4,1,2)])
freqsh$celltype <- stringr::str_split(freqsh$subset, " ", simplify = T)[,1]
freqsh$celltype <- factor(freqsh$celltype, levels=levels(factor(freqsh$celltype))[1:2])
freqsh$state <- stringr::str_split(freqsh$subset, " ", simplify = T)[,2]
freqsh$state <- factor(freqsh$state, levels=levels(factor(freqsh$state))[2:1])

# Calculate percentage from total cells
cn.totalcsh <- plyr::ddply(freqsh, c('celltype'), plyr::summarise, sum = sum(Freq))
for(i in 1:nrow(freqsh)) freqsh[i, "percent"] <- 100*freqsh[i, "Freq"] / cn.totalcsh[cn.totalcsh$celltype == freqsh$celltype[i], "sum"]
freqsh$subset <- factor(freqsh$subset, levels = levels(factor(freqsh$subset)))

# Plot percentages
barplot.pctsh <-  ggplot(freqsh, aes(x = celltype, y = percent)) + 
  geom_bar(aes(fill = state), stat = "identity", width = 0.7, alpha = 1, color = 'black') + 
  geom_text(data = cn.totalcsh[cn.totalcsh$celltype %in% levels(factor(cn.totalcsh$celltype)),], aes(y = 10, label = sum), size = 3, angle = 90, hjust = 1, color='white') + 
  geom_text(data = freqsh[freqsh$subset %in% levels(freqsh$subset)[c(1,3)],], color = 'white', 
            aes(y = percent-5, label =  paste0(round(percent, digits = 2),'%') ), size = 3, angle = 0, hjust = 0.5) + 
  #geom_text(data = cn.totalc.bd[cn.totalc.bd$subset %in% levels(cn.totalc.bd$subset)[6],], aes(y = 95, label = sum), 
  #          size = 3, angle = 90, hjust = 1, color = 'white') + 
  #geom_text(data = freqc.bd[freqc.bd$variable %in% levels(freqc.bd$variable)[c(2)] & freqc.bd$subset %in% levels(freqc.bd$subset)[5],], 
  #          aes(y = percent-15, label =  paste0(round(percent, digits = 1),'%') ), size = 3, angle = 0, hjust = 0.5, color = 'white') + 
  scale_x_discrete("", limits = levels(freqsh$variable)) + 
  scale_y_continuous("Percentage of total cells", breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_manual("State", values = c('white', 'black'), breaks = levels(freqsh$state)[1:2]) +
  ggtitle('HSCT patient', subtitle = "Cell state distribution by subset") + theme_custom + theme(aspect.ratio = NULL)
  #theme(axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1), aspect.ratio = NULL)
barplot.pctsh

#-----------------------------------------------------------------------------------------------------------------#
# Vizualization in tSNE embeddings
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data to visualize tSNE embeddings
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
             subset_host = pd[[x]]@meta.data$subset_host, 
             #subset_cluster = paste0(pd[[x]]@meta.data$cluster_new,' - ',pd[[x]]@meta.data$celltype),
             #subset_cluster = pd[[x]]@meta.data$subset_cluster,
             genotype = pd[[x]]@meta.data$genotype,
             male = ifelse(as.matrix(GetAssayData(pd[[x]])["RPS4Y1",])>0, 'male', 'n.d.'),
             female = ifelse(as.matrix(GetAssayData(pd[[x]])["XIST",])>0, 'female', 'n.d.'),
             row.names = 1:length(colnames(pd[[x]]))) })

# Plot tSNE for genotypes
ts.gt <- Map(x=names(pd), function(x) {
  ggplot(df.ts[[x]] %>% arrange(genotype), aes(x = tx, y = ty, color = genotype, fill = genotype)) + 
    geom_point(shape = 21, size = 0.8, alpha = 1, stroke = 0.5) + 
    scale_color_manual("Genotype", values = c("grey70", "grey40"), guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
    scale_fill_manual("Genotype", values = c("white", "black")) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + #facet_wrap(~variable, ncol = 4) + 
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample, subtitle = "SNV genotyping") })

# Plot tSNE for genotypes over subsets
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
  n <- length(levels(pd[[x]]@meta.data$subset_host))-1
  ggplot(df.ts[[x]] %>% arrange(subset_host), aes(x = tx, y = ty, color = subset_host, alpha = subset_host, size = subset_host)) + 
    geom_point(shape = 16, stroke = 0.5) + 
    scale_color_manual("Subset", values = c(coln$subset[[x]], 'black'), guide = guide_legend(override.aes = list(size = 4, alpha = c(rep(ifelse(x=='tp.sk',0.8,0.6),n),1)))) + 
    scale_alpha_manual("Subset", values = c(rep(ifelse(x=='tp.sk',0.8,0.3),n),1)) +
    scale_size_manual("Subset", values = c(rep(0.8,n),ifelse(x=='tp.sk',0.8,1.2))) + 
    scale_x_continuous(breaks = seq(-25,25,25)) + scale_y_continuous(breaks = seq(-25,25,25)) + #facet_wrap(~variable, ncol = 4) + 
    theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(pd[[x]]$sample, subtitle = "SNV genotyping") })

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Save barcodes
#write.csv(gb, file = paste0(dir.figs, "_barcodes.csv"), quote = F, row.names = T)
write.csv(bc.gt.sk, file = paste0(dir.figs.bc, "_barcodes_skin.csv"), quote = F, row.names = T)
write.csv(bc.gt.bd, file = paste0(dir.figs.bc, "_barcodes_blood.csv"), quote = F, row.names = T)
write.csv(bc.gt.mg, file = paste0(dir.figs.bc, "_barcodes_merged.csv"), quote = F, row.names = T)
write.csv(bc.gt.mg[bc.gt.mg$barcode %in% bc.bd,], file = paste0(dir.figs.bc, "_barcodes_merged_blood_only.csv"), quote = F, row.names = T)
write.table(bc.gt.sk[bc.gt.sk$genotype == "host", "barcode"], file = paste0(dir.figs.bc, "_barcodes_skin_host.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.sk[bc.gt.sk$genotype == "donor", "barcode"], file = paste0(dir.figs.bc, "_barcodes_skin_donor.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.bd[bc.gt.bd$genotype == "host", "barcode"], file = paste0(dir.figs.bc, "_barcodes_blood_host.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.bd[bc.gt.bd$genotype == "donor", "barcode"], file = paste0(dir.figs.bc, "_barcodes_blood_donor.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.mg[bc.gt.mg$genotype == "host", "barcode"], file = paste0(dir.figs.bc, "_barcodes_merged_host.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.mg[bc.gt.mg$genotype == "donor", "barcode"], file = paste0(dir.figs.bc, "_barcodes_merged_donor.tsv"), 
            quote = F, row.names = F, col.names = F)
write.table(bc.gt.mg[bc.gt.mg$barcode %in% bc.bd & bc.gt.mg$genotype == 'host',]$barcode, 
            file = paste0(dir.figs.bc, "_barcodes_merged_blood_only_host.tsv"), quote = F, row.names = F, col.names = F)
write.table(bc.gt.mg[bc.gt.mg$barcode %in% bc.bd & bc.gt.mg$genotype == 'donor',]$barcode, 
            file = paste0(dir.figs.bc, "_barcodes_merged_blood_only_donor.tsv"), quote = F, row.names = F, col.names = F)

# Save tables
#data.frame(tb.gs, `x`=rownames(tb.gs), row.names = 1:nrow(tb.gs))[c(length(tb.gs)+1,1:(length(tb.gs)))]
#tb <- function(x) htmlTable::htmlTable(data.frame(x, `x`=rownames(x), row.names = 1:nrow(x))[c(length(x)+1,1:(length(x)))], rnames = FALSE)
#tb <- function(x) data.frame(x, `x`=rownames(x), row.names = 1:nrow(x))[c(length(x)+1,1:(length(x)))]
tt <- ttheme_minimal(rowhead=list(fg_params=list(fontface="bold")))
pdf(paste0(dir.figs, "_stats.pdf"), height = 16, width = 6)
#grid.arrange(tableGrob(tb(tb.gs), theme=tt), tableGrob(tb(tb.gt), theme=tt), tableGrob(tb(tb.mp), theme=tt), tableGrob(tb(tb.nc), theme=tt), 
#             tableGrob(tb(tb.nc), theme=tt), tableGrob(tb(tb.ss), theme=tt), tableGrob(tb(tb.ss2), theme=tt), tableGrob(tb(tb.ss3), theme=tt), 
#             nrow = 8, ncol = 1)
grid.arrange(tableGrob(tb.gs, theme=tt), tableGrob(tb.gt, theme=tt), tableGrob(tb.mp, theme=tt), tableGrob(tb.nc, theme=tt), 
             tableGrob(tb.ss, theme=tt), tableGrob(tb.ss2, theme=tt), tableGrob(tb.ss3, theme=tt), 
             nrow = 7, ncol = 1)
dev.off()
#pdf(paste0(dir.figs, "_stats.pdf"), height = 6, width = 6)
#marrangeGrob(list(tableGrob(tb.mp, theme=tt), tableGrob(tb.nc, theme=tt), 
#             tableGrob(tb.gs, theme=tt), tableGrob(tb.ss, theme=tt)), nrow = 1, ncol = 1)
#dev.off()

# Save histograms
ggsave(plot = hist.cs.sk, file = paste0(dir.figs, "_hist_skin_cells_per_snv.pdf"), height = 4, width = 4.5)
ggsave(plot = hist.sc.sk, file = paste0(dir.figs, "_hist_skin_snv_per_cell.pdf"), height = 4, width = 4.5)
ggsave(plot = hist.cs.bd, file = paste0(dir.figs, "_hist_blood_cells_per_snv.pdf"), height = 4, width = 4.5)
ggsave(plot = hist.sc.bd, file = paste0(dir.figs, "_hist_blood_snv_per_cell.pdf"), height = 4, width = 4.5)
ggsave(plot = hist.cs.mg, file = paste0(dir.figs, "_hist_merged_cells_per_snv.pdf"), height = 4, width = 4.5)
ggsave(plot = hist.sc.mg, file = paste0(dir.figs, "_hist_merged_snv_per_cell.pdf"), height = 4, width = 4.5)
lapply(names(df), function(x) ggsave(plot = dot.sim[[x]], file = paste0(dir.figs, "_single_cell_similarity_",x,".pdf"), height = 4, width = 4.5))

# Save plots for frequencies of subsets between host and donor skin cells
ggsave(plot = barplot.pct, file = paste0(dir.figs, "_distribution.pdf"), height = 4, width = 3.3)
ggsave(plot = barplot.pctT, file = paste0(dir.figs, "_distribution_within_tcells.pdf"), height = 4, width = 3.3)
ggsave(plot = barplot.pctTnT, file = paste0(dir.figs, "_distribution_split_t_nont_cells.pdf"), height = 4.9, width = 4.5)
ggsave(plot = freqc.dist.sk, file = paste0(dir.figs, "_distribution_skin.pdf"), height = 4.35, width = 6)
ggsave(plot = freqc.dist.bd, file = paste0(dir.figs, "_distribution_merged_blood_only.pdf"), height = 4.45, width = 6.5)
ggsave(plot = barplot.pcts, file = paste0(dir.figs, "_distribution_CD4_vs_CD8.pdf"), height = 4.6, width = 4.5)
ggsave(plot = barplot.pctsh, file = paste0(dir.figs, "_distribution_cell_state_by_subset.pdf"), height = 4, width = 3.3)

#Save tSNE plots
lapply(names(df), function(x) ggsave(plot = ts.gt[[x]], file = paste0(dir.figs, "_tsne_",ifelse(x=='tp.sk','skin','blood'),"_genotype.pdf"), height = 4.1, width = 5))
lapply(names(df), function(x) ggsave(plot = ts.sg[[x]], file = paste0(dir.figs, "_tsne_",ifelse(x=='tp.sk','skin','blood'),"_genotype_subset.pdf"), height = 4.1, width = 5))

#pdf(paste0(dir.figs, "_plots_similarity_skin_to_blood.pdf"), height = 4, width = 4)
#plot(identSkinGeno2[ifSkinGeno1],identSkinGeno1[ifSkinGeno1],
#     xlim=c(0,0.8),ylim=c(0,0.8), xlab = "skin host similarity", ylab = "skin donor similarity", 
#     main = "HSCT skin\n(host vs. donor identity)", cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno2[ifSkinGeno2 & validCells],identSkinGeno1[ifSkinGeno2 & validCells], cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0.65, alpha = 0.5))
#abline(a=0,b=1, lty=2)
#legend("topleft", inset = 0.01, legend=c("skin donor", "skin host"), col=c(rgb(red = 0, green = 0, blue = 0, alpha = 0.5), rgb(red = 0, green = 0, blue = 0.65, alpha = 0.5)), cex=0.7, pch = 16)
#dev.off()
#
#pdf(paste0(dir.figs, "_plots_similarity.pdf"), height = 4, width = 4)
#plot(identSkinGeno1[ifSkinGeno1 & validCells],identBlood[ifSkinGeno1 & validCells],
#     xlim=c(0,0.8),ylim=c(0,0.8), xlab = "skin donor similarity", ylab = "blood similarity", 
#     main = "donor vs. blood identity\n(all cells)", cex = 1, pch = 16, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno1[ifSkinGeno2 & validCells],identBlood[ifSkinGeno2 & validCells], cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#legend("topleft", inset = 0.01, legend=c("skin donor", "skin host"), col=c("red", "black"), cex=0.7, pch = 16)
#
#plot(identSkinGeno2[ifSkinGeno2 & validCells],identBlood[ifSkinGeno2 & validCells],
#     xlim=c(0,0.8),ylim=c(0,0.8), xlab = "skin host similarity", ylab = "blood similarity", 
#     main = "host vs. blood identity\n(all cells)", cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno2[ifSkinGeno1 & validCells],identBlood[ifSkinGeno1 & validCells], cex = 1, pch = 16, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
#legend("topright", inset = 0.01, legend=c("skin donor", "skin host"), col=c("red", "black"), cex=0.7, pch = 16)
#
#plot(identSkinGeno1[ifSkinGeno1 & validCells]/identSkinGeno2[ifSkinGeno1 & validCells],identBlood[ifSkinGeno1 & validCells],
#     xlim=c(0,3),ylim=c(0,0.8), xlab = "donor similarity / host similarity (skin)", ylab = "blood similarity", 
#     main = "donor/host vs. blood identity\n(all cells)", cex = 1, pch = 16, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno1[ifSkinGeno2 & validCells]/identSkinGeno2[ifSkinGeno2 & validCells],identBlood[ifSkinGeno2 & validCells], cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#legend("topleft", inset = 0.01, legend=c("skin donor", "skin host"), col=c("red", "black"), cex=0.7, pch = 16)
#
#plot(identSkinGeno2[ifSkinGeno1 & validCells]/identSkinGeno1[ifSkinGeno1 & validCells],identBlood[ifSkinGeno1 & validCells],
#     xlim=c(0,3),ylim=c(0,0.8), xlab = "host similarity / donor similarity (skin)", ylab = "blood similarity", 
#     main = "host/donor vs. blood identity\n(all cells)", cex = 1, pch = 16, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno2[ifSkinGeno2 & validCells]/identSkinGeno1[ifSkinGeno2 & validCells],identBlood[ifSkinGeno2 & validCells], cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#legend("topright", inset = 0.01, legend=c("skin donor", "skin host"), col=c("red", "black"), cex=0.7, pch = 16)
#
##plot(identSkinGeno2[ifSkinGeno1 & validCells],identSkinGeno1[ifSkinGeno1 & validCells],
#plot(identSkinGeno2[ifSkinGeno1],identSkinGeno1[ifSkinGeno1],
#     xlim=c(0,0.8),ylim=c(0,0.8), xlab = "skin host similarity", ylab = "skin donor similarity", 
#     main = "HSCT skin\n(host vs. donor identity)", cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#points(identSkinGeno2[ifSkinGeno2 & validCells],identSkinGeno1[ifSkinGeno2 & validCells], cex = 1, pch = 16, col = rgb(red = 0, green = 0, blue = 0.65, alpha = 0.5))
#abline(a=0,b=1, lty=2)
#legend("topleft", inset = 0.01, legend=c("skin donor", "skin host"), col=c(rgb(red = 0, green = 0, blue = 0, alpha = 0.5), rgb(red = 0, green = 0, blue = 0.65, alpha = 0.5)), cex=0.7, pch = 16)
##xtest=seq(0,1,0.01)
##ytest=xtest*1
##points(xtest,ytest)
#
#plot(identSkinGeno1[ifSkinGeno1 & validCells],identBlood[ifSkinGeno1 & validCells],
#     xlim=c(0,1),ylim=c(0,1), xlab = "skin donor similarity", ylab = "blood similarity", 
#     main = "donor vs. blood identity\n(cells genotyped as donor)")
#plot(identSkinGeno1[validCells],identBlood[validCells],
#     xlim=c(0,1),ylim=c(0,1), xlab = "skin donor similarity", ylab = "blood similarity", 
#     main = "donor vs. blood identity\n(all cells)")
#plot(identSkinGeno2[ifSkinGeno2 & validCells],identBlood[ifSkinGeno2 & validCells],
#     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
#     main = "donor vs. blood identity\nc(ells genotyped as host)")
#plot(identSkinGeno2[validCells],identBlood[validCells],
#     xlim=c(0,1),ylim=c(0,1), xlab = "skin host similarity", ylab = "blood similarity", 
#     main = "donor vs. blood identity\n(all cells)")
#plot(identSkinGeno1[validCells],identSkinGeno2[ validCells],
#     xlim=c(0,1),ylim=c(0,1), xlab = "skin donor similarity", ylab = "skin host similarity", 
#     main = "donor vs. host identity\n(all cells)")
#dev.off()

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
#date     2019-02-14                  
#
#─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#package     * version date       lib source        
#assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)
#backports     1.1.2   2017-12-13 [1] CRAN (R 3.5.0)
#bindr         0.1.1   2018-03-13 [1] CRAN (R 3.5.0)
#bindrcpp      0.2.2   2018-03-29 [1] CRAN (R 3.5.0)
#callr         3.1.0   2018-12-10 [1] CRAN (R 3.5.0)
#cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)
#colorspace    1.3-2   2016-12-14 [1] CRAN (R 3.5.0)
#crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#desc          1.2.0   2018-05-01 [1] CRAN (R 3.5.0)
#devtools      2.0.1   2018-10-26 [1] CRAN (R 3.5.1)
#digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.0)
#dplyr         0.7.8   2018-11-10 [1] CRAN (R 3.5.0)
#fs            1.2.6   2018-08-23 [1] CRAN (R 3.5.0)
#ggplot2     * 3.1.0   2018-10-25 [1] CRAN (R 3.5.0)
#glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.0)
#gtable        0.2.0   2016-02-26 [1] CRAN (R 3.5.0)
#jsonlite      1.6     2018-12-07 [1] CRAN (R 3.5.0)
#labeling      0.3     2014-08-23 [1] CRAN (R 3.5.0)
#lattice       0.20-38 2018-11-04 [1] CRAN (R 3.5.0)
#lazyeval      0.2.1   2017-10-29 [1] CRAN (R 3.5.0)
#magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#Matrix        1.2-15  2018-11-01 [1] CRAN (R 3.5.1)
#memoise       1.1.0   2017-04-21 [1] CRAN (R 3.5.0)
#munsell       0.5.0   2018-06-12 [1] CRAN (R 3.5.0)
#openxlsx    * 4.1.0   2018-05-26 [1] CRAN (R 3.5.0)
#pillar        1.3.0   2018-07-14 [1] CRAN (R 3.5.0)
#pkgbuild      1.0.2   2018-10-16 [1] CRAN (R 3.5.0)
#pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.0)
#pkgload       1.0.2   2018-10-29 [1] CRAN (R 3.5.0)
#plyr          1.8.4   2016-06-08 [1] CRAN (R 3.5.0)
#prettyunits   1.0.2   2015-07-13 [1] CRAN (R 3.5.0)
#processx      3.2.1   2018-12-05 [1] CRAN (R 3.5.1)
#ps            1.2.1   2018-11-06 [1] CRAN (R 3.5.1)
#purrr         0.2.5   2018-05-29 [1] CRAN (R 3.5.0)
#R6            2.3.0   2018-10-04 [1] CRAN (R 3.5.0)
#Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.0)
#remotes       2.0.2   2018-10-30 [1] CRAN (R 3.5.0)
#reshape2      1.4.3   2017-12-11 [1] CRAN (R 3.5.0)
#rjson         0.2.20  2018-06-08 [1] CRAN (R 3.5.0)
#rlang         0.3.0.1 2018-10-25 [1] CRAN (R 3.5.0)
#rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.5.0)
#rstudioapi    0.8     2018-10-02 [1] CRAN (R 3.5.0)
#scales        1.0.0   2018-08-09 [1] CRAN (R 3.5.0)
#sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#stringi       1.2.4   2018-07-20 [1] CRAN (R 3.5.0)
#stringr       1.3.1   2018-05-10 [1] CRAN (R 3.5.0)
#tibble        1.4.2   2018-01-22 [1] CRAN (R 3.5.0)
#tidyselect    0.2.5   2018-10-11 [1] CRAN (R 3.5.0)
#usethis       1.4.0   2018-08-14 [1] CRAN (R 3.5.0)
#withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)
#yaml          2.2.0   2018-07-25 [1] CRAN (R 3.5.0)
#zip           1.0.0   2017-04-25 [1] CRAN (R 3.5.0)