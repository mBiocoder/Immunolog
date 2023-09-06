#=================================================================================================================#
# Exploratory analysis of differential expression between CD8+ CD45RA- T cells from healthy donors treated with 
# high salt compared to low salt concentrations
# Date: 18.01.2023
# Rscript purpose:
# - differential expression between the different conditions from raw read counts data (star mapping and salmon 
#   counting) with DESeq2
# - an overview of differential expression between the different conditions
# - clustering according to differential expression between the different conditions
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(stringr)
library(openxlsx)
library(reshape2)
library(plyr)
library(DESeq2)
library(tximport)
library(ggplot2)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)
library(gridExtra)
library(grid)
library(biomaRt)
library(stringr)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(pathview)
library(DT)
library(cowplot)
#library(bsselectR) #dropdown menu for named vectors

# Define custom theme for ggplot
theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = 'black'),
                      axis.text.y = element_text(size = 12.8, color = 'black'), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = 'black'), 
                      axis.ticks.y = element_line(size = 0.4, colour = 'black'),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour='black'),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = 'black', size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = 'bold', margin = margin(0,0,10,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,'cm'),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid <- c('EX0009')
procdatatid <- c('GA_PD0009')

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])


# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('cd8_salt')
an.descs <- c('cd8_salt')

#=================================================================================================================#
# Process data
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx('./ngs_sample_list.xlsx', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid, ]
#-----------------------------------------------------------------------------------------------------------------#
# Load raw read counts
#-----------------------------------------------------------------------------------------------------------------#

# Load transcript to gene map
tx2gene <- read.csv('tx2gene_gencode_pri.v32.annotation.csv')
# Remove version from gene and transcript names
tx2gene$TXNAME <- str_split(tx2gene$TXNAME, '\\.', n = 2, simplify = T)[, 1]
tx2gene$GENEID <- str_split(tx2gene$GENEID, '\\.', n = 2, simplify = T)[, 1]

# Defines paths for quantification data
pathdir <- list.dirs(paste0('.'))[-1]
pathfile <- list.files(pathdir, pattern = 'quant.sf')
path <- file.path(pathdir, pathfile)
names(path) <- nsl$Sample_name

#-----------------------------------------------------------------------------------------------------------------#
# Process raw read counts with DESeq2
#-----------------------------------------------------------------------------------------------------------------#

# Create sample information table
md <- data.frame(row.names = nsl$Sample_name, 
                 colsplit(nsl$Sample_name, '_', c('type', 'treatment', 'replicate')), 
                 condition = sub('_\\d', '', nsl$Sample_name), stringsAsFactors = T)


# Create list of contrasts
cn <- list(a = list(tr = levels(md$condition)[1], ct = levels(md$condition)[2]))
for(i in 1:length(cn)){
  names(cn)[i] <- paste(cn[[i]]$tr, 'vs', cn[[i]]$ct, sep = '_')
}

# Create lists to store data
dds <- sapply(names(cn), function(x) NULL)
res <- sapply(names(cn), function(x) NULL)
resh <- sapply(names(cn), function(x) NULL)
de <- sapply(names(cn), function(x) NULL)
desh <- sapply(names(cn), function(x) NULL)
#dds.rc <- sapply(names(cn), function(x) NULL)
nrc <- sapply(names(cn), function(x) NULL)
rl <- sapply(names(cn), function(x) NULL)
rl.bc <- sapply(names(cn), function(x) NULL)

# Loop over all contrasts
for(i in 1:length(cn)) {
  
  # Filter file paths
  pathf <- path[grep(paste0(cn[[i]]$ct, '|', cn[[i]]$tr), names(path))]
  
  # Import data with tximport
  txi <- tximport(pathf, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = T, ignoreAfterBar = T)
  
  # Filter metadata
  metadata <- md[rownames(md) %in% nsl[grep(paste0(cn[[i]]$ct, '|', cn[[i]]$tr), nsl$Sample_name), 'Sample_name'], ]
  
  # Choose a reference level for factors (comparison will be performed against the first level)
  metadata$type <- factor(metadata$type, levels = levels(factor(metadata$type))[c(1,2)])
  metadata$treatment <- factor(metadata$treatment, levels = levels(factor(metadata$treatment))[c(2,1)])
  metadata$replicate <- as.factor(metadata$replicate)
  metadata <- metadata[grep(paste0(cn[[i]]$ct, '|', cn[[i]]$tr), rownames(metadata)), ]
  metadata$condition <- relevel(factor(as.character(metadata$condition)), ref = cn[[i]]$ct)
  
  # Construct a DESeqDataSet containing all relevant information about the experiment
  dds[[i]] <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ replicate + condition)
  # Filter transcripts with with a minimum read count per million reads in at least one sample
  #dds <- dds[rowSums(cpm(counts(dds)) >= 0.5) >= 1, ]
  rm(txi)
  
  # Perform standard differential expression analysis using DESeq supplying raw read counts only
  dds[[i]] <- DESeq(dds[[i]])
  
  # Extract differential expression analysis using a two-group comparison
  res[[i]] <- results(dds[[i]], alpha = 0.05, lfcThreshold = lfc.cutoff, 
                      cooksCutoff = F, contrast = c('condition', cn[[i]]$tr, cn[[i]]$ct))
  # Perform shrinkage of effect size in parallel for better ranking and visualization of fold changes across groups
  resh[[i]] <- lfcShrink(dds[[i]], lfcThreshold = lfc.cutoff, 
                         coef = resultsNames(dds[[i]])[grep('condition', resultsNames(dds[[i]]))], type = 'apeglm')
  
  # Check statistics
  #summary(res)
  #summary(resh)
  # Extract differential expression table without and with shrinkage estimation for dispersion
  de[[i]] <- as.data.frame(res[[i]])
  de[[i]]$padj <- ifelse(is.na(de[[i]]$padj), 1, res[[i]]$padj)
  desh[[i]] <- as.data.frame(resh[[i]])
  desh[[i]]$padj <- if(lfc.cutoff == 0) ifelse(
    is.na(desh[[i]]$padj), 1, desh[[i]]$padj) else ifelse(is.na(desh[[i]]$svalue), 1, desh[[i]]$svalue)
  
  # Extract individual read counts for each sample (non-normalized and size-factor-normalized)
  #rc[[i]] <- counts(dds[[i]])
  #for(i in 1:length(cn)) { nrc[[i]] <- counts(dds[[i]], normalized = T) }
  nrc[[i]] <- counts(dds[[i]], normalized = T)
  #colnames(nrc[[i]]) <- sub('^', 'n', colnames(nrc[[i]]))
  
  # Transform data using log2 (and pseudocount 1) and rlog methods for visualization and clustering purposes
  #nt[[i]] <- normTransform(dds[[i]])
  #nt.rc[[i]] <- as.data.frame(assay(nt[[i]]))
  rl[[i]] <- rlog(dds[[i]], blind = F)
  #colnames(assay(rl[[i]])) <- sub('^', 'r', colnames(assay(rl[[i]])))
  #rl.rc <- as.data.frame(assay(rl[[i]]))
  
  # Remove batch effects due to paired samples (treated and untreated from same donor)
  #rl.bc[[i]] <- rl[[i]]
  #assay(rl.bc[[i]]) <- limma::removeBatchEffect(assay(rl.bc[[i]]), metadata$replicate)
  #limma::removeBatchEffect(assay(rl[[i]]), c(1,1,2,1,1,1))
}

#-----------------------------------------------------------------------------------------------------------------#
# Generate differential expression table for all contrasts
#-----------------------------------------------------------------------------------------------------------------#

# Filter for required values and categorize by sample
def <- de
for(i in 1:length(names(de))) {
  def[[i]][, 'sample'] <- names(de)[i]
  def[[i]][, 'geneid'] <- rownames(def[[i]])
  row.names(def[[i]]) <- 1:nrow(def[[i]])
  def[[i]] <- def[[i]][, c(8,1,2,6,3,7)]
}

# Merge data
dem <- def[[1]]

# Categorize according to significant differential expression calculated by DESeq2
dem[, 'DE'] <- c('unchanged')
dem[dem$log2FoldChange > lfc.cutoff & dem$padj < 0.05, 'DE'] <- c('upregulated')
dem[dem$log2FoldChange < lfc.cutoff & dem$padj < 0.05, 'DE'] <- c('downregulated')
colnames(dem)[2:5] <- c('meanExpression', 'lfc', 'fdr', 'se')

# Load database for conversion between different gene id types (generated with biomaRt)
genemap <- read.csv('./genemap_20191025.csv')
# Annotate genes according to Ensembl gene IDs
idx <- match(dem$geneid, genemap$ensembl_gene_id)
dem[, 'genename'] <- genemap$external_gene_name[idx]
dem[, 'description'] <- genemap$description[idx]


# Rearrange table
l <- length(names(cn))
der <- data.frame(dem[, c(1,8,9,2)], 
                  dcast(dem, geneid ~ sample, value.var = c('lfc'), fill = NA)[-1])
colnames(der)[-c(1:4)] <- sub('*', 'LFC_', colnames(der)[-c(1:4)])
der <- data.frame(der, dcast(dem, geneid ~ sample, value.var = c('DE'), fill = NA)[-1])
colnames(der)[-c(1:(l+4))] <- sub('*', 'DE_', colnames(der)[-c(1:(l+4))])
der <- data.frame(der, dcast(dem, geneid ~ sample, value.var = c('fdr'), fill = NA)[-1])
colnames(der)[-c(1:(2*l+4))] <- sub('*', 'FDR_', colnames(der)[-c(1:(2*l+4))])
der <- data.frame(der, dcast(dem, geneid ~ sample, value.var = c('se'), fill = NA)[-1])
colnames(der)[-c(1:(3*l+4))] <- sub('*', 'SE_', colnames(der)[-c(1:(3*l+4))])

# Annotate genes
dea <- de
for(i in 1:length(names(de))) {
  idxa <- match(rownames(dea[[i]]), genemap$ensembl_gene_id)
  dea[[i]][, 'genename'] <- genemap$external_gene_name[idxa]
  dea[[i]][, 'description'] <- genemap$description[idxa]
}
desha <- desh
for(i in 1:length(names(desh))) {
  idxs <- match(rownames(desha[[i]]), genemap$ensembl_gene_id)
  desha[[i]][, 'genename'] <- genemap$external_gene_name[idxs]
  desha[[i]][, 'description'] <- genemap$description[idxs]
}

#-----------------------------------------------------------------------------------------------------------------#
# Generate shrunken differential expression table for all contrasts
#-----------------------------------------------------------------------------------------------------------------#

# Filter for required values and categorize by sample
deshf <- desh
for(i in 1:length(names(desh))) {
  deshf[[i]][, 'sample'] <- names(desh)[i]
  deshf[[i]][, 'geneid'] <- rownames(deshf[[i]])
  row.names(deshf[[i]]) <- 1:nrow(deshf[[i]])
  deshf[[i]] <- deshf[[i]][, c(7,1,2,5,3,6)]
}

# Merge data
deshm <- deshf[[1]]

# Categorize according to significant differential expression calculated by deshSeq2
deshm[, 'DE'] <- c('unchanged')
deshm[deshm$log2FoldChange > lfc.cutoff & deshm$padj < 0.05, 'DE'] <- c('upregulated')
deshm[deshm$log2FoldChange < lfc.cutoff & deshm$padj < 0.05, 'DE'] <- c('downregulated')
colnames(deshm)[2:5] <- c('meanExpression', 'lfc', 'fdr', 'se')

# Annotate genes according to Ensembl gene IDs
deshm[, 'genename'] <- genemap$external_gene_name[match(deshm$geneid, genemap$ensembl_gene_id)]
deshm[, 'description'] <- genemap$description[match(deshm$geneid, genemap$ensembl_gene_id)]

# Rearrange table
deshr <- data.frame(deshm[, c(1,8,9,2)], 
                    dcast(deshm, geneid ~ sample, value.var = c('lfc'), fill = NA)[-1])
colnames(deshr)[-c(1:4)] <- sub('*', 'LFC_', colnames(deshr)[-c(1:4)])
deshr <- data.frame(deshr, dcast(deshm, geneid ~ sample, value.var = c('DE'), fill = NA)[-1])
colnames(deshr)[-c(1:(l+4))] <- sub('*', 'DE_', colnames(deshr)[-c(1:(l+4))])
deshr <- data.frame(deshr, dcast(deshm, geneid ~ sample, value.var = c('fdr'), fill = NA)[-1])
colnames(deshr)[-c(1:(2*l+4))] <- sub('*', 'FDR_', colnames(deshr)[-c(1:(2*l+4))])
deshr <- data.frame(deshr, dcast(deshm, geneid ~ sample, value.var = c('se'), fill = NA)[-1])
colnames(deshr)[-c(1:(3*l+4))] <- sub('*', 'SE_', colnames(deshr)[-c(1:(3*l+4))])

#=================================================================================================================#
# Overview of the data
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# PCA plot
#-----------------------------------------------------------------------------------------------------------------#

# Dimensionality reduction by Principal Component Analysis (PCA) for visualization of overall effects of 
# treatment and batch effects
## Remove zero variance rows PCA analysis
#pca <- limma::removeBatchEffect(assay(rl[[i]]), c(1,1,2,1,1,1))[apply(limma::removeBatchEffect(assay(rl[[i]]), c(1,1,2,1,1,1)), 1, var) != 0, ]
pca <- assay(rl[[1]])[apply(assay(rl[[1]]), 1, var) != 0, ]

## Filter top 500 genes with highest variance
pca <- pca[order(apply(pca, 1, var), decreasing = T), ][1:500, ]
pcad <- as.data.frame(prcomp(t(pca), scale. = T)$x)
pvar <- round(100 * summary(prcomp(t(pca), scale. = T))$importance[2, ])
pcad[, c("type", "treatment", "replicate")] <- colsplit(rownames(pcad), "_", c("type", "treatment", "replicate"))
pcad[, "condition"] <- gsub("_", " ", sub("_\\d$", "", rownames(pcad)))
pcad$condition <- factor(pcad$condition, levels = levels(factor(pcad$condition)))

pcad <- pcad[order(pcad$condition, decreasing = F), ]

limx.pcad <- c(1.1*min(pcad$PC1), 1.2*max(pcad$PC1))
limy.pcad <- c(1.1*min(pcad$PC2), 1.2*max(pcad$PC2))
dot.pca <- ggplot(pcad, aes(PC1, PC2, color = treatment, fill = treatment, label = replicate, shape = treatment)) + 
  geom_point(size = 8, alpha = 0.6) + 
  geom_text(size = 4, color = "black") + 
  scale_x_continuous(paste0("PC1: ", pvar[1], "% variance"), limits = limx.pcad, 
                     breaks = seq(-16,16,8), expand = c(0,0)) +
  scale_y_continuous(paste0("PC2: ", pvar[2], "% variance"), limits = limy.pcad, 
                     breaks = seq(-8,12,4), expand = c(0,0)) + 
  scale_shape_manual("Treatment:", values = c(21,21), 
                     guide = guide_legend(override.aes = list(size = 8, alpha = 0.6, 
                                                              fill = c('grey40', 'white'), 
                                                              color = c('black', 'black')))) + 
  scale_color_manual(values = c('black', 'black'), guide = F) + 
  scale_fill_manual(values = c('grey40', 'white'), guide = F) + 
  ggtitle("PCA (top 500 DE genes)", subtitle = "rlog(read counts)") + theme_custom
dot.pca

#-----------------------------------------------------------------------------------------------------------------#
# Differentially expressed genes
#-----------------------------------------------------------------------------------------------------------------#

# List of differentially expressed genes in each contrast
#for(i in 1:length(names(de))) { de[[i]][, "geneid"] <- rownames(de[[i]]) }
genes <- list(th17 = list(de = sapply(sub("vs_", "vs\n", names(de$tm)), function(x) NULL), 
                          up = sapply(sub("vs_", "vs\n", names(de$tm)), function(x) NULL), 
                          dn = sapply(sub("vs_", "vs\n", names(de$tm)), function(x) NULL)))
#for(i in 1:length(de)) { for(j in 1:length(de[[i]])) {
  genes[[i]]$de[j] <- list(rownames(subset(de[[i]][[j]], de[[i]][[j]]$padj < 0.05)))
  genes[[i]]$up[j] <- list(rownames(subset(de[[i]][[j]], de[[i]][[j]]$padj < 0.05 & 
                                             de[[i]][[j]]$log2FoldChange > lfc.cutoff)))
  genes[[i]]$dn[j] <- list(rownames(subset(de[[i]][[j]], de[[i]][[j]]$padj < 0.05 & 
                                             de[[i]][[j]]$log2FoldChange < lfc.cutoff)))
}}
genes <- list(de = subset(der, der$FDR_cd8_highsalt_vs_cd8_lowsalt < 0.05)$geneid, 
              up = subset(der, der$FDR_cd8_highsalt_vs_cd8_lowsalt < 0.05 & 
                            der$LFC_cd8_highsalt_vs_cd8_lowsalt > lfc.cutoff)$geneid,
              dn = subset(der, der$FDR_cd8_highsalt_vs_cd8_lowsalt < 0.05 & 
                            der$LFC_cd8_highsalt_vs_cd8_lowsalt < -lfc.cutoff)$geneid)

# Plot differential expression
#ded <- as.data.frame(cbind(upregulated = unlist(lapply(genes$up, length)), 
#                           downregulated = -unlist(lapply(genes$dn, length))))
ded <- as.data.frame(cbind(up = length(genes$up), 
                           down = -length(genes$dn)))
#rownames(ded) <- names(de[[1]])
#ded[, "cont"] <- gsub("_", " ", gsub("th17_", "", rownames(ded)))
ded$cont <- 'CD8'
ded <- melt(ded, id.vars = "cont")
ded$type <- 'CD8'
#ded[, c("type", "de")] <- str_split(ded$variable, "_", n = 2, simplify = T)
#ded$variable <- gsub("_", " ", ded$variable)
#ded$cont <- sub("RAneg", "RA-", sub("RApos", "RA+", ded$cont))
#ded$variable <- factor(ded$variable, levels = levels(factor(ded$variable))[c(2,4,1,3)])

# Plot data
#bar.ded <- ggplot(ded, aes(cont, value, group = type, fill = de, color = de, alpha = type)) +
#  geom_bar(position = position_dodge(width = 0.7), width = 0.6, stat = "identity") +
#  scale_x_discrete("") + #coord_flip() +
#  scale_y_continuous("Number of differentially\nexpressed genes", limits = c(-2200,1700),
#                     breaks = seq(-2000,1500,500), expand = c(0,0)) +
#  scale_alpha_manual("Type:", values = c(0.2,1)) +
#  scale_fill_manual("State:", values = c("royalblue", "firebrick3")) +
#  scale_color_manual("State:", values = c("royalblue", "firebrick3")) +
#  geom_hline(aes(yintercept = 0), colour = "black", linetype = "solid", size = 0.4) +
#  ggtitle("T cells", subtitle = "") + theme_custom +
#  theme(aspect.ratio = NULL, axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1, vjust = 1))# +
#bar.ded
#bar.ded <- ggplot(ded, aes(cont, value, group = type, fill = variable, color = variable, label = abs(value))) +
bar.ded <- ggplot(ded, aes(cont, value, fill = variable, group = type, color = variable, label = abs(value))) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6, stat = "identity") + 
  geom_text(aes(y = ifelse(value >= 0, (value + 200), (value - 200))), 
            position = position_dodge(width = 0.7), hjust = 0.5, show.legend = F) + 
  #geom_text(position = position_dodge(width = 0.7), hjust = ifelse(ded$value >= 0, -1, 1)) + 
  scale_x_discrete("", limits = unique(ded$cont)[1]) + #coord_flip() +
  scale_y_continuous("Number of significantly \n differentially expressed genes", limits = c(-2400,2400),
                     breaks = seq(-2000,2000,1000), labels = abs(seq(-2000,2000,1000)), expand = c(0,0)) +
  scale_fill_manual("Regulation", values = c("royalblue", "firebrick3", "#C6D2F6", "#E7BCBC"), guide = guide_legend(reverse = TRUE)) +
  scale_color_manual("Regulation", values = c("royalblue", "firebrick3", "royalblue", "firebrick3"), guide = guide_legend(reverse = TRUE)) + 
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "solid", size = 0.4) +
  ggtitle("T cells", subtitle = "Transcriptome statistics") + theme_custom + #coord_flip() + 
  theme(aspect.ratio = NULL, axis.text.x = element_text(size = 12.8, color = "black", angle = 0, hjust = 0.5, vjust = 0.5))# +
bar.ded

# Top 20 up- and down-regulated genes (fold change)
top <- dem
top <- top[!is.na(top$lfc), ]
top <- top[order(top$lfc, decreasing = T), ]
top <- rbind(head(top[top$fdr < 0.05, ], 30), tail(top[top$fdr < 0.05, ], 30))

#-----------------------------------------------------------------------------------------------------------------#
# Clustering by Pearson correlation coefficient
#-----------------------------------------------------------------------------------------------------------------#

# Pairwise correlation using Pearson correlation coefficient (r) to measure the strength of the linear 
# relationship between samples
corpc <- cor(assay(rl$cd8_highsalt_vs_cd8_lowsalt), method = "pearson")
rownames(corpc) <- gsub("_", " ", rownames(corpc))
corpc <- as.dist(1 - corpc)

## Dendrogram of rlog-transformed read counts for each sample, using the “complete" linkage funtion of the 
## unsupervised hierarchical clustering
par(mfrow=c(1,1), mar=c(5,2,4,9))
plot(as.dendrogram(hclust(corpc, method = "complete")), horiz = T, 
     main = "Distance by Pearson correlation")
## Plot dendogram
dend.pcdist <- recordPlot() #plot.new() or grid.newpage() before plotting?
par(mfrow=c(1,1), mar=c(5,4,4,2))
## Generate heatmap of the distance matrix to assess similarities and dissimilarities between samples
m.corpc <- as.matrix(corpc)
heat.pedist <- pheatmap::pheatmap(m.corpc, clustering_distance_rows = corpc, clustering_distance_cols = corpc, 
                                  col = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100), 
                                  main = "", legend_breaks = seq(0,0.0025,0.0005), silent = T)


#-----------------------------------------------------------------------------------------------------------------#
# Prepare list of ranked genes for enrichment analyses
#-----------------------------------------------------------------------------------------------------------------#

# Annotate genes according to Ensembl gene IDs
#ids <- bitr(dem$geneid, fromType='ALIAS', toType='ENSEMBL', OrgDb='org.Hs.eg.db')
dem[, 'genename'] <- genemap$external_gene_name[idx]
dem[, 'ensembl'] <- genemap$ensembl_gene_id[idx]
dem[, 'entrez'] <- genemap$entrezgene_id[idx]
dem[, 'uniprot'] <- genemap$uniprot_gn_id[idx]
dem[, 'description'] <- genemap$description[idx]
dem$description <- sub(' \\[.*','', dem$description)

# Rename empty cells to NA
dem[dem == ''] <- NA
#dem$dem <- factor(as.character(dem$dem), levels = levels(factor(as.character(dem$dem)))[c(2,3,1)])


#-----------------------------------------------------------------------------------------------------------------#
# Gene set enrichment analysis (GSEA) of GO gene sets
#-----------------------------------------------------------------------------------------------------------------#


# Read in data for bar plot
enriched_go_results <- read.table(file = "GA_AN0284_cd8_salt_overrep_ovrep_go.txt", sep = '\t', header = TRUE)

# Separate into up and downregulated and dysregulated sets
enriched_up_go <- enriched_go_results %>% filter(dif == "upregulated" & ontology== "BP")
enriched_down_go <- enriched_go_results %>% filter(dif == "downregulated" & ontology== "BP")
enriched_dys_go <- enriched_go_results %>% filter(dif == "dysregulated" & ontology== "BP")

# Order terms
#upregulated
up_order_levels <- enriched_up_go %>%
  arrange(p.adjust) %>%
  pull(Description)

enriched_up_go$Description <- factor(enriched_up_go$Description, 
                               levels = rev(up_order_levels))

#downregulated
down_order_levels <- enriched_down_go %>%
  arrange(p.adjust) %>%
  pull(Description)

enriched_down_go$Description <- factor(enriched_down_go$Description, 
                                     levels = rev(down_order_levels))


#Create fresh labels for those which are too long and plot with added features
bp_down_labels <- str_to_sentence(enriched_down_go[1:20,]$Description)
bp_down_labels[c(19, 20)] <- c("Activating cell surface receptor signaling pathway", "Regulating cell surface receptor signaling pathway" ) 

tbr2_down_plot <- ggplot(enriched_down_go[1:20, ]) +
  geom_col(aes(x = Description, 
               y = -log10(p.adjust)), 
           fill = "skyblue", 
           color = "black",
           width=0.6, 
           position = position_dodge(width=0.4)) + 
  coord_flip() +
  scale_x_discrete(name = "",
                   position = "top") +
  scale_y_continuous(name = "-Log(P.adjusted)", 
                     limits = c(0, 10)) +
  scale_y_reverse(name = "-Log(P.adjusted)",
                  expand = c(0, 0),
                  breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) +
  ggtitle("Down-regulated genes") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(color = 'black', size = 0.5),
        axis.line.y.right = element_line(colour = 'black', size = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8,
                                    face = "italic"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5)) +
  geom_text(aes(x = Description,
                y = -log10(p.adjust)/2,
                label = bp_down_labels),
            size = 6 / .pt) +
  theme(plot.margin = unit(c(0.3, 0, 0, 0), "cm"))

tbr2_down_plot



# Create the base of the bar plot
#upregulated
ggplot(enriched_up_go[1:20,]) +
  geom_col(aes(x = Description, 
               y = -log10(p.adjust))) + 
  coord_flip()


#Create fresh labels for those which are too long and plot with added features
bp_up_labels <- str_to_sentence(enriched_up_go[1:20,]$Description)
bp_up_labels[c(15, 16)] <- c("Nucleoside biosynthetic process", "Purine compound biosynthetic process" ) 

tbr2_up_plot <- ggplot(enriched_up_go[1:20,]) +
  geom_col(aes(x = Description, 
               y = -log10(p.adjust)), 
           fill = "lightcoral", 
           color = "black",
           width=0.6, 
           position = position_dodge(width=0.4)) + 
  coord_flip() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "-Log(P.adjusted)", 
                     limits = c(0, 25),
                     expand = c(0, 0)) +
  ggtitle("Up-regulated genes") +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black', 
                                 size = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 8,
                                    face = "italic"),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 8 , hjust = 0.5)) +
  geom_text(aes(x = Description,
                y = -log10(p.adjust)/2,
                label = bp_up_labels),
            size = 6 / .pt) +
  theme(plot.margin = unit(c(0.3, 1.5, 0, 0), "cm"))

tbr2_up_plot


#Combined plot for single enrichment plot showing up and downregulated pathways

# Tbr2 grid
tbr2 <- plot_grid(tbr2_down_plot,
                  tbr2_up_plot)
tbr2


# Adding annotations
annotated_tbr2 <- tbr2 +
  draw_label(label = "Number of \n genes", 
             x = 0.5,
             y = 0.96,
             size = 8) +
  draw_label(label = "46", 
             x = 0.49,
             y = 0.925,
             size = 7) +
  draw_label(label = "44", 
             x = 0.49,
             y = 0.89,
             size = 7) +
  draw_label(label = "45", 
             x = 0.49,
             y = 0.84,
             size = 7) +
  draw_label(label = "45", 
             x = 0.49,
             y = 0.795,
             size = 7) +
  draw_label(label = "59", 
             x = 0.49,
             y = 0.75,
             size = 7) +
  draw_label(label = "44", 
             x = 0.49,
             y = 0.701,
             size = 7) +
  draw_label(label = "93", 
             x = 0.49,
             y = 0.66,
             size = 7)  +
  draw_label(label = "94", 
             x = 0.49,
             y = 0.615,
             size = 7) +
  draw_label(label = "97", 
             x = 0.49,
             y = 0.57,
             size = 7) +
  draw_label(label = "46", 
             x = 0.49,
             y = 0.525,
             size = 7) +
  draw_label(label = "56", 
             x = 0.49,
             y = 0.48,
             size = 7) +
  draw_label(label = "77", 
             x = 0.49,
             y = 0.435,
             size = 7) +
  draw_label(label = "71", 
             x = 0.49,
             y = 0.385,
             size = 7) +
  draw_label(label = "94", 
             x = 0.49,
             y = 0.341,
             size = 7)  +
  draw_label(label = "39", 
             x = 0.49,
             y = 0.30,
             size = 7)  +
  draw_label(label = "92", 
             x = 0.49,
             y = 0.255,
             size = 7)  +
  draw_label(label = "87", 
             x = 0.49,
             y = 0.21,
             size = 7)  +
  draw_label(label = "65", 
             x = 0.49,
             y = 0.165,
             size = 7)  +
  draw_label(label = "62", 
             x = 0.49,
             y = 0.113,
             size = 7)  +
  draw_label(label = "66", 
             x = 0.49,
             y = 0.07,
             size = 7)  +
  draw_label(label = "97", 
             x = 0.51,
             y = 0.925,
             size = 7) +
  draw_label(label = "66", 
             x = 0.51,
             y = 0.89,
             size = 7) +
  draw_label(label = "95", 
             x = 0.51,
             y = 0.84,
             size = 7) +
  draw_label(label = "82", 
             x = 0.51,
             y = 0.795,
             size = 7) +
  draw_label(label = "120", 
             x = 0.51,
             y = 0.75,
             size = 7) +
  draw_label(label = "103", 
             x = 0.51,
             y = 0.701,
             size = 7) +
  draw_label(label = "71", 
             x = 0.51,
             y = 0.66,
             size = 7)  +
  draw_label(label = "115", 
             x = 0.51,
             y = 0.615,
             size = 7) +
  draw_label(label = "67", 
             x = 0.51,
             y = 0.57,
             size = 7) +
  draw_label(label = "30", 
             x = 0.51,
             y = 0.525,
             size = 7) +
  draw_label(label = "67", 
             x = 0.51,
             y = 0.48,
             size = 7) +
  draw_label(label = "26", 
             x = 0.51,
             y = 0.435,
             size = 7) +
  draw_label(label = "60", 
             x = 0.51,
             y = 0.385,
             size = 7) +
  draw_label(label = "78", 
             x = 0.51,
             y = 0.341,
             size = 7)  +
  draw_label(label = "78", 
             x = 0.51,
             y = 0.30,
             size = 7)  +
  draw_label(label = "65", 
             x = 0.51,
             y = 0.255,
             size = 7)  +
  draw_label(label = "48", 
             x = 0.51,
             y = 0.21,
             size = 7)  +
  draw_label(label = "63", 
             x = 0.51,
             y = 0.165,
             size = 7)  +
  draw_label(label = "86", 
             x = 0.51,
             y = 0.113,
             size = 7)  +
  draw_label(label = "47", 
             x = 0.51,
             y = 0.07,
             size = 7)  +
  draw_label(label = "Over???representation: GO term (Biological process)", 
             x = 0.5,
             y = 1,
             hjust = 0.5,
             vjust = 1,
             size = 9) 

annotated_tbr2

# Saving enrichment plot as figure
ggsave(plot = annotated_tbr2,
       device = "pdf",
       filename = "results/enrichment_go_result_up_and_down.pdf",
       width = 3.5,
       height = 2.25,
       units = "in",
       dpi = 300)


# Make small enrichment plot for dysregulated (unchanged) genes

data <- readRDS("GA_AN0284_cd8_salt_overrep_go.rds")

tbr2_dys_order_levels <- enriched_dys_go %>%
  arrange(p.adjust) %>%
  pull(Description)

enriched_dys_go$Description <- factor(enriched_dys_go$Description, 
                               levels = rev(tbr2_dys_order_levels))


# Create function to design dotplots
dotplotf <- function(data, color, title) {
  if(sum(data@result$qvalue < 0.1) > 0) {
    dotplot(data, showCategory = 20, color = color, title = title) + theme_custom
  } else ggplot() + ggtitle("No gene sets were significantly enriched!") + 
    theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5, margin = margin(5,0,20,0)))
}

dotplotf(data = data$bp$de, "p.adjust", "Over???representation: GO term
(Biological process)")



############################# Gustavo part ###########################


# Load object with pre-calculated GSEA
gsea <- readRDS("CD8_salt_gsea_all.rds")

# NES = normalized enrichment score (the normalized enrichment score accounts for differences in gene set size and 
# in correlations between gene sets and the expression dataset; therefore, it can be used to compare analysis 
# results across gene sets)
# padj = adjusted p-value that allows for comparison between different signatures when tested simultaneously
# qvalue = adjusted p-value using an opimised FDR approach and represents the chance of false positives between 
# all compounds with smaller or equal q-values (and not between all compounds independently of p-value)
# https://support.bioconductor.org/p/49864/
# http://genomics.princeton.edu/storeylab/papers/directfdr.pdf


# Create function to design dotplots
dotplotf <- function(data, color, title) {
  if(sum(data@result$qvalue < 0.1) > 0) {
    dotplot(data, showCategory = 50, color = color, title = title) + facet_grid(enriched ~ .) + theme_custom
  } else ggplot() + ggtitle("No gene sets were significantly enriched!") + 
    theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5, margin = margin(5,0,20,0)))
}

# # Calculate over-representation of GO gene sets over the list of genes ranked according to differential expression
# gsea <- list()
# gsea$go <- gseGO(geneList = rank.de$en, OrgDb = org.Hs.eg.db, ont = "ALL", nPerm = 10000, minGSSize = 10, 
#                  maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 1, verbose = FALSE, keyType = "ENSEMBL")
# gsea$go <- setReadable(gsea$go, 'org.Hs.eg.db', keytype = "ENSEMBL")
# gsea$go@result[, "enriched"] <- ifelse(gsea$go@result$NES > 0, tr, ct)
# colnames(gsea$go@result)[grep("qvalues", colnames(gsea$go@result))] <- c("qvalue")

# Filter data for qvalue < 0.1
gsea$go.cc.ct <- gsea$go.mf.ct <- gsea$go.bp.ct <- gsea$go.ct <- gsea$go
gsea$go.ct@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1, ]
gsea$go.bp.ct@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'BP', ]
gsea$go.mf.ct@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'MF', ]
gsea$go.cc.ct@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'CC', ]
gsea$go.cc.tr <- gsea$go.mf.tr <- gsea$go.bp.tr <- gsea$go.tr <- gsea$go
gsea$go.tr@result <- gsea$go@result[gsea$go@result$enriched == tr & gsea$go@result$qvalue < 0.1, ]
gsea$go.bp.tr@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'BP', ]
gsea$go.mf.tr@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'MF', ]
gsea$go.cc.tr@result <- gsea$go@result[gsea$go@result$enriched == ct & gsea$go@result$qvalue < 0.1 & 
                                         gsea$go@result$ONTOLOGY == 'CC', ]

# Plot results
dot <- list()
dot$gsea.go.ct <- dotplotf(gsea$go.ct, color = "qvalue", title = "GSEA: GO (All ontologies)")
dot$gsea.go.bp.ct <- dotplotf(gsea$go.bp.ct, color = "qvalue", title = "GSEA: GO (Biological Processes)")
dot$gsea.go.mf.ct <- dotplotf(gsea$go.mf.ct, color = "qvalue", title = "GSEA: GO (Molecular Functions)")
dot$gsea.go.cc.ct <- dotplotf(gsea$go.cc.ct, color = "qvalue", title = "GSEA: GO (Cellular Compartments)")
dot$gsea.go.tr <- dotplotf(gsea$go.tr, color = "qvalue", title = "GSEA: GO (All ontologies)")
dot$gsea.go.bp.tr <- dotplotf(gsea$go.bp.tr, color = "qvalue", title = "GSEA: GO (Biological Processes)")
dot$gsea.go.mf.tr <- dotplotf(gsea$go.mf.tr, color = "qvalue", title = "GSEA: GO (Molecular Functions)")
dot$gsea.go.cc.tr <- dotplotf(gsea$go.cc.tr, color = "qvalue", title = "GSEA: GO (Cellular Compartments)")

# Extract enrichment results to display table containing the data
table <- list()
table$go.all <- dotplot(gsea$go, showCategory = length(gsea$go@geneSets))$data[, c(13,1,3,14,4,16,6,7,8,9,2,12)]
table$go.all <- table$go.all[order(table$go.all$pvalue, decreasing = F), ]
table$go.all <- table$go.all[order(abs(table$go.all$NES), decreasing = T), ]
rownames(table$go.all) <- seq(length = nrow(table$go.all))
table$go.allf <- table$go.all[table$go.all$qvalue <= 0.1, ]
table$go.allf[sapply(table$go.allf, is.numeric)] <- lapply(table$go.allf[sapply(table$go.allf, is.numeric)], round, 3)
rownames(table$go.allf) <- seq(length = nrow(table$go.allf))
table$go.allf[, c(1,2,3,11)] <- lapply(table$go.allf[, c(1,2,3,11)], as.factor)



#-----------------------------------------------------------------------------------------------------------------#
# Save data and results
#-----------------------------------------------------------------------------------------------------------------#

# Save tables with differential expression
write.csv(der, file = "./output_files/diff_expression.csv", row.names = F)
write.table(der, file = "./output_files/diff_expression.txt", sep = '\t', row.names = F)
write.table(top[,c(1,8,9,2,3,4,5,7)], file = "./output_files/top30_fc.txt", sep = "\t", row.names = F)

# Save tables with differential expression
write.csv(deshr, file = './output_files/shdiff_expression.csv', row.names = F)

# Save rlog size factor normalized absolute expression
for(i in 1:length(names(rl))) {
  write.csv(assay(rl[[i]]), file = paste0(dir.figsr, '_', names(rl)[i], '_rlog_expression.csv'), row.names = T)
}
## Batch corrected
#for(i in 1:length(cn)) {
#  write.csv(assay(rl.bc[[i]]), file = paste0(dir.figsb, '_', names(rl.bc)[i], '_rlog_expression_bc.csv'), row.names = T)
#}

# Save size factor normalized absolute expression
for(i in 1:length(names(nrc))) {
  write.csv(nrc[[i]], file = paste0(dir.figsn, '_', names(nrc)[i], '_diff_expression.csv'), row.names = T)
}
# Save differential expression tables for each contrast
for(i in 1:(length(cn))) {
  write.csv(dea[[i]], file = paste0(dir.figsd, '_', names(dea)[i], '_diff_expression.csv'), row.names = T)
}
## Shrunken values
for(i in 1:(length(cn))) {
  write.csv(desha[[i]], file = paste0(dir.figss, '_', names(desha)[i], '_diff_expression_shrunken.csv'), row.names = T)
}

# Save read counts
#write.csv(counts(dds[[i]]), file = paste0(dir.figsc, '_', names(rl)[i], '_read_counts.csv'), row.names = T)
#counts <- read.csv(paste0(dir.figs, '_read_counts.csv'), row.names = 1)
#write.csv(format(assays(dds[[i]])[['avgTxLength']], digits=3), 
#          file = paste0(dir.figsc, '_', names(rl)[i], '_avgTxLength.csv'), row.names = T)
#avgTxLength <- as.matrix(read.csv(paste0(dir.figs, '_avgTxLength.csv'), row.names = 1))
#}

# Save plot
#ggsave(plot = dot.pca.ds, file = paste0(dir.figs.pca, '_pca_rlog.pdf'), height = 4.5, width = 6)
ggsave(plot = dot.pca, file = paste0(dir.figs.pca, '_pca_rlog.pdf'), height = 4.5, width = 6)
ggsave(plot = bar.ded, file = paste0(dir.figs.stat, "_diff_stat.pdf"), height = 5, width = 4)

pdf(paste0(dir.figs.cor, "_dend_pcdist.pdf"), height = 4, width = 6)
dend.pcdist
dev.off()
pdf(paste0(dir.figs.cor, "_heat_pedist.pdf"), height = 5, width = 5)
gridExtra::grid.arrange(grobs = list(heat.pedist[[4]]), nrow = 1, 
                        top = grid::textGrob("Clustering of pairwise distances\nby Pearson correlation", 
                                             gp = grid::gpar(fontsize = 16,font = 2)))
dev.off()
