#===================================================================================================================================#
# Title: Bulk RNAseq CD8+ CD45RA- T cells from healthy donors treated with high salt compared to low salt preprocessing and analysis
# Author: Mahima Arunkumar
# Figures: Supplementary Figure 5
#===================================================================================================================================#

#===================================================================================================================================#
# Prepare R session
#===================================================================================================================================#

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

experimentid <- c('EX0009')
procdatatid <- c('GA_PD0009')
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])
lfc.cutoff <- 0

an.desc <- c('cd8_salt')
an.descs <- c('cd8_salt')

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

#===================================================================================================================================#
# Process with DESeq2
#===================================================================================================================================#

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
  rm(txi)
  # Perform standard differential expression analysis using DESeq supplying raw read counts only
  dds[[i]] <- DESeq(dds[[i]])
  # Extract differential expression analysis using a two-group comparison
  res[[i]] <- results(dds[[i]], alpha = 0.05, lfcThreshold = lfc.cutoff, 
                      cooksCutoff = F, contrast = c('condition', cn[[i]]$tr, cn[[i]]$ct))
  # Perform shrinkage of effect size in parallel for better ranking and visualization of fold changes across groups
  resh[[i]] <- lfcShrink(dds[[i]], lfcThreshold = lfc.cutoff, 
                         coef = resultsNames(dds[[i]])[grep('condition', resultsNames(dds[[i]]))], type = 'apeglm')
  
  # Extract differential expression table without and with shrinkage estimation for dispersion
  de[[i]] <- as.data.frame(res[[i]])
  de[[i]]$padj <- ifelse(is.na(de[[i]]$padj), 1, res[[i]]$padj)
  desh[[i]] <- as.data.frame(resh[[i]])
  desh[[i]]$padj <- if(lfc.cutoff == 0) ifelse(
    is.na(desh[[i]]$padj), 1, desh[[i]]$padj) else ifelse(is.na(desh[[i]]$svalue), 1, desh[[i]]$svalue)
  
  nrc[[i]] <- counts(dds[[i]], normalized = T)
  rl[[i]] <- rlog(dds[[i]], blind = F)
}

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

# Generate shrunken differential expression table for all contrasts

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

#PCA
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
ded <- as.data.frame(cbind(up = length(genes$up), 
                           down = -length(genes$dn)))
ded$cont <- 'CD8'
ded <- melt(ded, id.vars = "cont")
ded$type <- 'CD8'

# Plot data
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
