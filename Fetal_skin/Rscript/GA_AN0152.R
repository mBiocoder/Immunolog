#=================================================================================================================#
# Exploratory analysis of 10x Genomics scRNA-Seq from CD45+ cells isolated from blood and skin of a healthy donor
# Date: 2019.01.22
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze 10x Genomics data processed with Cellranger
# - differential expression of each cluster
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(openxlsx)
library(Seurat)
library(dplyr)
#library(SingleR)
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
experimentid <- c("EX0004")

# Define Analysis_ID according to ngs_analysis_list
analysisid <- c("GA_AN0152")

# Define samples analyzed (all)
#sa <- c("cd")

# Define reference sample name
ct <- c("blood")
tr <- c("skin")

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c("sc_skin_blood_10x_markers")
an.descs <- c("sc_skin")

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx("../../ngs_sample_list.xlsm", sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = "NA", fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid, ]

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

# Copy "script" file to "analysis_old/scripts_old" folder adding run date to its name
file.copy(from = paste0(analysisid, ".R"), 
          to = paste0("analysis_old/", analysisid, format(Sys.time(), "_%Y%m%d_%H%M%S"), ".R"))

# Move "report" file to "analysis_old" folder adding creation date to its name
file.rename(from = paste0("", analysisid, ".html"), 
            to = paste0("analysis_old/", analysisid, 
                        format(file.info(paste0(analysisid, ".html"))$ctime, "_%Y%m%d_%H%M%S"), ".html"))

# Move "figures" folder to "analysis_old/figures_old" folder
if(dir.exists(paste0("figures/", analysisid, "_", an.desc))) {
  file.copy(from = paste0("figures/", analysisid, "_", an.desc), to = paste0("analysis_old/figures_old/"), 
            recursive = T, copy.date = T) }
file.rename(from = paste0("analysis_old/figures_old/", analysisid, "_", an.desc), 
            to = paste0("analysis_old/figures_old/", analysisid, "_", an.desc, 
                        format(file.info(paste0("analysis_old/figures_old/", 
                                                analysisid, "_", an.desc))$ctime, "_%Y%m%d_%H%M%S")))
unlink(paste0("figures/", analysisid, "_", an.desc), recursive = T)

# Create new folder to store figures from this analysis
dir.create(paste0("figures/", analysisid, "_", an.desc))
#dir.create(paste0("figures/", analysisid, "_", an.desc, "/blood"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/tsne"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/umap"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/heatmap"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/violin"))
dir.create(paste0("figures/", analysisid, "_", an.desc, "/skin/doublets"))
dir.figs <- paste0("figures/", analysisid, "_", an.desc, "/", analysisid, "_", an.descs)
#dir.figs.bd <- paste0("figures/", analysisid, "_", an.desc, "/blood/", analysisid, "_", an.desc)
dir.figs.sk <- paste0("figures/", analysisid, "_", an.desc, "/skin/", analysisid, "_", an.descs)
dir.figs.ts <- paste0("figures/", analysisid, "_", an.desc, "/skin/tsne/", analysisid, "_", an.descs)
dir.figs.um <- paste0("figures/", analysisid, "_", an.desc, "/skin/umap/", analysisid, "_", an.descs)
dir.figs.hm <- paste0("figures/", analysisid, "_", an.desc, "/skin/heatmap/", analysisid, "_", an.descs)
dir.figs.vp <- paste0("figures/", analysisid, "_", an.desc, "/skin/violin/", analysisid, "_", an.descs)
dir.figs.db <- paste0("figures/", analysisid, "_", an.desc, "/skin/doublets/", analysisid, "_", an.descs)

#=================================================================================================================#
# Analysis of processed data for skin sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Clustering and differential expression analysis with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Load Seurat object
sk <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_seurat_skin.rds"))
# Load SingleR object
sr.sk <- readRDS(paste0("../201812/figures/GA_AN0146_sc_skin_blood_10x/skin/GA_AN0146_sc_skin_blood_10x_singler_skin.rds"))
# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
sr.sk$meta.data$orig.ident <- sk@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.sk$meta.data$xy <- sk@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "tsne"))
sr.sk$meta.data$xy.um <- sk@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "umap"))
sr.sk$meta.data$clusters <- sk@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(sk))

# Find marker genes for each cluster
#sk.markers <- FindAllMarkers(sk, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
#                             test.use = "wilcox", return.thresh = 0.01, random.seed = 1)
#write.csv(sk.markers, file = paste0(dir.figs.sk, "_cluster_markers.csv"), row.names = T)
sk.markers <- read.csv(file = paste0(dir.figs.sk, "_cluster_markers.csv"), row.names = 1, stringsAsFactors = F)
top <- function(i=-Inf,j=Inf,k){ sk.markers[sk.markers$avg_logFC >= i & sk.markers$avg_logFC <= j,] %>% 
    group_by(cluster) %>% top_n(n = k, wt = abs(avg_logFC)) }
#sk.top1 <- sk.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
sk.top10 <- top(0,Inf,10)
sk.bot10 <- top(-Inf,0,10)
#mk <- c("IL7R", "CD4", "CD8A", "CD8B", "GNLY", "MS4A1", "FCGR3A", "NKG7", "GZMB")
#mk <- c('CD3E','CD4','CD8A', 'CD8B', 'CCR7', 'FCER1A', 'PPBP', 'FCGR3A', 'NKG7', 'GZMB', 'GZMA', 'GNLY', 'MS4A1', 'CD14', 'LYZ', 'CD34')
#treg <- c("CTLA4", "FOXO1", "FOXO3", "FOXP3", "IL2RA", "TIGIT", "ICOS")
#mk <- top(0,Inf,30)[top(0,Inf,30)$cluster == 6,]$gene
#mk <- sk.top10[sk.top10$cluster == 6,]$gene
mk <- top(0,Inf,30)
mn <- top(-Inf,0,30) #negative markers
#write.csv(mk, file = paste0(dir.figs.sk, "_cluster_markers_top30.csv"), row.names = T)
#write.csv(mn, file = paste0(dir.figs.sk, "_cluster_markers_bottom30.csv"), row.names = T)

# Overview of the data
VlnPlot(sk, features = c('CD3E','CD8A', 'CD8B', 'NKG7', 'GZMA', 'GNLY','MS4A1'), cols = NULL, pt.size = 0, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")
FeaturePlot(sk, features = mk[1], dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = NA, reduction = "tsne", blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = F)
hm.mk.sk <- DoHeatmap(sk, features = top(0,Inf,10)$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers (skin)")# + NoLegend()
hm.mn.sk <- DoHeatmap(sk, features = top(-Inf,0,10)$gene, disp.min = -2.5, disp.max = NULL, slot = "scale.data", label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T) + ggtitle("Heatmap of top 10 cluster markers (skin)")# + NoLegend()
DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1", split.by = "orig.ident", label = T, label.size = 8, pt.size = 1, cols = NULL) + theme_bw() + theme_custom
DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1", split.by = "orig.ident", label = T, label.size = 8, pt.size = 1, cols = NULL) + theme_bw() + theme_custom

# Visualization of marker gene expression (tSNE)
ts.mk.sk <- sapply(as.numeric(levels(sk@active.ident))+1, function(x) NULL)
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
df.sk <- data.frame(x = sk@reductions$tsne@cell.embeddings[,1],
                    y = sk@reductions$tsne@cell.embeddings[,2],
                    t(as.matrix(GetAssayData(sk, slot = "data")[mk[mk$cluster == i,]$gene,])))
df.sk <- reshape2::melt(df.sk, id.vars = c('x','y'))
ts.mk.sk[[i+1]] <- ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
  scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + facet_wrap(~variable, ncol = 6) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle(paste0("Expression of top 30 marker genes from cluster ", i))
}

# Visualization of marker gene expression (UMAP)
um.mk.sk <- sapply(as.numeric(levels(sk@active.ident))+1, function(x) NULL)
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
  df.sk <- data.frame(x = sk@reductions$umap@cell.embeddings[,1],
                      y = sk@reductions$umap@cell.embeddings[,2],
                      t(as.matrix(GetAssayData(sk, slot = "data")[mk[mk$cluster == i,]$gene,])))
  df.sk <- reshape2::melt(df.sk, id.vars = c('x','y'))
  um.mk.sk[[i+1]] <- ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
    scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
    scale_x_continuous(breaks = seq(-4,4,4)) + scale_y_continuous(breaks = seq(-4,4,4)) + facet_wrap(~variable, ncol = 6) + 
    theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("Expression of top 30 marker genes from cluster ", i))
}

# Visualization of marker gene expression (tSNE)
vp.mk.sk <- sapply(as.numeric(levels(sk@active.ident))+1, function(x) NULL)
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
  df.sk <- data.frame(cluster = sk@active.ident, 
                      t(as.matrix(GetAssayData(sk, slot = "data")[mk[mk$cluster == i,]$gene,])))
  df.sk <- reshape2::melt(df.sk, id.vars = c("cluster"))
  df.sk$value <- df.sk$value + rnorm(n = length(df.sk$value)) / 100000 #add noise
  vp.mk.sk[[i+1]] <- ggplot(df.sk, aes(x = cluster, y = value, fill = cluster)) + 
    geom_violin(scale = "width", color = "black", size = 0.2, show.legend = F) + 
    #geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2), color = "black", show.legend = F) + 
    scale_x_discrete(breaks = seq(0,12,1), labels = c(rbind(seq(0,12,2), ""))[-14]) + #scale_y_continuous(breaks = seq(-20,20,20)) + 
    xlab('') + ylab('Expression level') + ggtitle(paste0("Expression of top 30 marker genes from cluster ", i)) + 
    facet_wrap(~variable, ncol = 6) + theme_custom
}

s#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for everage expression of top cluster markers
#-----------------------------------------------------------------------------------------------------------------#

# Extract z-scores (scaled data)
hm <- GetAssayData(sk, slot = "scale.data")

# Add cluster information
hm <- reshape2::melt(hm)
hm[, "cluster"] <- hm$Var2
hm$cluster <- plyr::mapvalues(hm$cluster, from = levels(hm$cluster)[match(levels(hm$cluster), names(Idents(sk)))], 
                              to = as.character(Idents(sk)))

# Calculate average expression for each cluster
hm <- plyr::ddply(hm, c("Var1", "cluster"), summarise, mean = mean(value))
#write.csv(hm, file = paste0(dir.figs.hm, "_cluster_expression.csv"), row.names = T)
#hm <- read.csv(file = paste0(dir.figs.hm, "_cluster_expression.csv"), row.names = 1, stringsAsFactors = T)

# Transform data frame
hm$cluster <- as.numeric(as.character(hm$cluster))
hd <- reshape2::dcast(hm, Var1~cluster, value.var = "mean")
rownames(hd) <- hd$Var1
gn <- top(0,Inf,10)$gene
gn <- top(0,Inf,30)[top(0,Inf,30)$cluster == 0,]$gene

# Filter for cluster markers
hf <- hd[row.names(hd) %in% gn, -1]
hf <- hd[match(unique(gn), row.names(hd)), -1]

# Draw heatmap
#ann_row <- as.data.frame(top(0,Inf,10)[, 6:7])
#ann_row <- ann_row[-which(duplicated(ann_row$gene)),]
#ann_row <- as.data.frame(top(0,Inf,30)[top(0,Inf,30)$cluster == 0, 6:7])
#rownames(ann_row) <- ann_row$gene
#ann_row <- ann_row[-2]
#ann_colors <- list(cluster = c(`0` = "black", `1` = "black", `2` = "black", `3` = "black", `4` = "black", 
#                               `5` = "black", `6` = "black", `7` = "black", `8` = "black", `9` = "black", 
#                               `10` = "black", `11` = "black", `12` = "black", `13` = "black", `14` = "black"))
heat.sk <- pheatmap(hf, cluster_rows = T, show_rownames = T, cluster_cols = F, border_color = NA, #scale = "row", 
                    col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                    breaks = seq(-2,2,length.out = 101), 
                    #annotation_row = ann_row, annotation_colors = ann_colors, annotation_names_row = T, 
                    main = "Heatmap of top 10 cluster markers (skin)", #legend_breaks = seq(0,21,3), 
                    silent = F)

# Create function to plot heatmap
heat <- function(x,a,b,c){pheatmap(x, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = F, scale = c, 
                             cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                             border_color = NA, breaks = seq(a,b,length.out = 101), silent = T, 
                             col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100))}
heat.sk <- sapply(0:12, function(x) NULL)
heats.sk <- sapply(0:12, function(x) NULL)
for(i in 0:max(hm$cluster)){
  heat.sk[[i+1]] <- heat(hd[match(unique(top(0,Inf,30)[top(0,Inf,30)$cluster == i,]$gene), row.names(hd)), -1], -2,2,"none")
  heats.sk[[i+1]] <- heat(hd[match(unique(top(0,Inf,30)[top(0,Inf,30)$cluster == i,]$gene), row.names(hd)), -1], -2,2,"row")
}

# Plot heatmaps
grid.arrange(grobs = list(heat.sk[[1]][[4]], heat.sk[[2]][[4]], heat.sk[[3]][[4]], heat.sk[[4]][[4]], heat.sk[[5]][[4]], 
                          heat.sk[[6]][[4]], heat.sk[[7]][[4]], heat.sk[[8]][[4]], heat.sk[[9]][[4]], heat.sk[[10]][[4]], 
                          heat.sk[[11]][[4]], heat.sk[[12]][[4]], heat.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers: z-scores (skin)", gp = gpar(fontsize = 14, font = 2), vjust = 0.6))

grid.arrange(grobs = list(heats.sk[[1]][[4]], heats.sk[[2]][[4]], heats.sk[[3]][[4]], heats.sk[[4]][[4]], heats.sk[[5]][[4]], 
                          heats.sk[[6]][[4]], heats.sk[[7]][[4]], heats.sk[[8]][[4]], heats.sk[[9]][[4]], heats.sk[[10]][[4]], 
                          heats.sk[[11]][[4]], heats.sk[[12]][[4]], heats.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers: row-scaled z-scores (skin)", gp = gpar(fontsize = 14, font = 2), vjust = 0.6))

#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for everage expression of top cluster markers
#-----------------------------------------------------------------------------------------------------------------#

# Extract log transformed normalized data
ht <- as.matrix(GetAssayData(sk, slot = "data"))

# Add cluster information
ht <- reshape2::melt(ht)
ht[, "cluster"] <- ht$Var2
ht$cluster <- plyr::mapvalues(ht$cluster, from = levels(ht$cluster)[match(levels(ht$cluster), names(Idents(sk)))], 
                              to = as.character(Idents(sk)))

# Calculate average expression for each cluster only considering cells expressing the gene
hl <- plyr::ddply(ht, c("Var1", "cluster"), summarise, mean = mean(value))
#write.csv(hl, file = paste0(dir.figs.hm, "_cluster_expression_log.csv"), row.names = T)
#hl <- read.csv(file = paste0(dir.figs.hm, "_cluster_expression_log.csv"), row.names = 1, stringsAsFactors = T)

# Transform data frame
hl$cluster <- as.numeric(as.character(hl$cluster))
he <- reshape2::dcast(hl, Var1~cluster, value.var = "mean")
rownames(he) <- he$Var1

# Calculate average expression for each cluster only considering cells expressing the gene
hn <- plyr::ddply(ht, c("Var1", "cluster"), summarise, mean = sum(value)/sum(value != 0))
#write.csv(hn, file = paste0(dir.figs.hm, "_cluster_expression_log_sub.csv"), row.names = T)
#hn <- read.csv(file = paste0(dir.figs.hm, "_cluster_expression_log_sub.csv"), row.names = 1, stringsAsFactors = T)

# Transform data frame
hn$cluster <- as.numeric(as.character(hn$cluster))
hg <- reshape2::dcast(hn, Var1~cluster, value.var = "mean")
rownames(hg) <- hg$Var1
hg[is.na(hg)] <- 0

# Create function to plot heatmap
heatl <- function(x,a,b){pheatmap(x, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = F, #scale = "row", 
                             cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                             border_color = NA, breaks = seq(a,b,length.out = 101), silent = T, legend_breaks = seq(0,10,1), 
                             col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100))}
heatl.sk <- sapply(0:12, function(x) NULL)
heatn.sk <- sapply(0:12, function(x) NULL)
for(i in 0:max(hl$cluster)){
  heatl.sk[[i+1]] <- heatl(he[match(unique(top(0,Inf,30)[top(0,Inf,30)$cluster == i,]$gene), row.names(he)), -1], 0,3)
  heatn.sk[[i+1]] <- heatl(hg[match(unique(top(0,Inf,30)[top(0,Inf,30)$cluster == i,]$gene), row.names(hg)), -1], 0,4)
}

# Plot heatmaps
grid.arrange(grobs = list(heatl.sk[[1]][[4]], heatl.sk[[2]][[4]], heatl.sk[[3]][[4]], heatl.sk[[4]][[4]], 
                          heatl.sk[[5]][[4]], heatl.sk[[6]][[4]], heatl.sk[[7]][[4]], heatl.sk[[8]][[4]], 
                          heatl.sk[[9]][[4]], heatl.sk[[10]][[4]], heatl.sk[[11]][[4]], heatl.sk[[12]][[4]], 
                          heatl.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers: mean log expression (skin)", 
                            gp = gpar(fontsize = 14, font = 2), vjust = 0.6))
grid.arrange(grobs = list(heatn.sk[[1]][[4]], heatn.sk[[2]][[4]], heatn.sk[[3]][[4]], heatn.sk[[4]][[4]], 
                          heatn.sk[[5]][[4]], heatn.sk[[6]][[4]], heatn.sk[[7]][[4]], heatn.sk[[8]][[4]], 
                          heatn.sk[[9]][[4]], heatn.sk[[10]][[4]], heatn.sk[[11]][[4]], heatn.sk[[12]][[4]], 
                          heatn.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers: mean log expression of expressing cells (skin)", 
                            gp = gpar(fontsize = 14, font = 2), vjust = 0.6))

#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for everage expression of cell type markers
#-----------------------------------------------------------------------------------------------------------------#

# Create list of genes
#monocytes exclusion 'CD14', 'CD14.1', CCR2', 'CX3CR1' (survival)
#B exclusion  'CD19', 'CXCR5',  'CCR6', 'CD20'
#t 'CD3', 'CD4', 'CD8A', 'CD8B'
#mait and rorgtcd4 'CD161' (KLRB1), 
#gd 'TRDC', 'TRGC1', 'TRGC2'
#mem 'CD95' (FAS), 'CD45RO',
#naïve tcm 'CD27', 'CD28',  'CD45RA', 
#temra  'CD45RA', 
#trafficking 'CD29', 'CLA' (SELPLG), 'CX3CR1'
#trm 'CD103' (ITGAE), 'CD69', 'CXCR6', 'CD49a' (ITGA1),
#skin homing 'CCR10' (CCR4 and CCR8 too)
#nk (cd3-) nkt(cd3+) 'CD56', 
#treg 'CD25',  CD278 (ICOS check), 'CD127lo/neg', 
#activated-tcell 'CD25', 'HLA-DR', 'PD-1' (also tex), 
#teff 'CD127lo/neg', 'CXCR3', (and th1)
#naïve and trafficking 'CD62L' (SELL), 'CCR7', 'CXCR4',
#survival/naive/mem 'CD127', 
#naive cd4 young 'CD31', 
#distinguish naive(dn) tcm(-+) tem(dp) temra(+-) 'CD57'+ 'PD1'
#colon recruitment 'CCR9', 
#cutaneous T-cell lymphoma 'CCR4', 
#breast or prostate cancer cells 'CCR5'
# PDCD1 (PD1)(CD279), B3GAT1 (CD57), PECAM1 (CD31), CD127 (IL7R), CD25 (IL2RA), CD45R (PTPRC), CCR7 (CD197), CD107a (LAMP1), 
# KIT (CD117), CD11b (ITGAM), CD11c (ITGAX), CD122 (IL2RB), CD123 (IL3RA), CD141 (THBD), CD154 (CD40LG), CD16 (FCGR3A), 
# CD243 (MDR1, ABCB1), CD271 (NGFR), CD178 (FASLG), CD336 (NCR2), CD49B (ITGA2), CD56 (NCAM1), CRTH2 (PTGDR2), 
# CTLA4 (CD152), EMA (MUC1), GLUT1 (SLC2A1), Hobit (ZNF683), IL8 (CXCL8), Integrin VLA4 (ITGA4), p38a (MAPK14), 
# p38b (MAPK11), PD-L1 (CD274), Podoplanin (PDPN), RORγt (RORC), ST2/IL-33R (IL1RL1), T-bet (TBX21), GMCSF (CSF2), 

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
#               '', '', '', '', '', '', '', '', '', '')

# All TCR genes
#gn <- a[grep("TRG|TRD", a$`row.names(hd)`),1]
#gn %in% row.names(hd)

# Filter data for cluster markers
hd.gn <- hd[match(unique(gn), row.names(hd)), -1]
he.gn <- he[match(unique(gn), row.names(he)), -1]

# Draw heatmap
#ann_row <- as.data.frame(top(0,Inf,10)[, 6:7])
#ann_row <- ann_row[-which(duplicated(ann_row$gene)),]
#ann_row <- as.data.frame(top(0,Inf,30)[top(0,Inf,30)$cluster == 0, 6:7])
#rownames(ann_row) <- ann_row$gene
#ann_row <- ann_row[-2]
#ann_colors <- list(cluster = c(`0` = "black", `1` = "black", `2` = "black", `3` = "black", `4` = "black", 
#                               `5` = "black", `6` = "black", `7` = "black", `8` = "black", `9` = "black", 
#                               `10` = "black", `11` = "black", `12` = "black", `13` = "black", `14` = "black"))

# Draw heatmaps with gene annotation
ann_row <- data.frame(subset = names(gn), row.names = gn)
ann_colors <- list(subset = c('T cell' = "grey40", 
                              'CD4 T cell' = "royalblue2", 
                              'CD8 T cell' = "firebrick1", 
                              'Tn' = "grey70",
                              'DC' = "black", 
                              'Monocyte' = "black", 
                              'NK cell' = "darkorange3", 
                              'B cell' = "seagreen3", 
                              'HSC' = "black", 
                              'Tm' = "darkgoldenrod1", 
                              'Tem' = "steelblue1", 
                              'Treg' = "olivedrab", 
                              'Trm' = "firebrick3", 
                              'trafficking' = "wheat2", 
                              'RORgT cell' = "plum3", 
                              'Tn and Tcm' = "darkorchid2", 
                              'gdT cell' = "grey90"))
heat.sk.hd <- pheatmap(hd.gn, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = T, #scale = "row", 
                       cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                       border_color = NA, breaks = seq(-2,2,length.out = 101), silent = F, 
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                       #annotation_row = ann_row, annotation_names_row = T, annotation_colors = ann_colors, 
                       main = "Heatmap of marker genes: z-scores (skin)")
heat.sk.he <- pheatmap(he.gn, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = T, #scale = "row", 
                       cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                       border_color = NA, breaks = seq(0,3,length.out = 101), silent = F, 
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                       #annotation_row = ann_row, annotation_names_row = T, annotation_colors = ann_colors, 
                       main = "Heatmap of marker genes:\nmean log expression (skin)")

#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for everage expression of genes we have antibodies available
#-----------------------------------------------------------------------------------------------------------------#

# Genes with Abs available in the lab CD8????? IL12AorB??? IL17??? STAT5??? IL1RA???
ga <- c('AKT1', 'PYCARD', 'BCL6', 'ACTB', 'CASP8', 'CCR10', 'CCR3', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR8', 'CCR9', 
        'ITGAE', 'LAMP1', 'KIT', 'ITGAM', 'ITGAX', 'IL2RB', 'IL3RA', 'IL7R', 'CD14', 'THBD', 'CD40LG', 'FCGR3A', 
        'KLRB1', 'FASLG', 'CD19', 'CD1A', 'CD1C', 'ABCB1', 'IL2RA', 'CD27', 'NGFR', 'ICOS', 'CD28', 'CD3E', 'PECAM1', 
        'NCR2', 'CD34', 'CD4', 'CD40', 'PTPRC', 'CD46', 'ITGA1', 'ITGA2', 'CD5', 'NCAM1', 'SELL', 'CD69', 'CD8A', 
        'CD80', 'CD83', 'CD86', 'FAS', 'SELPLG', 'PTGDR2', 'CTLA4', 'CXCR3', 'MUC1', 'EOMES', 'FOXP3', 'GATA3', 
        'SLC2A1', 'CSF2', 'GNLY', 'GZMA', 'GZMB', 'GPR15', 'HLA-A', 'HLA-DRA', 'ZNF683', 'IFNG', 'IL10', 'IL12A', 
        'IL12RB2', 'IL13', 'IL17A', 'IL17F', 'IL1A', 'IL1B', 'IL1R1', 'IL1R2', 'IL1RN', 'IL2', 'IL21', 
        'IL22', 'IL26', 'IL4', 'IL4R', 'IL5', 'IL6', 'IL6R', 'CXCL8', 'IL9', 'ITGA4', 'MKI67', 'KLRG1', 'LAG3', 
        'MTOR', 'NLRP3', 'NFAT5', 'MAPK14', 'MAPK11', 'PDCD1', 'CD274', 'PRF1', 'FOXO1', 'FOXO1', 'MAPKAPK2', 'MAPK1', 
        'RPS6', 'SGK1', 'PDPN', 'RORC', 'S1PR1', 'IL1RL1', 'STAT1', 'STAT3', 'STAT4', 'STAT5A', 'STAT6', 'TBX21', 'TRBC1', 
        'TRDC', 'TGFB1', 'TLR2', 'TLR4', 'TLR9', 'TNF')

# Filter data for cluster markers
#hd.ga <- hd[match(unique(ga), row.names(hd)), -1]
he.ga <- he[match(unique(ga), row.names(he)), -1]
he.ga <- he.ga[rowSums(he.ga > 0.03) >= 1, ]
#he.ga <- he.ga[!(rowSums(he.ga) < 0.5), ]
hd.ga <- hd[match(row.names(he.ga), row.names(hd)), -1]

# Draw heatmaps
heat.sk.hz <- pheatmap(hd.ga, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = T, #scale = "row", 
                       cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                       border_color = NA, breaks = seq(-2,2,length.out = 101), silent = F, 
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                       #annotation_row = ann_row, annotation_names_row = T, annotation_colors = ann_colors, 
                       main = "Heatmap of marker genes: z-scores (skin)")
heat.sk.hl <- pheatmap(he.ga, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = T, #scale = "row", 
                       cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                       border_color = NA, breaks = seq(0,3,length.out = 101), silent = F, 
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                       #annotation_row = ann_row, annotation_names_row = T, annotation_colors = ann_colors, 
                       main = "Heatmap of marker genes:\nmean log expression (skin)")

#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for everage expression of genes analyzed by cyTOF in the atlas manuscript
#-----------------------------------------------------------------------------------------------------------------#

all <- c('CD4', 'CD8A', 'CD45RA', 'CD45', 'CD103', 'CD69', 'CXCR6', 'CD49a', 'CCR10', 'CCR4', 
         'CD14', 'CD14.1', 'CD3', 'CD161', 'CCR9', 'CXCR3', 'CD95', 'CXCR5', 'CD49d', 'CCR2', 'CCR5', 'CD19', 
         'CD27', 'CD56', 'PD-1', 'CD31', 'CD57', 'TCRgD', 'CXCR4', 'CD127', 'CD29', 'CD38', 'CD62L', 'CLA', 
         'HLA-DR', 'Cisplatin', 'CD25', 'ICOS', 'IntegrinB7', 'CD45', 'CCR6', 'CCR7', 'CX3CR1')

#-----------------------------------------------------------------------------------------------------------------#
# Heatmap for percentage of cells expressing the markers
#-----------------------------------------------------------------------------------------------------------------#

# Calculate percentage of cells in each cluster expressing each gene
hp <- plyr::ddply(ht, c("Var1", "cluster"), summarise, pct = (sum(value > 0)/length(value)))
write.csv(hp, file = paste0(dir.figs.hm, "_cluster_percent_expression.csv"), row.names = T)
#hp <- read.csv(file = paste0(dir.figs.hm, "_cluster_percent_expression.csv"), row.names = 1, stringsAsFactors = T)

# Transform data frame
hp$cluster <- as.numeric(as.character(hp$cluster))
hh <- reshape2::dcast(hp, Var1~cluster, value.var = "pct")
rownames(hh) <- hh$Var1
hh[is.na(hh)] <- 0

# Filter data for cluster markers
hh.gn <- hh[match(unique(gn), row.names(hh)), -1]

# Draw heatmaps with gene annotation
heat.sk.pc <- pheatmap(hh.gn, cluster_rows = T, show_rownames = T, show_colnames = T, cluster_cols = T, #scale = "row", 
                       cellwidth = 10, cellheight = 10, treeheight_row = 20, angle_col = 270, 
                       border_color = NA, breaks = seq(0,1,length.out = 101), silent = F, 
                       col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                       #annotation_row = ann_row, annotation_names_row = T, annotation_colors = ann_colors, 
                       main = "Heatmap for percentage of cells\nexpressing the marker genes (skin)")

#-----------------------------------------------------------------------------------------------------------------#
# Check for doublets
#-----------------------------------------------------------------------------------------------------------------#

# Identify cells with higher feature expression less stringently
dbl <- subset(sk, subset = nFeature_RNA > 1400 & nFeature_RNA < 2000)@assays[["RNA"]]@data@Dimnames[[2]]

# Plot cells
ts.dbl <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1", label = T, label.size = 6, pt.size = 1, 
        cols = NULL, cells = dbl) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("skin")
um.dbl <- DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1", label = T, label.size = 6, pt.size = 1, 
        cols = NULL, cells = dbl) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("skin")

# Identify cells with higher feature expression more strigently
dbm <- subset(sk, subset = nFeature_RNA > 1200 & nFeature_RNA < 2000)@assays[["RNA"]]@data@Dimnames[[2]]

# Plot cells
ts.dbm <- DimPlot(sk, reduction = "tsne", dims = c(1,2), group.by = "RNA_snn_res.1", label = T, label.size = 6, pt.size = 1, 
                 cols = NULL, cells = dbm) + theme_bw() + theme_custom + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle("skin")
um.dbm <- DimPlot(sk, reduction = "umap", dims = c(1,2), group.by = "RNA_snn_res.1", label = T, label.size = 6, pt.size = 1, 
                 cols = NULL, cells = dbm) + theme_bw() + theme_custom + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("skin")

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Save tSNE plots for skin from Seurat analysis of cluster markers
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
  ggsave(plot = ts.mk.sk[[i+1]], file = paste0(dir.figs.ts, "_tsne_cluster_", i, ".pdf"), height = 9, width = 10)
}

# Save UMAP plots for skin from Seurat analysis of cluster markers
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
  ggsave(plot = um.mk.sk[[i+1]], file = paste0(dir.figs.um, "_umap_cluster_", i, ".pdf"), height = 9, width = 10)
}

# Save violin plots plots for skin from Seurat analysis of cluster markers
for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
  ggsave(plot = vp.mk.sk[[i+1]], file = paste0(dir.figs.vp, "_violin_cluster_", i, ".pdf"), height = 9, width = 10)
}

# Save heatmaps from Seurat analysis
ggsave(plot = hm.mk.sk, file = paste0(dir.figs.hm, "_heatmap_cluster_markers.pdf"), height = 12, width = 15)
ggsave(plot = hm.mn.sk, file = paste0(dir.figs.hm, "_heatmap_cluster_markers_negative.pdf"), height = 7, width = 15)

# Save customized heatmaps
pdf(paste0(dir.figs.hm, "_heatmap_cluster_markers_mean_zscores.pdf"), height = 14, width = 18)
grid.arrange(grobs = list(heat.sk[[1]][[4]], heat.sk[[2]][[4]], heat.sk[[3]][[4]], heat.sk[[4]][[4]], heat.sk[[5]][[4]], 
                          heat.sk[[6]][[4]], heat.sk[[7]][[4]], heat.sk[[8]][[4]], heat.sk[[9]][[4]], heat.sk[[10]][[4]], 
                          heat.sk[[11]][[4]], heat.sk[[12]][[4]], heat.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers (skin): z-scores", gp = gpar(fontsize = 14, font = 2), vjust = 0.6))
dev.off()
pdf(paste0(dir.figs.hm, "_heatmap_cluster_markers_mean_zscores_scaled.pdf"), height = 14, width = 18)
grid.arrange(grobs = list(heats.sk[[1]][[4]], heats.sk[[2]][[4]], heats.sk[[3]][[4]], heats.sk[[4]][[4]], heats.sk[[5]][[4]], 
                          heats.sk[[6]][[4]], heats.sk[[7]][[4]], heats.sk[[8]][[4]], heats.sk[[9]][[4]], heats.sk[[10]][[4]], 
                          heats.sk[[11]][[4]], heats.sk[[12]][[4]], heats.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers (skin): row-scaled z-scores", gp = gpar(fontsize = 14, font = 2), vjust = 0.6))
dev.off()
pdf(paste0(dir.figs.hm, "_heatmap_cluster_markers_mean_log.pdf"), height = 14, width = 18)
grid.arrange(grobs = list(heatl.sk[[1]][[4]], heatl.sk[[2]][[4]], heatl.sk[[3]][[4]], heatl.sk[[4]][[4]], 
                          heatl.sk[[5]][[4]], heatl.sk[[6]][[4]], heatl.sk[[7]][[4]], heatl.sk[[8]][[4]], 
                          heatl.sk[[9]][[4]], heatl.sk[[10]][[4]], heatl.sk[[11]][[4]], heatl.sk[[12]][[4]], 
                          heatl.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers (skin): mean log expression", 
                            gp = gpar(fontsize = 14, font = 2), vjust = 0.6))
dev.off()
pdf(paste0(dir.figs.hm, "_heatmap_cluster_markers_mean_log_expressing.pdf"), height = 14, width = 18)
grid.arrange(grobs = list(heatn.sk[[1]][[4]], heatn.sk[[2]][[4]], heatn.sk[[3]][[4]], heatn.sk[[4]][[4]], 
                          heatn.sk[[5]][[4]], heatn.sk[[6]][[4]], heatn.sk[[7]][[4]], heatn.sk[[8]][[4]], 
                          heatn.sk[[9]][[4]], heatn.sk[[10]][[4]], heatn.sk[[11]][[4]], heatn.sk[[12]][[4]], 
                          heatn.sk[[13]][[4]]), nrow = 3, 
             top = textGrob("Heatmap of top 30 cluster markers (skin): mean log expression of expressing cells", 
                            gp = gpar(fontsize = 14, font = 2), vjust = 0.6))
dev.off()

# Save heatmaps from customized targets
ggsave(plot = heat.sk.hd, file = paste0(dir.figs.hm, "_heatmap_markers_zscores.pdf"), height = 8, width = 5)
ggsave(plot = heat.sk.he, file = paste0(dir.figs.hm, "_heatmap_markers_mean_log.pdf"), height = 8, width = 5)
ggsave(plot = heat.sk.hz, file = paste0(dir.figs.hm, "_heatmap_abavail_zscores.pdf"), height = 14, width = 5)
ggsave(plot = heat.sk.hl, file = paste0(dir.figs.hm, "_heatmap_abavail_mean_log.pdf"), height = 14, width = 5)
ggsave(plot = heat.sk.pc, file = paste0(dir.figs.hm, "_heatmap_markers_percentage.pdf"), height = 8, width = 5)

# Save plot with putiative doublets (or cells with highest number of genes detected)
ggsave(plot = ts.dbl, file = paste0(dir.figs.db, "_tsne_doublets_less.pdf"), height = 5, width = 5)
ggsave(plot = um.dbl, file = paste0(dir.figs.db, "_umap_doublets_less.pdf"), height = 5, width = 5)
ggsave(plot = ts.dbm, file = paste0(dir.figs.db, "_tsne_doublets_more.pdf"), height = 5, width = 5)
ggsave(plot = um.dbm, file = paste0(dir.figs.db, "_umap_doublets_more.pdf"), height = 5, width = 5)

##-----------------------------------------------------------------------------------------------------------------#
## Exclusive marker genes (the top ones also appear within the top 30 of the general cluster markers)
##-----------------------------------------------------------------------------------------------------------------#
#
##
## Find marker genes for each cluster, which are exclusively expressed in the cluster
#em <- NULL
#for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
#  for(j in as.numeric(levels(sk@active.ident))[as.numeric(levels(sk@active.ident)) != i]) {
#emi <- FindMarkers(sk, ident.1 = i, ident.2 = j, min.pct = 0.1, logfc.threshold = 0.25, only.pos = T, 
#                  test.use = "wilcox", min.diff.pct = -Inf, random.seed = 1)
#emi[, c("genes", "cluster", "ref")] <- list(row.names(emi), i, j)
#em <- rbind(em, emi)
#}}
#write.csv(em, file = paste0(dir.figs.sk, "_cluster_exclusive_markers.csv"), row.names = T)
#ems <- em[order(em$avg_logFC, decreasing = T),]
#ems <- ems[ems$pct.2 < 0.2, ]
#ems <- ems[ems$genes %in% ems[count(ems, genes)[count(ems, genes)$n == 12, ]$genes, ]$genes, ]
#em.list <- unique(ems$genes)
#
# Visualization of marker gene expression (UMAP)
#ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
#  scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
#  scale_x_continuous(breaks = seq(-4,4,4)) + scale_y_continuous(breaks = seq(-4,4,4)) + facet_wrap(~variable, ncol = 6) + 
#  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("Expression of top 30 marker genes from cluster "))
#
#um.em.sk <- sapply(as.numeric(levels(sk@active.ident))+1, function(x) NULL)
#for(i in 0:max(as.numeric(as.character(sk@active.ident)))) {
#  df.sk <- data.frame(x = sk@reductions$umap@cell.embeddings[,1],
#                      y = sk@reductions$umap@cell.embeddings[,2],
#                      t(as.matrix(GetAssayData(sk, slot = "data")[em.list[1:36],])))
#  df.sk <- reshape2::melt(df.sk, id.vars = c('x','y'))
#  um.em.sk[[i+1]] <- ggplot(df.sk, aes(x = x, y = y, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
#    scale_color_gradient(low = "gray", high = "blue", na.value = "blue", limits = c(0,3)) + 
#    scale_x_continuous(breaks = seq(-4,4,4)) + scale_y_continuous(breaks = seq(-4,4,4)) + facet_wrap(~variable, ncol = 6) + 
#    theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0("Expression of top 30 marker genes from cluster ", i))
#}
#um.em.sk[[10]]
#
#FeaturePlot(sk, features = em.list[1:5], dims = c(1,2), cols = c("lightgrey", "blue"), pt.size = NULL, 
#            min.cutoff = NA, max.cutoff = NA, reduction = "umap", blend = FALSE, blend.threshold = 0.5, 
#            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = F)
#VlnPlot(sk, features = em.list, cols = NULL, pt.size = 0, adjust = 1,  y.max = NULL, ncol = NULL, slot = "data")

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
#date     2019-01-31                  
#
#─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#package      * version    date       lib source                           
#assertthat     0.2.0      2017-04-11 [1] CRAN (R 3.5.0)                   
#backports      1.1.2      2017-12-13 [1] CRAN (R 3.5.0)                   
#bibtex         0.4.2      2017-06-30 [1] CRAN (R 3.5.0)                   
#bindr          0.1.1      2018-03-13 [1] CRAN (R 3.5.0)                   
#bindrcpp       0.2.2      2018-03-29 [1] CRAN (R 3.5.0)                   
#bitops         1.0-6      2013-08-17 [1] CRAN (R 3.5.0)                   
#callr          3.1.0      2018-12-10 [1] CRAN (R 3.5.0)                   
#caTools        1.17.1.1   2018-07-20 [1] CRAN (R 3.5.0)                   
#cli            1.0.1      2018-09-25 [1] CRAN (R 3.5.0)                   
#cluster        2.0.7-1    2018-04-13 [1] CRAN (R 3.5.1)                   
#codetools      0.2-15     2016-10-05 [1] CRAN (R 3.5.1)                   
#colorspace     1.3-2      2016-12-14 [1] CRAN (R 3.5.0)                   
#cowplot        0.9.3      2018-07-15 [1] CRAN (R 3.5.0)                   
#crayon         1.3.4      2017-09-16 [1] CRAN (R 3.5.0)                   
#data.table     1.11.8     2018-09-30 [1] CRAN (R 3.5.0)                   
#desc           1.2.0      2018-05-01 [1] CRAN (R 3.5.0)                   
#devtools       2.0.1      2018-10-26 [1] CRAN (R 3.5.1)                   
#digest         0.6.18     2018-10-10 [1] CRAN (R 3.5.0)                   
#dplyr        * 0.7.8      2018-11-10 [1] CRAN (R 3.5.0)                   
#fitdistrplus   1.0-11     2018-09-10 [1] CRAN (R 3.5.0)                   
#fs             1.2.6      2018-08-23 [1] CRAN (R 3.5.0)                   
#future         1.10.0     2018-10-17 [1] CRAN (R 3.5.0)                   
#future.apply   1.0.1      2018-08-26 [1] CRAN (R 3.5.0)                   
#gbRd           0.4-11     2012-10-01 [1] CRAN (R 3.5.0)                   
#gdata          2.18.0     2017-06-06 [1] CRAN (R 3.5.0)                   
#ggplot2      * 3.1.0      2018-10-25 [1] CRAN (R 3.5.0)                   
#ggrepel        0.8.0      2018-05-09 [1] CRAN (R 3.5.0)                   
#ggridges       0.5.1      2018-09-27 [1] CRAN (R 3.5.0)                   
#globals        0.12.4     2018-10-11 [1] CRAN (R 3.5.0)                   
#glue           1.3.0      2018-07-17 [1] CRAN (R 3.5.0)                   
#gplots         3.0.1      2016-03-30 [1] CRAN (R 3.5.0)                   
#gridExtra    * 2.3        2017-09-09 [1] CRAN (R 3.5.0)                   
#gtable         0.2.0      2016-02-26 [1] CRAN (R 3.5.0)                   
#gtools         3.8.1      2018-06-26 [1] CRAN (R 3.5.0)                   
#htmltools      0.3.6      2017-04-28 [1] CRAN (R 3.5.0)                   
#htmlwidgets    1.3        2018-09-30 [1] CRAN (R 3.5.0)                   
#httr           1.4.0      2018-12-11 [1] CRAN (R 3.5.1)                   
#ica            1.0-2      2018-05-24 [1] CRAN (R 3.5.0)                   
#igraph         1.2.2      2018-07-27 [1] CRAN (R 3.5.0)                   
#irlba          2.3.2      2018-01-11 [1] CRAN (R 3.5.0)                   
#jsonlite       1.6        2018-12-07 [1] CRAN (R 3.5.0)                   
#KernSmooth     2.23-15    2015-06-29 [1] CRAN (R 3.5.1)                   
#lattice        0.20-38    2018-11-04 [1] CRAN (R 3.5.0)                   
#lazyeval       0.2.1      2017-10-29 [1] CRAN (R 3.5.0)                   
#listenv        0.7.0      2018-01-21 [1] CRAN (R 3.5.0)                   
#lmtest         0.9-36     2018-04-04 [1] CRAN (R 3.5.0)                   
#lsei           1.2-0      2017-10-23 [1] CRAN (R 3.5.0)                   
#magrittr       1.5        2014-11-22 [1] CRAN (R 3.5.0)                   
#MASS           7.3-51.1   2018-11-01 [1] CRAN (R 3.5.1)                   
#Matrix         1.2-15     2018-11-01 [1] CRAN (R 3.5.1)                   
#memoise        1.1.0      2017-04-21 [1] CRAN (R 3.5.0)                   
#metap          1.0        2018-07-25 [1] CRAN (R 3.5.0)                   
#munsell        0.5.0      2018-06-12 [1] CRAN (R 3.5.0)                   
#npsurv         0.4-0      2017-10-14 [1] CRAN (R 3.5.0)                   
#openxlsx     * 4.1.0      2018-05-26 [1] CRAN (R 3.5.0)                   
#pbapply        1.3-4      2018-01-10 [1] CRAN (R 3.5.0)                   
#pheatmap     * 1.0.12     2019-01-04 [1] CRAN (R 3.5.2)                   
#pillar         1.3.0      2018-07-14 [1] CRAN (R 3.5.0)                   
#pkgbuild       1.0.2      2018-10-16 [1] CRAN (R 3.5.0)                   
#pkgconfig      2.0.2      2018-08-16 [1] CRAN (R 3.5.0)                   
#pkgload        1.0.2      2018-10-29 [1] CRAN (R 3.5.0)                   
#plotly         4.8.0      2018-07-20 [1] CRAN (R 3.5.1)                   
#plyr           1.8.4      2016-06-08 [1] CRAN (R 3.5.0)                   
#png            0.1-7      2013-12-03 [1] CRAN (R 3.5.0)                   
#prettyunits    1.0.2      2015-07-13 [1] CRAN (R 3.5.0)                   
#processx       3.2.1      2018-12-05 [1] CRAN (R 3.5.1)                   
#ps             1.2.1      2018-11-06 [1] CRAN (R 3.5.1)                   
#purrr          0.2.5      2018-05-29 [1] CRAN (R 3.5.0)                   
#R.methodsS3    1.7.1      2016-02-16 [1] CRAN (R 3.5.0)                   
#R.oo           1.22.0     2018-04-22 [1] CRAN (R 3.5.0)                   
#R.utils        2.7.0      2018-08-27 [1] CRAN (R 3.5.0)                   
#R6             2.3.0      2018-10-04 [1] CRAN (R 3.5.0)                   
#RANN           2.6        2018-07-16 [1] CRAN (R 3.5.0)                   
#RColorBrewer * 1.1-2      2014-12-07 [1] CRAN (R 3.5.0)                   
#Rcpp           1.0.0      2018-11-07 [1] CRAN (R 3.5.0)                   
#Rdpack         0.10-1     2018-10-04 [1] CRAN (R 3.5.0)                   
#remotes        2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
#reticulate     1.10       2018-08-05 [1] CRAN (R 3.5.0)                   
#rlang          0.3.0.1    2018-10-25 [1] CRAN (R 3.5.0)                   
#ROCR           1.0-7      2015-03-26 [1] CRAN (R 3.5.0)                   
#rprojroot      1.3-2      2018-01-03 [1] CRAN (R 3.5.0)                   
#rstudioapi     0.8        2018-10-02 [1] CRAN (R 3.5.0)                   
#rsvd           1.0.0      2018-11-06 [1] CRAN (R 3.5.0)                   
#Rtsne          0.15       2018-11-10 [1] CRAN (R 3.5.0)                   
#scales         1.0.0      2018-08-09 [1] CRAN (R 3.5.0)                   
#SDMTools       1.1-221    2014-08-05 [1] CRAN (R 3.5.0)                   
#sessioninfo    1.1.1      2018-11-05 [1] CRAN (R 3.5.1)                   
#Seurat       * 3.0.0.9000 2018-12-12 [1] Github (satijalab/seurat@aa9ffda)
#stringi        1.2.4      2018-07-20 [1] CRAN (R 3.5.0)                   
#stringr        1.3.1      2018-05-10 [1] CRAN (R 3.5.0)                   
#survival       2.43-3     2018-11-26 [1] CRAN (R 3.5.0)                   
#tibble         1.4.2      2018-01-22 [1] CRAN (R 3.5.0)                   
#tidyr          0.8.2      2018-10-28 [1] CRAN (R 3.5.0)                   
#tidyselect     0.2.5      2018-10-11 [1] CRAN (R 3.5.0)                   
#tsne           0.1-3      2016-07-15 [1] CRAN (R 3.5.0)                   
#usethis        1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
#viridisLite    0.3.0      2018-02-01 [1] CRAN (R 3.5.0)                   
#withr          2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
#yaml           2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
#zip            1.0.0      2017-04-25 [1] CRAN (R 3.5.0)                   
#zoo            1.8-4      2018-09-19 [1] CRAN (R 3.5.0)                   
#
#[1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library