#=================================================================================================================#
# Data description: (EX0033) scRNAseq CD8 Tm high salt and low salt 
# Date: 02.02.2023
# Author: Mahima Arunkumar
# Rscript purpose:
# - Preprocessing
# - Downstream analysis
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(DoubletFinder)
library(ggplot2)
library(RColorBrewer)
library(ArchR)
library(leiden)
library(EnhancedVolcano)
library( "DESeq2" )
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(cowplot)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(presto)
library(ComplexHeatmap)
library(enrichplot)

################################# Preprocessing ################################
# CD8 T cells
# 0170: CD8 Tm (CD8+CD45RA-) low salt
# 0171: CD8 Tm + NaCl (CD8+CD45RA-) high salt
# 0172: CD8 Tn (CD8+CD45RA+)
# 0173: CD8 Tn + NaCl (CD8+CD45RA+)
# 0174: CD8 (CD8+)

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



NVAR = 8000   # number of variable features to be computed
var_cutoff = 0.8 # to choose number of PCA components
# Load UMI counts from 10x data

seurat_0170 <- Read10X(data.dir = "raw_data/0170/")
#seurat_0170
seurat_0171 <- Read10X(data.dir = "raw_data/0171/")
#seurat_0171

seurat_0170 <- CreateSeuratObject(counts = seurat_0170$`Gene Expression`, min.cells = 0, min.features = 0, 
                                  project = "low salt",assay="RNA")
seurat_0171 <- CreateSeuratObject(counts = seurat_0171$`Gene Expression`, min.cells = 0, min.features = 0, 
                                  project = "high salt",assay="RNA")

# adding conditions
seurat_0170$condition = "low salt"
seurat_0171$condition = "high salt"

combined.object <- merge(seurat_0170, c(seurat_0171) , 
                         add.cell.ids = c("low salt", "high salt"), 
                         project = "salt")

#Find mtgenes
#combined.object[["percent.mt"]] <- PercentageFeatureSet(combined.object, pattern = "^MT-")
#VlnPlot(combined.object, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

#Filter dataset
combined.object <- subset(combined.object, 
                          nFeature_RNA < 9000 &
                            nCount_RNA < 80000)

#ncol(combined.object) ##3731
#combined.object

## normalize and scale the data
combined.object <- NormalizeData(combined.object, normalization.method = "LogNormalize", scale.factor = 1e4)
combined.object <- FindVariableFeatures(combined.object, selection.method = "vst", nfeatures = NVAR)
combined.object <- ScaleData(combined.object,vars.to.regress=c("nCount_RNA", "nFeature_RNA"))
combined.object <- RunPCA(combined.object)
#p1<-DimPlot(combined.object, reduction = "pca", group.by = "condition",pt.size=0.1) + theme_custom + ggtitle("PCA")
#p1
#ggsave("./figures/PCA.pdf",p1, units = "cm", width = 15, height = 15)

# choosing dimension of the PCA
pca = combined.object[["pca"]]
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
npca <- which(cumsum(varExplained)>var_cutoff)[1]


combined.object <- RunUMAP(combined.object, reduction = "pca", dims = 1:npca,
                           seed.use = 1, reduction.name = "umap",umap.method="uwot")
p2<-DimPlot(combined.object, reduction = "umap", group.by = "condition",pt.size=0.1) + theme_custom + ggtitle("UMAP")
p2

combined.object <- RunTSNE(combined.object, reduction = "pca", dims = 1:npca,
                           seed.use = 1, reduction.name = "tsne", tsne.method = "Rtsne")
p0t<-DimPlot(combined.object, reduction = "tsne", group.by = "condition",pt.size=0.1) + theme_custom + ggtitle("tSNE")
p0t

#ggsave("../scripts/figures/25-01-23/UMAP.pdf",p2)

combined.object <- FindNeighbors(combined.object, reduction = "pca", k.param = 20, dims = 1:npca, do.plot = F)
combined.object<- FindClusters(combined.object, algorithm = 4, random.seed = 0)
col.ls <- ArchRPalettes[1]
col.ls <- as.vector(col.ls$stallion)
p3 <- DimPlot(combined.object, reduction = "umap", group.by = "RNA_snn_res.0.8",pt.size=0.1,label=TRUE,cols=col.ls) + theme_custom + ggtitle("leiden clusters")
ggsave("../scripts/figures/25-01-23/Leiden.pdf",p3)

# identify the dead-cells
mito_genes <- rownames(combined.object)[grep(pattern = "^MT-",rownames(combined.object))]

p4<-DotPlot(object = combined.object, features = mito_genes,col.min = -1,col.max = 1) + 
  theme_custom + theme(axis.text.x = element_text(size = 8.8, color = "black",angle = 90)) + 
  scale_x_discrete(name ="Mitochondrial genes")+scale_y_discrete(name ="Leiden clusters")

#ggsave("../scripts/figures/25-01-23/Dotplot.pdf",p4)

## ---- DETECTION OF DOUBLETS USING DOUBLET FINDER ----

# doublet finder
sweep.res <- paramSweep_v3(combined.object, PCs = 1:npca, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2,xlab="pK",ylab="BCmetric") + theme_custom
#Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
#Optimal pK is the max of the bomodality coefficent (BCmvn) distribution

bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

# optimal.pk <- 0.04
## Homotypic doublet proportion estimate

annotations <- combined.object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal.pk * nrow(combined.object@meta.data)) 
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
combined.object <- doubletFinder_v3(seu = combined.object, 
                                 PCs = 1:npca, 
                                 pK = optimal.pk,
                                 nExp = nExp.poi.adj)
DF.name = colnames(combined.object@meta.data)[grepl("DF.classification", colnames(combined.object@meta.data))]

p8 <- DimPlot(combined.object, reduction = "umap",group.by = DF.name,pt.size=0.1) + theme_custom + ggtitle("Doublets")
p8
ggsave("./figures/doubletfinder_UMAP.pdf",p8, units = "cm", width = 15, height = 15)

p9 <- DimPlot(combined.object, reduction = "tsne",group.by = DF.name,pt.size=0.1) + theme_custom + ggtitle("Doublets")
p9
ggsave("./figures/doubletfinder_tSNE.pdf",p9, units = "cm", width = 15, height = 15)


dim(combined.object)
combined.object = combined.object[, combined.object@meta.data[, DF.name] == "Singlet"]
dim(combined.object)

#Save and load object
saveRDS(combined.object,"./output_files/combined_object_no_doublets")
combined.object <- readRDS("./output_files/combined_object_no_doublets")

## Clustering ##
combined.object <- FindNeighbors(combined.object, reduction = "pca", k.param = 20, do.plot = F)
combined.object<- FindClusters(combined.object, algorithm = 4, random.seed = 0, resolution = 0.8)
#col.ls <- ArchRPalettes[1]
#col.ls <- as.vector(col.ls$stallion)

combined.object.final <- RunUMAP(combined.object, reduction = "pca", dims = 1:24,
                                        seed.use = 1, 
                                        reduction.name = "umap",umap.method="uwot")

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6)) #concat two color palettes for more than 12 colors
p10<-DimPlot(combined.object.final, reduction = "umap", 
            group.by = "orig.ident",pt.size=0.1) + theme_custom + ggtitle("UMAP")
p10
ggsave("./figures/UMAP.pdf",p10, units = "cm", width = 15, height = 15)

p11<-DimPlot(combined.object.final, reduction = "tsne", 
             group.by = "orig.ident",pt.size=0.1) + theme_custom + ggtitle("tSNE")
p11
ggsave("./figures/tSNE.pdf",p11, units = "cm", width = 15, height = 15)

#Plot leiden clusters
p12 <- DimPlot(combined.object.final, reduction = "umap", pt.size=0.1, label = TRUE ) + theme_custom + ggtitle("UMAP")
p12
ggsave("./figures/UMAP_leiden.pdf",p12, units = "cm", width = 15, height = 15)

p13 <- DimPlot(combined.object.final, reduction = "tsne", pt.size=0.1, cols = mycolors) + theme_custom + ggtitle("tSNE")
p13
ggsave("./figures/tSNE_leiden.pdf",p13, units = "cm", width = 15, height = 15)

#Save and load object
saveRDS(combined.object.final,"./output_files/combined_object_after_clustering")
combined.object.final <- readRDS("./output_files/MA_combined_object_after_clustering")

#################################################################################

# Stacked bar plot (high salt/low salt percentage per cluster)
cell_clusters = dplyr::select(combined.object.final@meta.data,seurat_clusters,orig.ident)
df_rel_freq=table(cell_clusters)/as.vector(table(combined.object.final$orig.ident))*100
par(mfrow=c(1, 1), mar=c(5, 4, 5, 6))
cell_numbers = as.vector(table(combined.object.final$seurat_clusters))
no_of_clusters <- length(unique(combined.object.final$seurat_clusters))
ticks <- sprintf("%s",seq(1:no_of_clusters))
a1<-as.data.frame(df_rel_freq)
a1_low <- a1[a1$orig.ident=="low salt",]
a1_high <- a1[a1$orig.ident=="high salt",]
df= cbind(a1_high$Freq,a1_low$Freq)

#pdf(file = "../scripts/figures/25-01-23/stacked_bar_wo_xy_labels.pdf",bg="transparent",width=8,height = 4)
p3<-barplot(t(df_rel_freq),legend = c("High Salt","Low Salt"),col=c("brown1","#30D5C8"),ylab="Percentage of cells (split by high and low salt) per cluster", xlab="Leiden clusters",ylim=c(0,40),args.legend = list(bty = "n",x = "topright"),bg = 'transparent', space= 0.5, width = c(0.3,0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)) 
text(p3,apply(df, 1, sum)+1.5 , labels = cell_numbers)
ggsave("./figures/stacked_barplot.pdf",text(p3,apply(df, 1, sum)+1.5 , labels = cell_numbers), units = "cm", width = 15, height = 15)
#p3
dev.off()

#Heatmap per leiden cluster #change idents for heatmap per cluster (only.pos = TRUE)
combined_object_final.markers <- FindAllMarkers(combined.object.final, only.pos = TRUE, min.pct = 0.25, 
                                                test.use="wilcox",logfc.threshold = 0.25) #default is wilcox anyway

combined_object_final.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

p14 <- DoHeatmap(combined.object.final, features = top10$gene, group.colors = mycolors)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p14
ggsave("./figures/Heatmap_per_cluster_only_positive.pdf",p14, units = "cm", width = 60, height = 25)

#Heatmap per leiden cluster #change idents for heatmap per cluster (only.pos = FALSE)
combined_object_final_mixed.markers <- FindAllMarkers(combined.object.final, only.pos = FALSE, min.pct = 0.25, 
                                                test.use="wilcox",logfc.threshold = 0.25) #default is wilcox anyway

combined_object_final_mixed.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

p14_a <- DoHeatmap(combined.object.final, features = top10$gene, group.colors = mycolors)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p14_a
ggsave("./figures/Heatmap_per_cluster_all.pdf",p14_a, units = "cm", width = 15, height = 25)


#Heatmap without Seurat dropping genes 
markers<- presto::wilcoxauc(combined.object.final, 'seurat_clusters', assay = 'data')
markers<- top_markers(markers, n = 10)
#reorder columns
markers <- markers[, c("rank", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")]
markers

all_markers <- markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>% .[!is.na(.)]

mat<- combined.object.final[["RNA"]]@data[all_markers, ] %>% as.matrix()

## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- combined.object.final@meta.data$seurat_clusters

# what's the value range in the matrix
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

cheatmap <- Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 9),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 6),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = FALSE,
        raster_quality = 4)

cheatmap

ggsave("./figures/Heatmap_per_leiden_all_top_10_genes_per_cluster.pdf",cheatmap, units = "cm", width = 15, height = 20)


############# Differential Gene Expression analysis usnig FindMarkers ###########

#Find all markers distinguishing high salt from low salt
Idents(object = combined.object.final ) <- "orig.ident"
markers_high_vs_low <- FindMarkers(combined.object.final, ident.1 = "high salt", ident.2 = "low salt")

sum(markers_high_vs_low$avg_log2FC < 0) # downregulated genes (292)
sum(markers_high_vs_low$avg_log2FC > 0) # upregulated genes (653)

#Write DEG table to file
write.table(markers_high_vs_low, file="./output_files/DEG_high_salt_vs_low_salt.tsv", quote=FALSE, sep='\t', col.names = NA)

#Write in Excel files (one for upreg. and one for downreg.)
markers_high_vs_low <- tibble::rownames_to_column(markers_high_vs_low, "gene")
markers_high_vs_low <- markers_high_vs_low %>% rowwise() %>% mutate(differential_expression = if_else(avg_log2FC > 0,'upregulated','downregulated'))

upregulated_all_DEGs <- markers_high_vs_low[which(markers_high_vs_low$differential_expression == "upregulated"),] 
downregulated_all_DEGs <- markers_high_vs_low[which(markers_high_vs_low$differential_expression == "downregulated"),]

write.xlsx(upregulated_all_DEGs, './output_files/Upregulated_all_DEGs.xlsx')
write.xlsx(downregulated_all_DEGs, './output_files/Downregulated_all_DEGs.xlsx')

symbols <- mapIds(org.Hs.eg.db, keys = rownames(markers_high_vs_low), keytype = "SYMBOL", column="ENSEMBL")

#Heatmap per group (high salt/low salt)

#Filter markers by p-value and log-fold change
up_markers <- markers_high_vs_low[markers_high_vs_low$p_val_adj < 0.05 & markers_high_vs_low$avg_log2FC > 0, ] %>% arrange(desc(avg_log2FC)) %>% slice(1:20) #647
down_markers <- markers_high_vs_low[markers_high_vs_low$p_val_adj < 0.05 & markers_high_vs_low$avg_log2FC < 0, ] %>% arrange(avg_log2FC) %>% slice(1:20) #290

# Get gene expression values for the selected markers
#heatmap_data_up <- as.data.frame(combined.object.final@assays$RNA@data[rownames(up_markers), ])
# Plot heatmap
#p16_a <- pheatmap::pheatmap(heatmap_data_up, scale = "row", color= c("blue", "white", "red"), show_colnames = FALSE)
#p16_a
#up_markers %>% top_n(n = 20, wt = avg_log2FC) -> up_top20

p16_a <- DoHeatmap(combined.object.final, group.by = "orig.ident", features = rownames(up_markers))  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_a
ggsave("./figures/Heatmap_per_salt_identity_upregulated.pdf",p16_a, units = "cm", width = 10, height = 15)

p16_b <- DoHeatmap(combined.object.final, group.by = "orig.ident", slot= "counts", features = rownames(down_markers))  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_b
ggsave("./figures/Heatmap_per_salt_identity_downregulated.pdf",p16_b, units = "cm", width = 10, height = 15)


all_markers <- c(rownames(up_markers), rownames(down_markers))
p16_c <- DoHeatmap(combined.object.final, group.by = "orig.ident", slot= "counts", features = all_markers)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_c
ggsave("./figures/Heatmap_per_salt_identity_up_and_downregulated.pdf",p16_c, units = "cm", width = 10, height = 15)


#Volcanoplot
p15 <- EnhancedVolcano(markers_high_vs_low,
                      lab = rownames(markers_high_vs_low),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 10e-32,
                      labSize = 2,
                      pointSize = c(ifelse(markers_high_vs_low$avg_log2FC>2, 2, 0.5)),
                      FCcutoff = 0.5,
                      title = 'High salt versus Low salt') + theme_custom
p15
ggsave("./figures/Volcanoplot_high_vs_low_salt.pdf",p15, units = "cm", width = 20, height = 15)


###################### Enrichment analysis high vs. low salt ####################
DEG_genes = subset(markers_high_vs_low, p_val_adj < 0.05) #(937)

# we want the log2 fold change 
original_gene_list <- DEG_genes$avg_log2FC

# name the vector
names(original_gene_list) <- rownames(DEG_genes)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
write.csv(names(gene_list),file="./output_files/genelist_sorted_for_enrichment.csv")

#gene_list <- read.csv("./output_files/genelist_sorted_for_enrichment.csv")
#seperate list for upregulated and downregulated

#All genes from the list
gene_list_up = gene_list[gene_list > 0]
gene_list_dn = gene_list[gene_list < 0]

#upregulated significant genes
enrich_GO_up_res <- enrichGO(gene = names(gene_list_up), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                      ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#downregulated significant genes
enrich_GO_dn_res <- enrichGO(gene = names(gene_list_dn), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                      ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff= 0.05, readable = T)

#Write enrichment tables to file
write.table(enrich_GO_up_res, file="./output_files/enrich_GO_up_res.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(enrich_GO_dn_res, file="./output_files/enrich_GO_dn_res.tsv", quote=FALSE, sep='\t', col.names = NA)

#Read in enrichment tables
enrich_GO_up_res <- read.delim("./output_files/enrich_GO_up_res.tsv", sep = "\t")
enrich_GO_dn_res <- read.delim("./output_files/enrich_GO_dn_res.tsv", sep = "\t")
  

#Enrichment dotplot for only upregulated genes
selected_pathways <- c("regulation of T cell activation", "regulation of cell-cell adhesion" , "mononuclear cell differentiation" , "extrinsic apoptotic signaling pathway", "regulation of antigen receptor-mediated signaling pathway",
                       "positive regulation of nitric oxide metabolic process" , "cellular response to interleukin-1", "positive regulation of cellular amide metabolic process" , "positive regulation of leukocyte activation",
                       "negative regulation of phosphate metabolic process", "regulation of carbohydrate metabolic process", "regulation of vitamin metabolic process", "pyrimidine ribonucleoside metabolic process", "nitric oxide metabolic process", "regulation of generation of precursor metabolites and energy")

dp_up <- dotplot(enrich_GO_up_res, showCategory=selected_pathways) + ggtitle("Overrepresentation: GO term \n (Biological process) for upregulated genes")
ggsave("./figures/Enrichment_for_upreg_signif_genes_BP_GO.pdf",dp_up, units = "cm", width = 15, height = 20)

#Enrichment dotplot for only downregulated genes
dp_dn <- dotplot(enrich_GO_dn_res, showCategory=15) + ggtitle("Overrepresentation: GO term \n (Biological process) for downregulated genes")
ggsave("./figures/Enrichment_for_downreg_signif_genes_BP_GO.pdf",dp_dn, units = "cm", width = 15, height = 20)


#Metabolism specific GO terms with groupGO() for upregulated DEGs
level_4_GO_metabolism <- groupGO(
  names(gene_list_up),
  'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  level = 4,
  readable = FALSE
)
write.table(level_4_GO_metabolism, file="./output_files/high_vs_low_salt_level_4_metabolism_GOterms.tsv", quote=FALSE, sep='\t', col.names = NA)


#Enrichment dotplot 
categories <- c("protein metabolic process", "generation of precursor metabolites and energy", "1regulation of macromolecule metabolic process", 
                "positive regulation of nitrogen compound metabolic process",
                "positive regulation of catabolic process", "cellular biosynthetic process", 
                "positive regulation of macromolecule metabolic process", "positive regulation of cellular metabolic process", "positive regulation of metabolic process", "regulation of cellular metabolic process", "regulation of primary metabolic process",
                "cellular aromatic compound metabolic process", "phosphorus metabolic process", "heterocycle metabolic process", "positive regulation of lipid metabolic process", "organic cyclic compound metabolic process")

enrich_upregulated_level_4 <- barplot(level_4_GO_metabolism, showCategory=categories, drop = FALSE) + ggtitle("GO terms related to metabolism \n for all significant upregulated DEGs (high salt vs. low salt)")
enrich_upregulated_level_4

ggsave("./figures/Enrichment__upregulated_level_4.pdf",enrich_upregulated_level_4, units = "cm", width = 19, height = 22)

#Subset to various GO term levels
#level 3
enrich_GO_up_res_filtered <- gofilter(enrich_GO_up_res, level = 3)
write.table(enrich_GO_up_res_filtered, file="./output_files/MA_high_vs_low_salt_level_3_all_GOterms.tsv", quote=FALSE, sep='\t', col.names = NA)


#level 4
enrich_GO_up_res_filtered <- gofilter(enrich_GO_up_res, level = 4)
write.table(enrich_GO_up_res_filtered, file="./output_files/MA_high_vs_low_salt_level_4_all_GOterms.tsv", quote=FALSE, sep='\t', col.names = NA)

#level 5
enrich_GO_up_res_filtered <- gofilter(enrich_GO_up_res, level = 5)
write.table(enrich_GO_up_res_filtered, file="./output_files/MA_high_vs_low_salt_level_5_all_GOterms.tsv", quote=FALSE, sep='\t', col.names = NA)


# Order terms
#upregulated
up_order_levels <- enrich_GO_up_res %>%
  arrange(p.adjust) %>%
  pull(Description)

enrich_GO_up_res$Description <- factor(enrich_GO_up_res$Description, 
                                     levels = rev(up_order_levels))

#Create fresh labels for those which are too long and plot with added features
bp_up_labels <- str_to_sentence(enrich_GO_up_res[1:20,]$Description)
#bp_up_labels[c(15, 16)] <- c("Nucleoside biosynthetic process", "Purine compound biosynthetic process" ) 

tbr2_up_plot <- ggplot(enrich_GO_up_res[1:20,]) +
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


#downregulated
down_order_levels <- enrich_GO_dn_res %>%
  arrange(p.adjust) %>%
  pull(Description)

enrich_GO_dn_res$Description <- factor(enrich_GO_dn_res$Description, 
                                       levels = rev(down_order_levels))


#Create fresh labels for those which are too long and plot with added features
bp_down_labels <- str_to_sentence(enrich_GO_dn_res[1:20,]$Description)
bp_down_labels[c(15, 18, 19, 20)] <- c("RNA splicing with bulged adenosine", "RNA splicing (transester. reactions)", "Pos. reg. of chr. organization", "Pos. reg. of DNA biosynthetic process") 

tbr2_down_plot <- ggplot(enrich_GO_dn_res[1:20, ]) +
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



#Combined plot for single enrichment plot showing up and downregulated pathways
# Tbr2 grid
tbr2 <- plot_grid(tbr2_down_plot,
                  tbr2_up_plot)
tbr2


# Adding annotations

#enrich_GO_dn_res[1:20, ]$Count
#enrich_GO_up_res[1:20, ]$Count

annotated_tbr2 <- tbr2 +
  draw_label(label = "Number of \n genes", 
             x = 0.5,
             y = 0.96,
             size = 8) +
  draw_label(label = "41", 
             x = 0.49,
             y = 0.925,
             size = 7) +
  draw_label(label = "45", 
             x = 0.49,
             y = 0.89,
             size = 7) +
  draw_label(label = "29", 
             x = 0.49,
             y = 0.84,
             size = 7) +
  draw_label(label = "24", 
             x = 0.49,
             y = 0.795,
             size = 7) +
  draw_label(label = "24", 
             x = 0.49,
             y = 0.75,
             size = 7) +
  draw_label(label = "16", 
             x = 0.49,
             y = 0.701,
             size = 7) +
  draw_label(label = "23", 
             x = 0.49,
             y = 0.66,
             size = 7)  +
  draw_label(label = "22", 
             x = 0.49,
             y = 0.615,
             size = 7) +
  draw_label(label = "24", 
             x = 0.49,
             y = 0.57,
             size = 7) +
  draw_label(label = "24", 
             x = 0.49,
             y = 0.525,
             size = 7) +
  draw_label(label = "13", 
             x = 0.49,
             y = 0.48,
             size = 7) +
  draw_label(label = "28", 
             x = 0.49,
             y = 0.435,
             size = 7) +
  draw_label(label = "16", 
             x = 0.49,
             y = 0.385,
             size = 7) +
  draw_label(label = "19", 
             x = 0.49,
             y = 0.341,
             size = 7)  +
  draw_label(label = "23", 
             x = 0.49,
             y = 0.30,
             size = 7)  +
  draw_label(label = "23", 
             x = 0.49,
             y = 0.255,
             size = 7)  +
  draw_label(label = "12", 
             x = 0.49,
             y = 0.21,
             size = 7)  +
  draw_label(label = "23", 
             x = 0.49,
             y = 0.165,
             size = 7)  +
  draw_label(label = "13", 
             x = 0.49,
             y = 0.113,
             size = 7)  +
  draw_label(label = "12", 
             x = 0.49,
             y = 0.07,
             size = 7)  +
  draw_label(label = "44", 
             x = 0.51,
             y = 0.925,
             size = 7) +
  draw_label(label = "48", 
             x = 0.51,
             y = 0.89,
             size = 7) +
  draw_label(label = "43", 
             x = 0.51,
             y = 0.84,
             size = 7) +
  draw_label(label = "46", 
             x = 0.51,
             y = 0.795,
             size = 7) +
  draw_label(label = "37", 
             x = 0.51,
             y = 0.75,
             size = 7) +
  draw_label(label = "40", 
             x = 0.51,
             y = 0.701,
             size = 7) +
  draw_label(label = "45", 
             x = 0.51,
             y = 0.66,
             size = 7)  +
  draw_label(label = "47", 
             x = 0.51,
             y = 0.615,
             size = 7) +
  draw_label(label = "39", 
             x = 0.51,
             y = 0.57,
             size = 7) +
  draw_label(label = "46", 
             x = 0.51,
             y = 0.525,
             size = 7) +
  draw_label(label = "40", 
             x = 0.51,
             y = 0.48,
             size = 7) +
  draw_label(label = "25", 
             x = 0.51,
             y = 0.435,
             size = 7) +
  draw_label(label = "25", 
             x = 0.51,
             y = 0.385,
             size = 7) +
  draw_label(label = "38", 
             x = 0.51,
             y = 0.341,
             size = 7)  +
  draw_label(label = "40", 
             x = 0.51,
             y = 0.30,
             size = 7)  +
  draw_label(label = "37", 
             x = 0.51,
             y = 0.255,
             size = 7)  +
  draw_label(label = "40", 
             x = 0.51,
             y = 0.21,
             size = 7)  +
  draw_label(label = "27", 
             x = 0.51,
             y = 0.165,
             size = 7)  +
  draw_label(label = "40", 
             x = 0.51,
             y = 0.113,
             size = 7)  +
  draw_label(label = "32", 
             x = 0.51,
             y = 0.07,
             size = 7)  +
  draw_label(label = "Overrepresentation: GO term (Biological process)", 
             x = 0.5,
             y = 1,
             hjust = 0.5,
             vjust = 1,
             size = 9) 

annotated_tbr2

# Saving enrichment plot as figure
ggsave(plot = annotated_tbr2,
       device = "pdf",
       filename = "./figures/enrichment_go_result_up_and_down.pdf",
       width = 15.5,
       height = 6.25,
       units = "in",
       dpi = 600)


gse <- gseGO(geneList=gene_list, 
             ont ="BP",
             nPerm = 10000,
             keyType = "SYMBOL",
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH")

pdot <- dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign,scales="free",
                                                                  space="free_y") + theme(axis.text.y = element_text(size = 4, color = "black"),
                                                                                          axis.text.x = element_text(size = 3, color = "black"),axis.title.x = element_text(size = 10, margin = margin(6,0,0,0)),
                                                                                          panel.background = element_rect(fill='transparent'),
                                                                                          plot.background = element_rect(fill='transparent', color=NA),) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5))

pdot
head(gse)
#ggsave("../scripts/figures/31-01-2023/cluster_prof_gseGO_high_vs_low_30.pdf",pdot,dpi=320,bg='transparent')
#barplot(gse)



########################### Enrichment for top 50 significant genes #############

#Top 50 genes from the list
gene_list_up_50 = gene_list[gene_list > 0][1:50]
gene_list_dn_50 = gene_list[gene_list < 0][1:50]

#upregulated significant genes
enrich_GO_up_res_50 <- enrichGO(gene = names(gene_list_up_50), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#downregulated significant genes
enrich_GO_dn_res_50 <- enrichGO(gene = names(gene_list_dn_50), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff= 0.05, readable = T)

#Write enrichment tables to file
write.table(enrich_GO_up_res_50, file="./output_files/enrich_GO_up_res_50.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(enrich_GO_dn_res_50, file="./output_files/enrich_GO_dn_res_50.tsv", quote=FALSE, sep='\t', col.names = NA)

p17 <- dotplot(enrich_GO_up_res, showCategory=20) + ggtitle("Overrepresentation: GO term \n (Biological process) for top 50 signif. upregulated genes")
ggsave("./figures/Enrichment_top_50_upregulated.pdf",p17, units = "cm", width = 20, height = 25)

p18 <- dotplot(enrich_GO_dn_res, showCategory=20) + ggtitle("Overrepresentation: GO term \n (Biological process) for top 50 signif. downregulated genes")
ggsave("./figures/Enrichment_top_50_downregulated.pdf",p18, units = "cm", width = 20, height = 25)




#Do enrichment anaylsis and check if cell death is enriched
markers_cluster_12 <- FindMarkers(object = combined.object.final, ident.1 = "12")
markers_cluster_12 <- markers_cluster_12[markers_cluster_12$p_val_adj < 0.05, ]
markers_cluster_12_fin <- markers_cluster_12[markers_cluster_12$p_val_adj < 0.05, ] %>% slice(1:100)

#Filter markers by p-value and log-fold change
up_markers_cluster_12 <- markers_cluster_12[markers_cluster_12$p_val_adj < 0.05 & markers_cluster_12$avg_log2FC > 0, ] %>% arrange(desc(avg_log2FC)) %>% slice(1:20) #647
down_markers_cluster_12 <- markers_cluster_12[markers_cluster_12$p_val_adj < 0.05 & markers_cluster_12$avg_log2FC < 0, ] %>% arrange(avg_log2FC) %>% slice(1:20) #290


#Convert DEGs to a list of Entrez IDs
#gene_list_fin <- as.list(rownames(markers_cluster_12_fin))
gene_list <- as.list(rownames(markers_cluster_12))
gene_list_up <- as.list(rownames(up_markers_cluster_12))
gene_list_down <- as.list(rownames(down_markers_cluster_12))

cluster_12_enrich_GO_res_all <- enrichGO(gene = gene_list, keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                                         ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)


#cluster_12_enrich_GO_res_upregulated <- enrichGO(gene = gene_list_up, keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
#ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#cluster_12_enrich_GO_res_downregulated <- enrichGO(gene = gene_list_down, keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
#ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#Write enrichment tables to file
write.table(cluster_12_enrich_GO_res_all, file="./output_files/cluster_12_enrich_GO_res_all.tsv", quote=FALSE, sep='\t', col.names = NA)
#write.table(cluster_12_enrich_GO_res_upregulated, file="./output_files/cluster_12_enrich_GO_res_upregulated.tsv", quote=FALSE, sep='\t', col.names = NA)
#write.table(cluster_12_enrich_GO_res_downregulated, file="./output_files/cluster_12_enrich_GO_res_downregulated.tsv", quote=FALSE, sep='\t', col.names = NA)


#Enrichment dotplot
p20 <- dotplot(cluster_12_enrich_GO_res_all, showCategory=30) + ggtitle("Overrepresentation: GO term (Biological process) for all \n significant cluster 12 marker genes")
p20
ggsave("./figures/Enrichment_cluster_12_result_all_signif_DEGs_30_GOterms.pdf",p20, units = "cm", width = 20, height = 45)




# Enrichment analysis for cluster 2,3,10,11 vs rest
Idents(object = combined.object.final ) <- "seurat_clusters"
#cluster_group <- c("2","3","10","11") #high salt clusters
#rest_clusters <- setdiff(1:length(Idents(combined.object.final)), cluster_group)

# Perform DEG analysis
cluster_markers_high_vs_low <- FindMarkers(combined.object.final, ident.1 = c("2","3","10","11"), ident.2 = c("1", "4", "5", "6", "7", "8", "9", "12", "13", "14"))
cluster_DEG_genes = subset(cluster_markers_high_vs_low, p_val_adj < 0.05) #(781)

# we want the log2 fold change 
cluster_original_gene_list <- cluster_DEG_genes$avg_log2FC

# name the vector
names(cluster_original_gene_list) <- rownames(cluster_DEG_genes)

# omit any NA values 
gene_list<-na.omit(cluster_original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
write.csv(names(gene_list),file="./output_files/genelist_sorted_for_enrichment_for_clusters.csv")


#All genes from the list
gene_list_up = gene_list[gene_list > 0]
gene_list_dn = gene_list[gene_list < 0]

sum(cluster_DEG_genes$avg_log2FC < 0) # downregulated genes (103)
sum(cluster_DEG_genes$avg_log2FC > 0) # upregulated genes (678)

#upregulated significant genes
enrich_GO_up_res <- enrichGO(gene = names(gene_list_up), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#downregulated significant genes
enrich_GO_dn_res <- enrichGO(gene = names(gene_list_dn), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff= 0.05, readable = T)

#Write enrichment tables to file
write.table(enrich_GO_up_res, file="./output_files/enrich_clusters_GO_up_res.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(enrich_GO_dn_res, file="./output_files/enrich_clusters_GO_dn_res.tsv", quote=FALSE, sep='\t', col.names = NA)

#Enrichment dotplot for only upregulated genes
p21 <- dotplot(enrich_GO_up_res, showCategory=15) + ggtitle("Overrepresentation: GO term (Biological process) \n for upregulated genes ([2,3,10,11] vs. rest)")
ggsave("./figures/Enrichment_clusters_for_upreg_signif_genes_BP_GO.pdf",p21, units = "cm", width = 17, height = 20)

#Enrichment dotplot for only downregulated genes
p22 <- dotplot(enrich_GO_dn_res, showCategory=15) + ggtitle("Overrepresentation: GO term (Biological process) \n for downregulated genes ([2,3,10,11] vs. rest)")
ggsave("./figures/Enrichment_clusters_for_downreg_signif_genes_BP_GO.pdf",p22, units = "cm", width = 17, height = 20)


#GO term level 5 -> looking specifically at cell death terms given DEG for clusters [2,3,10,11] vs. rest
cluster_level_5_GO <- groupGO(
  rownames(cluster_DEG_genes),
  'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  level = 5,
  readable = FALSE
)
write.table(cluster_level_5_GO, file="./output_files/cluster12_vs_rest_level_5_GOterms.tsv", quote=FALSE, sep='\t', col.names = NA)
#cluster_level_5_GO <- read.delim(file="./output_files/cluster12_vs_rest_level_5_GOterms.tsv")

#Enrichment dotplot 
categories <- c("T cell apoptotic process", "necrotic cell death", "autophagic cell death", "cell death in response to oxidative stress", "neuron death", "regulation of extrinsic apoptotic signaling pathway via death domain receptors", "cellular response to reactive oxygen species", 
                "extrinsic apoptotic signaling pathway in absence of ligand", "negative regulation of extrinsic apoptotic signaling pathway", "regulation of mitochondrial outer membrane permeabilization involved in apoptotic signaling pathway", "intrinsic apoptotic signaling pathway",
                "regulation of apoptotic signaling pathway", "regulation of intrinsic apoptotic signaling pathway", "negative regulation of apoptotic signaling pathway", "positive regulation of apoptotic signaling pathway", "extrinsic apoptotic signaling pathway", 
                "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress", "intrinsic apoptotic signaling pathway in response to DNA damage", "regulation of cysteine-type endopeptidase activity involved in apoptotic process", 
                "negative regulation of intrinsic apoptotic signaling pathway", "intrinsic apoptotic signaling pathway by p53 class mediator", "leukocyte apoptotic process", "apoptotic mitochondrial changes", "lymphocyte apoptotic process", "positive regulation of intrinsic apoptotic signaling pathway", 
                "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator", "positive regulation of cysteine-type endopeptidase activity involved in apoptotic process", "regulation of oxidative stress-induced cell death", "mitochondrial outer membrane permeabilization involved in programmed cell death", 
                "programmed necrotic cell death", "cell death in response to hydrogen peroxide", "negative regulation of oxidative stress-induced cell death", "autophagic cell death", "regulation of oxidative stress-induced cell death", "mitochondrial outer membrane permeabilization involved in programmed cell death", "cell death",
                "negative regulation of neuron death", "regulation of neuron apoptotic process", "response to tumor necrosis factor", "cellular response to tumor necrosis factor", "programmed necrotic cell death", "necroptotic process", "tumor necrosis factor-mediated signaling pathway", "cellular response to tumor necrosis factor",
                "negative regulation of stress-activated MAPK cascade", "negative regulation of stress-activated protein kinase signaling cascade", "negative regulation of oxidative stress-induced intrinsic apoptotic signaling pathway", "cellular response to chemical stress")

enrich_cluster_level_5 <- barplot(cluster_level_5_GO, showCategory=categories, drop = FALSE) + ggtitle("GO terms related to cell death \n for all significant DEGs ([2,3,10,11] vs. rest)")
enrich_cluster_level_5

ggsave("./figures/Enrichment_clusters_for_celldeath_GO_terms_all_signif_DEGs.pdf",enrich_cluster_level_5, units = "cm", width = 17, height = 22)


#Heatmap for DEG for clusters [2,3,10,11] vs. rest
# Check if cells belong to high salt cluster or low salt cluster
combined.object.final$salt_label <- ifelse(test = combined.object.final$seurat_clusters %in% c(2,3,10,11), yes = "high salt cluster", no = "low salt cluster")

#Filter markers by p-value and log-fold change
clusters_up_markers <- cluster_markers_high_vs_low[cluster_markers_high_vs_low$p_val_adj < 0.05 & cluster_markers_high_vs_low$avg_log2FC > 0, ] %>% arrange(desc(avg_log2FC)) %>% slice(1:20)
clusters_down_markers <- cluster_markers_high_vs_low[cluster_markers_high_vs_low$p_val_adj < 0.05 & cluster_markers_high_vs_low$avg_log2FC < 0, ] %>% arrange(avg_log2FC) %>% slice(1:20) 

# Get gene expression values for the selected markers
clusterheat_a <- DoHeatmap(combined.object.final, group.by = "salt_label", features = rownames(clusters_up_markers), size = 3)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
clusterheat_a
ggsave("./figures/Heatmap_per_salt_cluster_identity_upregulated.pdf",clusterheat_a, units = "cm", width = 15, height = 15)

clusterheat_b <- DoHeatmap(combined.object.final, group.by = "salt_label", features = rownames(clusters_down_markers), size = 3)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
clusterheat_b
ggsave("./figures/Heatmap_per_salt_cluster_identity_downregulated.pdf",clusterheat_b, units = "cm", width = 15, height = 15)

##################################################################################

################# Violinplot and Featureplots for genesets of interest ###########

#Cytotoxicity test with AddModuleScore
#Cytotoxic geneset from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/cards/BIOCARTA_TCYTOTOXIC_PATHWAY)
cytotoxic_genes <- c("THY1","CD3G","ICAM1","ITGB2","CD3D","CD3E","CD247","ITGAL","CD8A","CD28","PTPRC","CD2")
Idents(object = combined.object.final ) <- "seurat_clusters"
combined.object.final <- AddModuleScore(combined.object.final,
                       features = list(cytotoxic_genes),
                       name="Cytotoxicity")
#Plot scores
#head(combined.object.final[[]]) # look at metadata
FeaturePlot(combined.object.final,
            features = "Cytotoxicity1", label = TRUE, repel = TRUE, pt.size = 0.6) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


DotPlot(combined.object.final, features = cytotoxic_genes) + RotatedAxis()


#Cluster 12 check
#Check cluster 12 (which is enriched for MT genes) for cell death from MSigDB
cell_death_genes <- c("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(cell_death_genes),
                                        name="Cell_death")
#Plot scores
head(combined.object.final[[]]) # look at metadata
FeaturePlot(combined.object.final,
            features = "Cell_death1", label = TRUE, repel = TRUE, pt.size = 0.6) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#DotPlot(combined.object.final, features = cell_death_genes) + RotatedAxis()


#Check cluster 12 (which is enriched for MT genes) for cell death from GO term cell death (http://amigo.geneontology.org/amigo/term/GO:0008219)
cell_death_GO_genes <- c("BCL7C","CCN1","CIAPIN1","SOX9","LGALS12","POLB","PYGL","S100A9","YBX3","HSPA5","ZFAND6","PHLDA2","GSN","CLC","CIB1","MYD88","AEN","SMO","PRUNE2","HELLS","TMEM14A","TRAF5","LMBR1L","DNM1L","P2RX4","DYNLL1","FAM3B","FNIP1","UBE2B","SLFN12","GPR65","HOXA13","GPER1","HMGB2","HIPK1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(cell_death_GO_genes),
                                        name="GO_Cell_death")
#Plot scores
head(combined.object.final[[]]) # look at metadata
p19 <- FeaturePlot(combined.object.final,
            features = "GO_Cell_death1", label = TRUE, repel = TRUE, pt.size = 0.6) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p19
ggsave("./figures/UMAP_GO_celldeath_geneset.pdf",p19, units = "cm", width = 10, height = 10)


#Violin plots and Featureplots (UMAP)

#Check Cytotoxicity (Positive regulation of T cell mediated cytotoxicity: GO:0001916, http://amigo.geneontology.org/amigo/term/GO:0001916)
cytotoxicity_GO_genes <- c("HLA-C","CD1D","CD1E","SLC22A13","HLA-F","MR1","P2RX7","CD1A","ULBP3","XCL1","HFE","ULBP1","ULBP2","HLA-A","HLA-DRB1","RAET1E","B2M","NECTIN2","TAP2","PTPRC","IL12B","IL12A","RAET1L","IL23A","HLA-G","CD1C","CD1B","MICA","MICB","HLA-E","HLA-DRA","HLA-H","HLA-B","PVR","STX7","FADD","RAET1G","CYRIB","IL23R","IL12RB1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(cytotoxicity_GO_genes),
                                        name="GO_cytotoxicity")

# Create a feature plot and a violin plot side by side
p22 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_cytotoxicity1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_cytotoxicity1", 
             group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
             + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5))

p22
ggsave("./figures/Feature_and_Violinplot_Cytotoxicity_high_vs_low_salt.pdf",p22, units = "cm", width = 27, height = 15)


p22_leiden <- VlnPlot(object = combined.object.final, features = "GO_cytotoxicity1", 
        group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p22_leiden
ggsave("./figures/Violinplot_Cytotoxicity_leiden.pdf",p22_leiden, units = "cm", width = 25, height = 15)


#Check Activation (Positive regulation of T cell activation: GO:0050870, http://amigo.geneontology.org/amigo/term/GO:0050870)
activation_GO_genes <- c("ZP4","SMARCB1","IL6ST","CD46","RASGRP1","TESPA1","CD274","CD74","THY1","CD160","DNAJA3","NFKBID","TNFSF14","LILRB1","PRKCQ","PIK3CA","SELENOK","RARA","HLX","HLA-DQB2","JAK3","ICOS","ARID1B","RAG1","SLC7A1","CCR7","LILRB4","HSPD1","SMARCD1","BCL10","ABL1","ABL2","DPP4","KLHL25","DHPS","CD28","NCK2","PTPN11","PNP","SOCS5","DOCK8","ICOSLG","MALT1","TMIGD2","CD1D","PCK1","TNFSF13B","GLI3","GLI2","CD86","PTPN22","VNN1","ZBTB1","HLA-DRB3","IHH","HES1","ZMIZ1","IL6","VAV1","HHLA2","PYCARD","CD27","RASAL3","TRAF6","IGFBP2","KITLG","BTN2A2","CD80","CARD11","BCL6","RPS3","NKAP","CCDC88B","CD55","CCL19","SHH","SHB","KLRK1","FYN","LCK","MDK","PHF10","KAT5","CLECL1P","IL7R","SMARCD2","BAD","ZBTB16","PRKCZ","TNFRSF14","AKT1","EFNB3","IGF2","GPAM","ZAP70","ZP3","NFKBIZ","NCKAP1L","WNT10B","IL2RA","EPO","SMARCC1","DUSP10","SYK","CCR2","CD5","IL15","XCL1","EFNB2","ACTL6A","RHOA","HLA-DPB1","HLA-A","PBRM1","HLA-DQA2","HLA-DRB1","HLA-DQB1","SMARCA4","SMARCA2","B2M","MAP3K8","TNFSF9","LGALS9","CSK","HSPH1","NOD2","CD4","IL21","HLA-DRB5","STAT5B","HAVCR2","CORO1A","AGER","SIRPB1","AP3B1","CBFB","IL1B","IL1A","IFNG","FOXP3","LILRB2","CCL21","AIF1","HLA-DOA","IL2RG","RHOH","CD47","HLA-DMB","HLA-DMA","LEP","IL1RL2","RUNX3","PTPRC","IL12B","IL12A","CD70","TGFBR2","RUNX1","CD3E","CD83","PIK3R6","PPP3CA","SPTA1","NLRP3","EFNB1","ADA","RIPK2","SOX4","VSIR","LYN","YES1","VCAM1","CD209","ITPKB","PTPN6","CR1","PDCD1LG2","TYK2","EGR3","BRD7","SASH3","IL23A","XBP1","HLA-G","CD81","GATA3","FCHO1","SLAMF1","ARID1A","JAK2","IL4","HLA-DPA1","SOX13","CD6","FLOT2","IGF1","EBI3","AMBRA1","SMARCE1","HLA-DOB","HLA-DRB4","HLA-E","IL2","TNFSF11","ANXA1","ACTL6B","NCK1","SIRPG","HLA-DQA1","HLA-DRA","SART1","TNFSF4","VTCN1","IL7","IL18","CD24","SMARCC2","BMI1","IL36B","CAV1","IL4I1","CD276","TNFRSF13C","FADD","TFRC","CD40LG","HMGB1","IL4R","ZBTB7B","SRC","SOX12","ARID2","AP3D1","LGALS1","ACTB","SMARCD3","CYRIB","ADAM8","SIRPA","SOCS1","CCL5","CCL2","LEF1","FOXO3","IL23R","IL12RB1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(activation_GO_genes),
                                        name="GO_activation")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p23 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_activation1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_activation1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p23
ggsave("./figures/Feature_and_Violinplot_Activation_high_vs_low_salt.pdf",p23, units = "cm", width = 27, height = 15)

#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
Activation_sub <- combined.object.final_sub@meta.data$GO_activation1
Activation_test_result <- wilcox.test(Activation_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
Activation_p_value <- Activation_test_result$p.value 
Activation_p_value #1.971527e-151


p23_leiden <- VlnPlot(object = combined.object.final, features = "GO_activation1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p23_leiden
ggsave("./figures/Violinplot_Activation_leiden.pdf",p23_leiden, units = "cm", width = 25, height = 15)


#Check Effector function (Reference Effector genes: "effector genes" in Debdas list)
effctor_GO_genes <- c("GZMB","GZMK","IFNG","EOMES","ITGAD","ITGAX","ITGB7","CXCR3","CCR5")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(effctor_GO_genes),
                                        name="GO_effector_function")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p24 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_effector_function1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_effector_function1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p24
ggsave("./figures/Feature_and_Violinplot_Effector_function_high_vs_low_salt.pdf",p24, units = "cm", width = 27, height = 15)


p24_leiden <- VlnPlot(object = combined.object.final, features = "GO_effector_function1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p24_leiden
ggsave("./figures/Violinplot_Effector_function_leiden.pdf",p24_leiden, units = "cm", width = 25, height = 15)



#Check TCR signalling function (GO:0050862 )
TCR_sig_GO_genes <- c("RELA","CCR7","BCL10","CD81","TESPA1","RAB29","CARD11","RPS3","TRAT1","IKBKG","LCK","CD226","PRKD2","NECTIN2","ADA","KCNN4")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(TCR_sig_GO_genes),
                                        name="GO_TCR_signalling")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p25 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_TCR_signalling1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_TCR_signalling1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p25
ggsave("./figures/Feature_and_Violinplot_TCR_signalling_high_vs_low_salt.pdf",p25, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
TCRsig_sub <- combined.object.final_sub@meta.data$GO_TCR_signalling1
TCRsig_test_result <- wilcox.test(TCRsig_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
TCRsig_p_value <- TCRsig_test_result$p.value 
TCRsig_p_value #3.53706e-11

p25_leiden <- VlnPlot(object = combined.object.final, features = "GO_TCR_signalling1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p25_leiden
ggsave("./figures/Violinplot_TCR_signalling_leiden.pdf",p25_leiden, units = "cm", width = 25, height = 15)


#Check stess activated MAPK signalling function (GO:0051403)
MAPK_GO_genes <- c("MYD88","MAPK10","ARHGEF6","LGALS9","ERCC6","SH2D3C","WNT5A","NOD2","MAP4K3","CASR","IRAK1","AGT","MAPK11","MAPK8IP2","LRRK2","MAPK1","ZFP36L1","MINK1","DUSP10","MAP4K1","ZFP36","TNF","CDC42EP5","DUSP9","PAFAH1B1","IL1B","MAP3K5","NOX1","IKBKB","MAP3K12","MAP4K2","MAP2K3","MAPK8IP1","PBK","ATF2","MAP2K4","MAPK9","MAPK8","TNFSF11","MAP2K7","TNFRSF19","TLR4","TLR3","STRADB","IRAK4","MAPK13","STK3","PTGER4","CRYAB","SH2D3A","RIPK2","NPHS1","GPS1","TLR7","SMAD3","MAP3K20","CARD9","MAP3K11","MAPK14","CCM2","MAP2K6","MAP3K7","MAPK3","DAXX","MAPKAPK2","HACD3","NFKB1","GPS2","CRKL","TRIB1","MAP3K10","MAP3K13","NOD1","TAOK2")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(MAPK_GO_genes),
                                        name="GO_Stress_activated_MAPK_cascade")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p26 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p26
ggsave("./figures/Feature_and_Violinplot_MAPK_cascade_high_vs_low_salt.pdf",p26, units = "cm", width = 27, height = 15)

#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
MAPK_sub <- combined.object.final_sub@meta.data$GO_Stress_activated_MAPK_cascade1
MAPK_test_result <- wilcox.test(MAPK_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
MAPK_p_value <- MAPK_test_result$p.value #9.779305e-157

p26_leiden <- VlnPlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p26_leiden
ggsave("./figures/Violinplot_MAPK_cascade_leiden.pdf",p26_leiden, units = "cm", width = 25, height = 15)


#Check positive regulation of calcineurin-NFAT signaling cascade (GO:0070886)
NFAT_GO_genes <- c("CIB1","PPP3CC","PTBP1","SPPL3","ERBB3","STIMATE","PPP3CB","PPP3R1","TNF","IGF1","CEFIP","PPP3R2","LACRT","PPP3CA","SLC9A1","AKAP6","CAMTA1","AKAP5","CHP2","CHERP","LMCD1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(NFAT_GO_genes),
                                        name="GO_NFAT_signalling")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p27 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p27
ggsave("./figures/MA_Feature_and_Violinplot_NFAT_signalling_high_vs_low_salt.pdf",p27, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
NFAT_sub <- combined.object.final_sub@meta.data$GO_NFAT_signalling1
NFAT_test_result <- wilcox.test(NFAT_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
NFAT_p_value <- NFAT_test_result$p.value 
NFAT_p_value #1.887247e-223

p27_leiden <- VlnPlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p27_leiden
ggsave("./figures/MA_Violinplot_NFAT_cascade_leiden.pdf",p27_leiden, units = "cm", width = 25, height = 15)


#Check TRM signature from Gustavo's Heatmap from Science Immunology Paper 2022 (only upregulated)
TRM_GO_genes <- c("XIST","UBC","LGALS3","MT-CO2","VIM","ANKRD28","RGS1","RGCC","HSPA1B","MT-ND4","HSP90AB1","PPP1R15A")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(TRM_GO_genes),
                                        name="GO_TRM")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p28 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_TRM1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_TRM1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p28
ggsave("./figures/MA_Feature_and_Violinplot_TRM_high_vs_low_salt.pdf",p28, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
TRM_sub <- combined.object.final_sub@meta.data$GO_TRM1
TRM_test_result <- wilcox.test(TRM_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
TRM_p_value <- TRM_test_result$p.value 
TRM_p_value #6.731356e-57

p28_leiden <- VlnPlot(object = combined.object.final, features = "GO_TRM1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p28_leiden
ggsave("./figures/MA_Violinplot_TRM_leiden.pdf",p28_leiden, units = "cm", width = 25, height = 15)



#Check Exhaustion geneset (Genelist 1 fromScience paper)
Exhaustion_GO_genes <- c("TOX","SNX9","ITGAE","HAVCR2","IRF4","UBE2F","LAG3","STMN1","NDFIP2","ENTPD1","PRDM1","CD63","PDCD1","CD2BP2","FKBP1A","CTLA4","RAB27A","TPI1","CD38","DUSP4","CDCA8","TIGIT","PHLDA1","NCAPG2","VCAM1","ITM2A","CDKN3","CD27","IFI35","SNAP47","ISG15","IGFLR1","STAT3","RAD51","WARS","CCNB1","SYNGR2","BUB1","GBP2","SIRPG","LYST","SEMA4A","BST2","CXCR6","PARK7","ACP5","CCL4L2","HLA-DRA","RGS2","FCRL3","NAB1","OSBPL3","ID3","ICOS","CCR5","FAM3C","GOLIM4","PTPN11","ACP5","CKS2","FUT8","GALM","HLA-DMA")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(Exhaustion_GO_genes),
                                        name="Exhaustion")

head(combined.object.final[[]])

# Create a feature plot and a violin plot side by side
p29 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "Exhaustion1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "Exhaustion1", 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p29
ggsave("./figures/MA_Feature_and_Violinplot_Exhaustion_high_vs_low_salt.pdf",p29, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
Exhaustion_sub <- combined.object.final_sub@meta.data$Exhaustion1
Exhaustion_test_result <- wilcox.test(Exhaustion_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
Exhaustion_p_value <- Exhaustion_test_result$p.value 
Exhaustion_p_value #7.609517e-76

p29_leiden <- VlnPlot(object = combined.object.final, features = "Exhaustion1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p29_leiden
ggsave("./figures/MA_Violinplot_Exhaustion_leiden.pdf",p28_leiden, units = "cm", width = 25, height = 15)



#Make violinplot and UMAP for specific genes of interest (ICOS, PDCD1, CTLA4, ITGAE)
features <- c("ICOS", "PDCD1", "CTLA4", "ITGAE")

p30 <- FeaturePlot(object = combined.object.final, features = features) 
p30
ggsave("./figures/MA_Featureplot_genes_of_interest_high_vs_low_salt.pdf",p30, units = "cm", width = 20, height = 15)

#ICOS
p31 <- VlnPlot(object = combined.object.final, features = c("ICOS"), 
          group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p31
ggsave("./figures/MA_Violinplot_ICOS_gene_high_vs_low_salt.pdf",p31, units = "cm", width = 13, height = 13)

#PDCD1
p32 <- VlnPlot(object = combined.object.final, features = c("PDCD1"), 
               group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p32
ggsave("./figures/MA_Violinplot_PDCD1_gene_high_vs_low_salt.pdf",p32, units = "cm", width = 13, height = 13)

#CTLA4
p33 <- VlnPlot(object = combined.object.final, features = c("CTLA4"), 
               group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p33
ggsave("./figures/MA_Violinplot_CTLA4_gene_high_vs_low_salt.pdf",p33, units = "cm", width = 13, height = 13)


#ITGAE
p34 <- VlnPlot(object = combined.object.final, features = c("ITGAE"), 
               group.by = "orig.ident") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p34
ggsave("./figures/MA_Violinplot_ITGAE_gene_high_vs_low_salt.pdf",p34, units = "cm", width = 13, height = 13)







library(SeuratDisk)

Convert("./output_data/high_salt/high_salt.h5ad", ".h5seurat")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat("./output_data/high_salt/high_salt.h5Seurat")


