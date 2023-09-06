# Loading test data.
library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")

# Loading reference data with Ensembl annotations.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData()

table(hpca$main_types)
table(hpca$types)
table(blueprint_encode$main_types)
table(blueprint_encode$types)

# Performing predictions.
library(SingleR)
predictions <- SingleR(test=new.data, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)

#visualize
predictions$scores[1:10,]

#heatmap
plotScoreHeatmap(predictions)

#violin diatribution
plotDeltaDistribution(predictions)


all.markers <- metadata(predictions)$de.genes
beta.markers <- unique(unlist(all.markers$beta))
new.data$labels <- predictions$labels

library(scater)
plotHeatmap(new.data, order_columns_by="labels", features=beta.markers)



################################################################################

library(scRNAseq)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 


adult <- readRDS("./output/adult_results.rds") #skin adult
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(list(counts=counts))
sce

pretend.cell.labels <- sample(letters, ncol(counts), replace=TRUE)
pretend.gene.lengths <- sample(10000, nrow(counts))

sce <- SingleCellExperiment(list(counts=counts),
                            colData=DataFrame(label=pretend.cell.labels),
                            rowData=DataFrame(length=pretend.gene.lengths),
                            metadata=list(study="GSE111111")
)
sce

#-----------------------------------------------------------------------------------------------------------------#
# Biological annotation of clusters from blood sample using SingleR
#-----------------------------------------------------------------------------------------------------------------#

library(Seurat)

# Create the SingleR object
singlec <- CreateSinglerObject(GetAssayData(adult, slot = "data"), 
                               annot = NULL, project.name = "blood_ann", min.genes = 200, #500
                               technology = "10X", species = "Human", citation = "",
                               ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                               fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
                               reduce.file.size = T, numCores = 4)#SingleR.numCores)
#saveRDS(singlec, file = paste0(dir.figs.bd, "_singler_blood_clust.rds"))




###################################################################################

#Fetal
train_single <- trainSingleR(GetAssayData(adult, slot = 'data'),adult@active.ident, genes = 'de',  sd.thresh = 1,
                             de.method = c("classic", "wilcox", "t"),
                             de.n = NULL,
                             de.args = list(),
                             aggr.ref = FALSE,
                             aggr.args = list(),
                             recompute = TRUE,
                             restrict = NULL,
                             assay.type = "logcounts",
                             check.missing = TRUE  )

sr.bd <- classifySingleR(
  GetAssayData(adult, slot = 'data'),
  train_single,
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  sd.thresh = NULL,
  prune = TRUE,
  assay.type = "logcounts",
  check.missing = TRUE)

saveRDS(sr.bd, file = "./output/SingleR/adult_singleR_results.rds")
# Load Seurat object
sk <- readRDS("./output/fetal_qc.rds")

sr.sk <- readRDS('fetal_singleR_results.rds')


# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
sr.sk$meta.data$orig.ident <- sk@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.sk$meta.data$xy <- sk@reductions$umap@cell.embeddings[,1] # the tSNE coordinates (Embeddings(sk, reduction = "tsne"))
sr.sk$meta.data$xy.um <- sk@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(sk, reduction = "umap"))
sr.sk$meta.data$clusters <- sk@active.ident # the Seurat clusters, if 'clusters' not provided (Idents(sk))


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
                                                                   "mB cell", "other", "mB cell", "other", "nB cell", "NK cell", "Treg"))

#####################################################################################

#sce_read_h5(file = "../data/Paper/GSE156972_raw_gene_bc_matrices_h5.h5", assay.name = "RNA")

library(DropletUtils)
sce10x_fetus <- read10xCounts("../data/Paper/GSE156972_raw_gene_bc_matrices_h5.h5")
sce10x <- read10xCounts("../data/EX0004/raw_feature_bc_matrix/")

#reference data from Human cell atlas 
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

#SingleR new predictions
library(SingleR)
predictions_new <- SingleR(test=sce10x, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions_new$labels)


##################################################################################

adult <- readRDS("./output/adult_results.rds") #skin adult
fetus <- readRDS("./output/fetal_result.rds") #fetal skin

library(Seurat)
library(SingleR)

# Generating SingleR+Seurat object
singler = CreateSinglerSeuratObject(GetAssayData(adult, slot = "data"),
                                    annot = NULL, project.name = 'adult_ann', min.genes = 200,
                                    technology = "10X", species = "Human", citation = "",
                                    ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                    fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
                                    reduce.file.size = T, numCores = 4)

singler <- readRDS('./GA_AN0146_sc_skin_blood_10x_singler_skin.rds')

singler.new = convertSingleR2Browser(singler)
#saveRDS(singler.new,file=paste0(singler.new@project.name,'.rds')




############################################ New try #####################################

# Loading reference data with Ensembl annotations.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData()

table(hpca$main_types)
table(hpca$types)
table(blueprint_encode$main_types)
table(blueprint_encode$types)

adult <- readRDS("./output/adult_results.rds") #skin adult
fetus <- readRDS("./output/fetal_result.rds") #fetal skin

#Run SingleR
library(Seurat)
library(SingleR)
library(ggplot2)
library(cowplot)

counts <- GetAssayData(adult)

singler <- CreateSinglerObject(counts=counts,
                               project.name="adult_ann", # choose
                               min.genes = 200, # ignore cells with fewer than 200 transcripts
                               technology = "10X", # choose
                               species = "Human",
                               citation = "EX0004 healthy adult", # choose
                               ref.list = list(hpca=hpca, bpe=blueprint_encode),
                               normalize.gene.length = FALSE,        # needed for full-length platforms (e.g. smartseq)
                               variable.genes = "de",  # see vignette
                               fine.tune = FALSE, # TRUE would take very long
                               reduce.file.size = TRUE, # leave out less-often used fields 
                               do.signatures = FALSE,
                               do.main.types = TRUE,
                               numCores = SingleR.numCores)

#inspect SingleR results
show(names(singler$singler))

for (ref.set in names(singler$singler) ) {
  types <- singler$singler[[ref.set]]$SingleR.single.main$labels[,1]
  cat("==== ", ref.set, ": ====\n")
  show(sort(table(types), decreasing=TRUE))
}

## For hpca and blueprint_encode also show the 
## detailed cell typings (as opposed to main_types results) : 
for (ref.set in c("hpca", "bpe") ) {
  types <- singler$singler[[ref.set]]$SingleR.single$labels[,1]
  subrefset <- paste(ref.set, "subtypes", sep="_") 
  cat("==== ", subrefset, ": ====\n")
  show(sort(table(types), decreasing=TRUE))
}

#We will stick to the main_types from now on, for brevity. 
#In order to easily visualize the various classifications in the tSNE plots we have to add them to meta.data slot of the Seurat object
for (ref.set in names(singler$singler) ) {
  types <- singler$singler[[ref.set]]$SingleR.single.main$labels[,1]
  adult <- AddMetaData(adult,
                         metadata=types,
                         col.name=paste0(ref.set,"_type" ) )
}

## Check if this worked and get an impression of the concordance of classification
interesting.columns <- c("celltypes", "hpca_type", "bpe_type")

## repeat the following a few times:
random.rows <- sort(sample(ncol(adult), size=20))
adult@meta.data[ random.rows,  interesting.columns]

#Plotting in UMAP
panel.labels <- c('hpca','bpe') #shorthand labels

p1 <- DimPlot(adult, group.by="seurat_clusters", reduction="umap")
p2 <- DimPlot(adult, group.by='hpca_type', reduction="umap")
p3 <- DimPlot(adult, group.by='bpe_type', reduction="umap")
plot_grid(p1, p2, p3, nrow=2, ncol=2, labels=panel.labels)


#Highlight the cell type of interest from the UMAP
findCells <- function(obj, column, values, name=NULL) {
  ## Given a Seurat OBJ, return a list with the names of the cells where
  ## the specified meta.data COLUMN equals any of the strings specified
  ## in VALUES (both must be characters or factors). Name of the list member
  ## must be specified using the NAME argument if length(values)>1
  stopifnot(is(obj, "Seurat"))
  stopifnot(is.character(column))
  stopifnot(column %in% names(obj@meta.data))
  col <- obj@meta.data[[column]]
  stopifnot(is.character(col) || is.factor(col))
  values <- unique(values)
  stopifnot(is.character(values) || is.factor(values))
  if (length(values)>1 && is.null(name))
    stop("findCells: specify a name to be used for the selection")
  if(is.null(name))
    name <- values
  stopifnot(is.character(name))
  rem <- setdiff(c(values), col)
  if(length(rem)>0)stop("findCells: requested value(s) never occurs in this column: ", rem)
  l <- list(colnames(obj)[ col %in% values ])
  names(l) <- name
  l
}         #find cells

## "Let's look at individual cell types
#Monocytes
p1 <- DimPlot(adult, group.by='hpca_type', reduction="umap",
              cells.highlight=findCells(adult, 'hpca_type', 'Monocyte'))
p2 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'Monocytes'))
plot_grid(p1, p2, nrow=2, ncol=2, labels=paste(panel.labels, "Monocyte"))

#DC
p1 <- DimPlot(adult, group.by='hpca_type', reduction="umap",
              cells.highlight=findCells(adult, 'hpca_type', 'DC'))
p2 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'DC'))
plot_grid(p1, p2, nrow=2, ncol=2, labels=paste(panel.labels, "DC"))

#B-cells
p1 <- DimPlot(adult, group.by='hpca_type', reduction="umap",
              cells.highlight=findCells(adult, 'hpca_type', 'B_cell'))
p2 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'B-cells'))
plot_grid(p1, p2, nrow=2, ncol=2, labels=paste(panel.labels, "B-cells"))

#NK cells
p1 <- DimPlot(adult, group.by='hpca_type', reduction="umap",
              cells.highlight=findCells(adult, 'hpca_type', 'NK_cell'))
p2 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'NK cells'))
plot_grid(p1, p2, nrow=2, ncol=2, labels=paste(panel.labels, "NK-cells"))

#T cells
p1 <- DimPlot(adult, group.by='hpca_type', reduction="umap",
              cells.highlight=findCells(adult, 'hpca_type', 'T_cells'))
p2 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'CD4+ T-cells'))
p3 <- DimPlot(adult, group.by='bpe_type', reduction="umap",
              cells.highlight=findCells(adult, 'bpe_type', 'CD8+ T-cells'))
plot_grid(p1, nrow=2, ncol=2, labels=paste(panel.labels, "T-cells"))
plot_grid(p2, nrow=2, ncol=2, labels=paste("bpe", "T-cells"))
plot_grid(p3, nrow=2, ncol=2, labels=paste("bpe", "T-cells"))


#Discrepancies
DimPlot(adult, split.by='seurat_clusters', group.by='hpca_type', reduction="umap")

