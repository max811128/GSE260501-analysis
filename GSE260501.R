library(Seurat)
library(tidyverse)

UTI<- Read10X(data.dir = "F:/GSE260501_RAW/")
UTI<- CreateSeuratObject(counts = UTI, project = "UTI", min.cells = 3, min.features = 200)
UTI<- PercentageFeatureSet(UTI, pattern = "^mt-", col.name = "percent.mt")
UTI<- subset(UTI, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
#VlnPlot
VlnPlot(UTI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

UTI <- NormalizeData(UTI)
UTI <- FindVariableFeatures(UTI, selection.method = "vst", nfeatures = 2000)

UTI <- ScaleData(UTI, features = VariableFeatures(UTI))

UTI <- RunPCA(UTI, features = VariableFeatures(UTI))
UTI <- RunUMAP(UTI, reduction = "pca", dims = 1:30)
UTI <- FindNeighbors(UTI, reduction = "pca", dims = 1:30)
UTI <- FindClusters(UTI, resolution = 0.25)

# Dimplot
DimPlot(UTI, reduction = "umap", label = TRUE)

## SingleR annotation
library(SingleR)
library(SingleCellExperiment)
library(celldex)
ref<-celldex::MouseRNAseqData()
data.input <- GetAssayData(UTI, assay = "RNA", layer = "data")
meta <- data.frame(row.names = colnames(UTI))
sce <- SingleCellExperiment(assays = list(logcounts = data.input), colData = meta)
results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
UTI$SingleR.label <- results$labels

# Dimplot
DimPlot(UTI, group.by = "SingleR.label", label = TRUE, reduction = "umap")


UTI <- RenameIdents(UTI, `0` = "Monocytes",`1` = "Monocytes", `2` = "Neutrophils", `3` = "T cells", `4` = "Monocytes", `5` = "Macrophages", `6` = "NK cells", `7` = "T cells", `8` = "T cells", `9` = "Macrophages", `10` = "Dendritic cells")
DimPlot(UTI, label = TRUE, reduction = "umap")

### Cell Neumbers
cluster_counts <- table(Idents(UTI))
print(cluster_counts, split.by = "orig.ident")

# Expression of NETosis associated mRNA
# Violin plots can also be split on some variable. Simply add the splitting variable to object
# metadata and pass it to the split.by argument
VlnPlot(UTI, features = c("Il1b"))
VlnPlot(UTI, features = c("Gla"))
VlnPlot(UTI, features = c("Anxa1"))
VlnPlot(UTI, features = c("Egr1"))

# CellChat
library(CellChat)
library(reticulate)
library(patchwork)
library(circlize)
meta <- data.frame(group = UTI$SingleR.label, row.names = names(UTI$SingleR.label))
data.input <- GetAssayData(UTI, assay = "RNA", layer = "data")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Load the ligand-receptor interaction database
CellChatDB.mouse <- CellChatDB.mouse
# ShowDatabaseCategory(CellChatDB.mouse)
CellChatDB.use <- CellChatDB.mouse
cellchat@DB <- CellChatDB.use

# Subset and pre-processing the expression data 
# Subset the expression data to use less RAM
cellchat <- subsetData(cellchat)

# Pre-processing the expression data
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight

#Circle plot
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
