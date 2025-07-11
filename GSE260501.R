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

#Dimplot
DimPlot(UTI, reduction = "umap", label = TRUE)





# Dimplot
DimPlot(UTI, group.by = "
", label = TRUE, reduction = "umap")


