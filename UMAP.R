library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


# load the data
# Data format
# row: features
# column: Nuclei
mydata <- read.csv("/home/j.fermin/EEL5934-cell-map/new_feature_transpose.csv", header = TRUE)

# exclude first column
mydata <- mydata[, -1]


myseurat <- CreateSeuratObject(counts = mydata)

# Apply Log Normalization
myseurat <- NormalizeData(myseurat)


# Find the most variable features
myseurat <- FindVariableFeatures(myseurat,selection.method = "vst", nfeatures = 15)


myseurat <- ScaleData(myseurat)

# apply PCA
myseurat <- RunPCA(myseurat, features = VariableFeatures(object = myseurat))


# Do KNN
myseurat <- FindNeighbors(myseurat, dims = 1:14)#, annoy.metric = 'cosine')


# Unsupervised Clustering
myseurat <- FindClusters(myseurat, resolution = 0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(myseurat), 5)

# Perform UMAP
myseurat <- RunUMAP(myseurat, dims = 1:14)

# Plot UMAP
DimPlot(myseurat, reduction = "umap")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(myseurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(myseurat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
myseurat.markers <- FindAllMarkers(myseurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


x <- myseurat.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(myseurat, features = x$gene[1:4])

FeaturePlot(myseurat, features = x$gene[5:8])


new.cluster.ids <- levels(myseurat) 
names(new.cluster.ids) <- levels(myseurat)

myseurat <- RenameIdents(myseurat, new.cluster.ids)
DimPlot(myseurat, reduction = "pca", label = TRUE, pt.size = 0.5)



plot1 <- DimPlot(myseurat, reduction = "umap", label = TRUE, pt.size = 0.1)
plot2 <- DimPlot(myseurat, reduction = "pca", label = TRUE, pt.size = 0.5)
plot1
plot2

# Get the UMAP coordinates and cluster ID
umap_coords <- as.data.frame(Embeddings(object = myseurat, reduction = "umap"))
cluster_id <- Idents(myseurat)

# Combine the data into a single data frame
umap_data <- cbind(umap_coords, cluster_id)

# Save the data to a CSV file
write.csv(umap_data, file = "/home/j.fermin/EEL5934-cell-map/umap_data.csv", row.names = FALSE)
