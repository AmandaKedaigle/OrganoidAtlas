## Usage: DA_Analysis_Green.R <QC_seurat_object> <Grouping_Variable> <output_basename> <Cluster>
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(openxlsx)

args <- commandArgs(trailingOnly=TRUE)

seur <- readRDS(args[1])
groupingVar <- args[2]
outName <- args[3]
celltype <- args[4]

#Differentially accessible peaks -----
DefaultAssay(seur) <- "ATAC"
Idents(seur) <- seur@meta.data[,groupingVar]

print(paste("Starting analysis of:",celltype))
myDARs <- FindMarkers(seur, ident.1=celltype, min.pct=0.05, logfc.threshold = 0, test.use = 'LR', latent.vars = 'nCount_ATAC')
myDARs$peak = rownames(myDARs)
myDARs$gene = NULL
closest_genes = ClosestFeature(seur, regions = rownames(myDARs))
myDARs$closest_gene = closest_genes$gene_name
myDARs$closest_gene_distance = closest_genes$distance
print(celltype)
print(nrow(myDARs))

saveRDS(myDARs, paste0(outName, "_", celltype, ".RDS"))

