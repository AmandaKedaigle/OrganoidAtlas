## Usage: DA_Analysis_Green.R <QC_seurat_object> <output_basename> <celltype>
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

args <- c("Trevino_Shared_SUV3.RDS", "DARS_Trevino3", "Early_RG")
args <- commandArgs(trailingOnly=TRUE)

seur <- readRDS(args[1])
DefaultAssay(seur) <- "ATAC"
Idents(seur) <- seur$CellSubType

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(seur) <- annotations

celltype <- args[3]
celltype <- gsub("_"," ",celltype)
print(paste("Running:", celltype, "in",strsplit(args[1],"_")[[1]][1]))
myDARs <- FindMarkers(seur, ident.1=celltype, min.pct=0.05, logfc.threshold = 0, test.use = 'LR', latent.vars = 'nCount_ATAC')
myDARs$peak = rownames(myDARs)
myDARs$gene = NULL
closest_genes = ClosestFeature(seur, regions = rownames(myDARs))
myDARs$closest_gene = closest_genes$gene_name
myDARs$closest_gene_distance = closest_genes$distance
celltype <- gsub(" ","_",celltype)
saveRDS(myDARs, paste1(args[2],"_",celltype,".RDS"))
