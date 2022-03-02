## Run DecontX locally
library(Seurat)
library(celda)
library(gridExtra)
library(ggplot2)

seur <- readRDS("6mo_harm_111120.rds")
age <- "6 months"
seur$Organoid <- paste(seur$dataset, seur$orig.ident, sep="_")
seur$Age <- "6mo"

Decontaminated <- decontX(x = seur@assays$RNA@counts, batch=seur$Organoid)
saveRDS(Decontaminated, "DecontX_output_6mo.rds")
seur$Contamination <- Decontaminated$contamination
#p1 <- VlnPlot(seur, "Contamination", group.by="FinalName", pt.size = 0) + ggtitle("Estimated ambient RNA by cell type") + NoLegend()
#p2 <- VlnPlot(seur, "Contamination", group.by="Organoid", pt.size = 0) + ggtitle("Estimated ambient RNA by organoid") + NoLegend()
#grid.arrange(p1, p2, ncol=2)

seur[["DecontX"]] <- CreateAssayObject(counts=Decontaminated$decontXcounts)
DefaultAssay(seur) <- "DecontX"
saveRDS(seur, "6mo_ambientRemoved.rds")

# p1 <- DimPlot(seur, group.by="FinalName", label=T) + NoLegend() + ggtitle(paste("Original UMAP\nat", age))
# p2 <- FeaturePlot(seur, "Contamination") + ggtitle(paste("Estimated ambient RNA content\nat", age))
# seur <- SCTransform(seur, conserve.memory=T)
# seur <- RunPCA(seur, npcs=30)
# seur <- RunUMAP(seur, dims=1:30)
# p3 <- DimPlot(seur, group.by="FinalName", label=T) + NoLegend() + ggtitle(paste("Ambient RNA removed UMAP\nat",age))
# p1+p2+p3

  
  
  
  