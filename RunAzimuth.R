## Usage: Azimuth.R <Query_Seurat> <Reference_Seurat> <Ref_MetaColumn> <output_Directory>
## Run Azimuth on SUV organoid ATAC data and try to map it to Trevino clusters
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
args <- commandArgs(trailingOnly=TRUE)
#args <- c("2Months.RDS", "CZI_RNA.RDS", "NameBroad", "org2mo_to_CZI")
queryName <- args[1]
refName <- args[2]
refMeta <- args[3]
outDir <- args[4]

ref <- readRDS(refName)
set.seed(11)
for(i in unique(ref[[refMeta]])[,1]){
  mycells <- c(mycells, sample(which(ref[[refMeta]]==i), min(450,sum(ref[[refMeta]]==i))))
}
ref <- subset(ref, cells=mycells)
ref$Sample <- tail(strsplit(refName, "/")[[1]], n=1)
query <- readRDS(queryName)
query$Sample <- tail(strsplit(refName, "/")[[1]], n=1)
if(!"NameBroad" %in% names(query@meta.data)){
  query$NameBroad <- query$FinalName
}
ref <- RunUMAP(ref, reduction="pca", dims=1:30, return.model=TRUE)
Idents(ref) <- ref[[refMeta]]
query <- SCTransform(query, conserve.memory = T)
#ref <- FindTopFeatures(ref, min.cutoff=10)
#ref <- RunUMAP(ref, reduction="pca", dims=1:30, return.model=T)
print("Finding achors!")
anchors1 <- FindTransferAnchors(
  reference = ref,
  query = query,
  reference.reduction = "pca",
  reduction="pcaproject",
  dims = 1:30,
  k.anchor=7,
  max.features=250,
  k.filter=NA
)
print("Mapping query onto reference!")
query <- MapQuery(
  anchorset = anchors1,
  reference = ref,
  query = query,
  refdata = ref[[refMeta]][,1],
  reference.reduction = "pca",
  new.reduction.nam = "ref.pca",
  reduction.model = "umap"
)
dir.create(outDir)
print("Plotting!")
p1 = DimPlot(query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
p2 = FeaturePlot(query, reduction = "ref.umap", features = "predicted.id.score")
png(paste0(outDir,"/Azimuth_Predicted.png"))
p1 + p2
dev.off()
Idents(query) = "predicted.id"
query$confident.pred = query$predicted.id
query$confident.pred[which(query$predicted.id.score < .7)] = 'unsure'
library(scales)
cols = hue_pal()(length(table(query$confident.pred))-1)
names(cols) <- names(table(query$confident.pred))[-length(table(query$confident.pred))]
cols[["unsure"]] = "gray"
p1 <- DimPlot(query, group.by="confident.pred", reduction='ref.umap', pt.size=.1, cols=cols, label=T) + ggtitle("prediction score over .7 only") + NoLegend()
p2 <- DimPlot(query, group.by="NameBroad", reduction="ref.umap", label=F)+ ggtitle("Original Subtypes")
png(paste0(outDir,"/Azimuth_Confident.png"))
p1+p2
dev.off()
saveRDS(query, paste0(outDir,"/Azimuth_Results_Subtype.RDS"))

myTab <- table(query$NameBroad, query$predicted.id)
saveRDS(myTab, paste0(outDir,"/Table_orig_v_predicted.RDS"))
myTab <- t(t(myTab)/apply(myTab,2,sum)*100)
myTab <- myTab/apply(myTab,1,sum)*100
library(reshape2)
myMelt <- melt(myTab)
myMelt$Var1 <- factor(myMelt$Var1)
myMelt$Var2 <- factor(myMelt$Var2)

ggplot(myMelt, aes(x=Var2,
                   y=Var1,
                   size=value, color=value)) +
  geom_point() +
  scale_color_gradient(low="grey", high="purple4") +
  NoLegend() +
  xlab("CZI Cell Types") +
  ylab("Organoid Cell Types") +
  ggtitle(outDir)
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold", size=14, hjust=0.5))


ggsave(paste0(outDir,"/DotPlot.pdf"), width=6, height=5)

saveRDS(query@meta.data, paste0(outDir, "/QueryMetadata.RDS"))

