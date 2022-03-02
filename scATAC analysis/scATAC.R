library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(writexl)
library(patchwork)
library(reshape2)
library(ComplexHeatmap)

#Load cellranger output -------
counts <- Read10X_h5("filtered_peak_bc_matrix.h5") #from Agg CellRanger
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = 'fragments.tsv.gz',
  min.cells = 1
)

seur <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
seur$orig.ident = substr(rownames(seur@meta.data), nchar(rownames(seur@meta.data)), nchar(rownames(seur@meta.data)))

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(seur) <- annotations

#QC --------
# compute nucleosome signal score per cell
seur <- NucleosomeSignal(object = seur)

# compute TSS enrichment score per cell
seur <- TSSEnrichment(object = seur, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
seur$pct_reads_in_peaks <- seur$peak_region_fragments / seur$passed_filters * 100
seur$blacklist_ratio <- seur$blacklist_region_fragments / seur$peak_region_fragments

seur$high.tss <- ifelse(seur$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(seur, group.by = 'high.tss') + NoLegend()
seur$nucleosome_group <- ifelse(seur$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
FragmentHistogram(object = seur, group.by = 'nucleosome_group')
png("1m-atacQC.png", res=80, width=1000, height=400)
VlnPlot(seur,features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,ncol = 5)
dev.off()

#Set my QC thresholds based on above and on cellranger summaries
seur <- subset(seur, subset = peak_region_fragments > 1500 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 40 &
    #below from Signac recommendataions - nearly all our cells pass already
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

#Normalization & dim reduction --------
seur <- RunTFIDF(seur)
seur <- FindTopFeatures(seur, min.cutoff = 'q0') #using all peaks
seur <- RunSVD(seur)

#Check correlation of LSI components with sequencing depth
DepthCor(seur) #first component is super highly correlated so we leave it out downstream

#UMAP & cluster ------
seur <- RunUMAP(object = seur, reduction = 'lsi', dims = 2:30)
seur <- FindNeighbors(object = seur, reduction = 'lsi', dims = 2:30)
seur <- FindClusters(object = seur, verbose = FALSE, algorithm = 3)

png("umap.png", res=100, width=600, height=400)
DimPlot(object = seur, label = TRUE)
dev.off()
png("umap.per.org.png", res=100, width=900, height=300)
DimPlot(seur, split.by = "orig.ident", label=T) + NoLegend()
dev.off()

#Gene activity ----------
gene.activities <- GeneActivity(seur)
seur[['RNA']] <- CreateAssayObject(counts = gene.activities)
seur <- NormalizeData(seur, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(seur$nCount_RNA))

DefaultAssay(seur) <- 'RNA'
g.list = c("FOXG1","SNAP25","BCL11B","SATB2","GAD1","DLX2","VIM","PAX6","SPARCL1","HOPX","EOMES","TOP2A")
png('geneActivity.feature.plot.umap.png', width = 1500, height=400*ceiling(length(g.list)/3), res=180)
FeaturePlot(seur,features = g.list, pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
dev.off()

saveRDS(seur, "clusteredSign.rds")


#Integrate with RNA! ------
seur_rna = readRDS("mito210c1_d28_111520_final.rds")
transfer.anchors <- FindTransferAnchors(reference = seur_rna, query = seur, reduction = 'cca')
saveRDS(transfer.anchors, "transfer.anchors.rds")

predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = seur_rna$FinalName,weight.reduction = seur[['lsi']],dims = 2:30)

seur <- AddMetaData(object = seur, metadata = predicted.labels)

png("predicted.celltypes.png",res=100, width=600, height=400)
DimPlot(seur,group.by = 'predicted.id',label = TRUE,repel = TRUE) + ggtitle('Predicted Cell Types')
dev.off()

saveRDS(seur, "clusteredSign.rds")


#Differentially accessible peaks -----
DefaultAssay(seur) <- 'peaks'

da_peaks <- FindAllMarkers(seur, min.pct = 0.2, test.use = 'LR', latent.vars = 'peak_region_fragments')
da_peaks$peak = da_peaks$gene
da_peaks$gene = NULL
closest_genes = ClosestFeature(seur, regions = rownames(da_peaks))
da_peaks$closest_gene = closest_genes$gene_name
da_peaks$closest_gene_distance = closest_genes$distance
write_xlsx(da_peaks, path="DA_Peaks_LR.xlsx")

#Assign final cell types -----
#to the majority predicted ID per cluster, except where Ana determines otherwise from above excel
t = table(seur$seurat_clusters, seur$predicted.id)
calls = colnames(t)[apply(t, which.max, MARGIN = 1)]
seur$CellType = NA
for (i in 0:(length(calls)-1)) { seur$CellType[seur$seurat_clusters==i] = calls[i+1] }
seur$CellType[seur$seurat_clusters==8] = "IP"
Idents(seur) = "CellType"
saveRDS(seur, "clusteredSign.rds")

cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown')
cols = cols[match(levels(seur),levels)]

tiff("CellTypes.colors.1mATAC.tiff", res=120, width=650, height=400)
DimPlot(object = seur, label = F, cols=cols) + NoAxes()
dev.off()

cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown')
cols = cols[match(levels(factor(seur$predicted.id)),levels)]
names(cols) = levels(factor(seur$predicted.id))


png("predicted.celltypes.colors.png",res=120, width=650, height=400)
DimPlot(seur,group.by = 'predicted.id',label = F, cols = cols) + ggtitle('Predicted Cell Types') + NoAxes()
dev.off()

scores = seur@meta.data[,startsWith(colnames(seur@meta.data), "prediction.score")]
scores = scores[,!colnames(scores) %in% c("prediction.score.max","prediction.score.Subcortical")] #subcortical (without neurons/prog) is left over from previous RNA cell typing
colnames(scores) = substr(colnames(scores), 18, nchar(colnames(scores)))
ha = HeatmapAnnotation(FinalCellType = seur$CellType, col= list(FinalCellType=cols), which="row")
pdf("prediction.scores.heatmap.pdf", width=5, height=6)
Heatmap(as.matrix(scores), column_title = "Prediction Scores", show_row_names = F,
        show_column_dend = F, right_annotation = ha, col = c("white", "darkmagenta"))
dev.off()
tiff("prediction.scores.heatmap.tiff", res=100, width=500, height=600)
Heatmap(as.matrix(scores), column_title = "Prediction Scores", show_row_names = F,
        show_column_dend = F, right_annotation = ha, col = c("white", "darkmagenta"))
dev.off()

#Coembed UMAPs ------
transfer.anchors = readRDS("transfer.anchors.rds")
seur = readRDS("clusteredSign.rds")
seur_rna = readRDS("mito210c1_d28_111520_final.rds")
# note that we restrict the imputation to variable genes from scRNA-seq, could use all
genes.use <- VariableFeatures(seur_rna)
refdata <- GetAssayData(seur_rna, assay = "RNA", slot = "data")[genes.use, ]
# imputation will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seur[["lsi"]])
seur[["RNA"]] <- imputation
coembed <- merge(x = seur_rna, y = seur)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$jointcelltype <- ifelse(!is.na(coembed$FinalName), coembed$FinalName, coembed$CellType)
coembed$tech <- ifelse(!is.na(coembed$FinalName), "RNA", "ATAC")
p1 <- DimPlot(coembed, group.by = "tech") + NoAxes()

cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown')
cols = cols[match(levels(factor(coembed$jointcelltype)),levels)]
names(cols) = levels(factor(coembed$jointcelltype))
p2 <- DimPlot(coembed, group.by = "jointcelltype", label = F, cols=cols) + NoAxes()

pdf("ATACandRNA.CoEmbedded.pdf", width=10, height=4)
p1+p2
dev.off()
saveRDS(coembed,"coembed.rds")


# Proportions ----
#Percentage of cells in each cluster per organoid
counts = as.matrix(table(seur$CellType,seur$orig.ident))
counts = t(t(counts)/colSums(counts))
counts = melt(counts,varnames = c('CellType','Org'))
counts$Org = paste0("Org",as.numeric(counts$Org))
tiff('SUV1m.ATAC.percent.Cells.in.CellTypes.tiff',res=125,width=750,height=800)
ggplot(counts, aes(fill=CellType, y=value, x=Org)) + 
  geom_bar(stat="identity", position="fill", show.legend=T)+
  #geom_text(aes(label=paste(round(value*100), "%")), position= position_fill(vjust = 0.5)) +
  scale_fill_manual(values = cols) +
  labs(y="Percentage of Cells", x = "") +
  theme_minimal(base_line_size = 0) + NoLegend() + NoAxes()
dev.off()

#Calculate mutual information
library(mpmi)
assigns = seur@meta.data[,c("CellType","orig.ident")]
real_mi = dmi(assigns)$mi[1,2]
randomized_mi = list()
for (i in 1:1000){
  randomClusts = sample(assigns$CellType)
  assigns$CellType = randomClusts
  randomized_mi[i] = as.numeric(dmi(assigns)$mi[1,2])
}
randomized_mi = as.numeric(randomized_mi)
z = (real_mi - mean(randomized_mi)) / sd(randomized_mi)
#print MI score and z-score
real_mi
z


# Motif analysis ----
library(JASPAR2016) #want to update to 2020 but it wont install
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

seur = readRDS("clusteredSign.rds")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(JASPAR2016,opts = list(species = "Homo sapiens"))

# add motif information
seur <- AddMotifs(seur, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
DefaultAssay(seur) = "peaks"
Idents(seur) = "CellType"

# get top differentially accessible peaks
da_peaks <- FindAllMarkers(seur, test.use = 'LR', latent.vars = 'peak_region_fragments')
saveRDS(da_peaks, "da_peaks_perCellType.rds")
dir.create("homer_motifs_1m")
for (i in levels(factor(da_peaks$cluster))) { 
  top.da.peak <- da_peaks[da_peaks$cluster==i & da_peaks$p_val_adj < 0.1 & da_peaks$avg_logFC>0, "gene" ]
  print(paste(i, length(top.da.peak)))
  chr = sapply(strsplit(top.da.peak, "-"), "[", 1)
  start = sapply(strsplit(top.da.peak, "-"), "[", 2)
  end = sapply(strsplit(top.da.peak, "-"), "[", 3)
  strand = rep("+",length(end))
  t = data.frame(row.names = top.da.peak, "chr"=chr, "start"=start, "end"=end, "strand"=strand)
  write.table(t, file = paste0("homer_motifs_1m/peaks.",i,".txt"), row.names = T,col.names = F, quote=F, sep="\t")
}

#homer:
#PATH=$PATH:/stanley/levin_dr_storage/akedaigle/homer/bin/
#for i in *.txt; do /stanley/levin_dr_storage/akedaigle/homer/bin/findMotifsGenome.pl $i hg38 ${i%.txt} -size 300 -mask; done



# Cicero -----
library(SeuratWrappers)
library(cicero) #HERE not installing and monocle had issues too, come back to this

seur = readRDS("clusteredSign.rds")

cds <- as.cell_data_set(seur)
cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(seur)
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(cicero, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns, "conns.rds")

ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seur) <- links

#Coverage plot ----
#need fragments.tsv files
seur = readRDS("clusteredSign.rds")
Idents(seur) = "CellType"
cols = c('#8dd3c7','#bebada', '#fb9a99', '#80b1d3', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown')
cols = cols[match(levels(seur),levels)]


cov = CoveragePlot(seur,region = "chr3-62355000-62410000",links = F)
link = LinkPlot(seur, region = "chr3-62355000-62410000", min.cutoff = 0.1)
CombineTracks(plotlist = list(cov, link)) &
  scale_fill_manual(values=cols)
ggsave("covFEZF2.pdf")