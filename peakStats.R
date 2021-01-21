library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(reshape2)
library(data.table)

setwd("~/Documents/scATAC/mergeWT/")

#Just cell-type specific % in each age seperately ----
seur1 = readRDS("../SUV_Mito210_1m/wt/Signac/clusteredSign.rds")
seur3 = readRDS("../SUV_Mito210_3m/wt/Signac/clusteredSign.rds")
seur6 = readRDS("../SUV_Mito210_6m/wt/Signac/clusteredSign.rds")
da_pks1 = readRDS("../SUV_Mito210_1m/wt/Signac/da_peaks_perCellType.rds")
da_pks3 = readRDS("../SUV_Mito210_3m/wt/Signac/da_peaks_perCellType.rds")
da_pks6 = readRDS("../SUV_Mito210_6m/wt/Signac/da_peaks_perCellType.rds")

length(unique(da_pks1$gene)) / nrow(seur1)
length(unique(da_pks3$gene)) / nrow(seur3)
length(unique(da_pks6$gene)) / nrow(seur6)

specificity = sapply(rownames(seur6), function (x) {
  if (nrow(da_pks6[da_pks6$gene==x,]) == 0) "Nonspecific"
  else if (nrow(da_pks6[da_pks6$gene==x,]) == 1) as.character(da_pks6[da_pks6$gene==x, "cluster"])
  else "Up in >1 CellType"})

t = table(specificity)
tm = melt(t)
tm$specificity = relevel(factor(tm$specificity), ref = "Nonspecific")
levels = c("Nonspecific", 'aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown', "Up in >1 CellType")
cols = c("gray70",'#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9', "black")
cols = cols[match(levels(tm$specificity),levels)]
ggplot(tm, aes("", value, fill=specificity)) + 
  geom_bar(stat="identity") + coord_polar("y") +
  theme_minimal() +
  scale_fill_manual(values=cols) +
  labs(x="", y="", fill="Peak Specificity")
ggsave("PeakSpecificity-6m.pie.pdf")


#Compare to enhancer and promoter regions
# read in peak sets
peaks.1<- read.table("../SUV_Mito210_1m/wt/Signac/peaks.bed",col.names = c("chr", "start", "end"))
peaks.3<- read.table("../SUV_Mito210_3m/wt/Signac/peaks.bed",col.names = c("chr", "start", "end"))
peaks.6<- read.table("../SUV_Mito210_6m/wt/Signac/peaks.bed",col.names = c("chr", "start", "end"))
gr.1 <- makeGRangesFromDataFrame(peaks.1)
gr.3 <- makeGRangesFromDataFrame(peaks.3)
gr.6 <- makeGRangesFromDataFrame(peaks.6)

#Enhancers from genehancer
genehancer = read.table("GeneHancer_version_4-4.txt", sep="\t", header = T)
gh = genehancer[,c(1,4,5)]
colnames(gh) = c("chr", "start", "end")
gr.gh = makeGRangesFromDataFrame(gh)

ol1 = overlapsAny(gr.1, gr.gh)
ol3 = overlapsAny(gr.3, gr.gh)
ol6 = overlapsAny(gr.6, gr.gh)

#Promoters & gene bodies as calculated for gene activity by Signac

CollapseToLongestTranscript <- function(ranges) {  #Function from Signac code!
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
    ]
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
transcripts <- CollapseToLongestTranscript(ranges = annotations)
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"]
# extend to include promoters
transcripts <- Extend(transcripts,upstream = 2000,downstream = 0) #default values

pro1 = overlapsAny(gr.1, transcripts)
pro3 = overlapsAny(gr.3, transcripts)
pro6 = overlapsAny(gr.6, transcripts)

fracOverlaps = data.frame(time = c("1m", "3m", "6m"), fracEnh = c(sum(ol1)/ length(ol1), sum(ol3)/ length(ol3), sum(ol6)/ length(ol6)),
                          fracPro = c(sum(pro1)/length(pro1), sum(pro3)/length(pro3), sum(pro6)/length(pro6)))
f = melt(fracOverlaps)
ggplot(f, aes(time, value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() + labs(x="", y="Fraction of Peaks")+
  scale_fill_manual(values=c("#bb5566", "#004488"), name="Overlap", labels=c("Enhancers", "Promoters & Gene Bodies"))

ggsave("peakPromotersAndEnhancers.pdf")