##complex heatmap from DESeq2 excels

library(ComplexHeatmap)
library(readxl)
library(viridis)

#reading in files
dataset_names = c("name1", "name2", "name3")
dataset_filenames = c("name1.xlsx", "name2.xlsx", "name3.xlsx")
dataset_files = lapply(dataset_filenames, read_xlsx)

#Setting thresholds for what to include in heatmap -- this can be easily changed
padj = 0.0015
absFC = 1.5

#Get list of all genes that pass those thresholds in any dataset
genes = unique(unlist(lapply(dataset_files, function(x) x[x$padj<padj & abs(x$log2FoldChange)>absFC,"gene"])))

#Get genes and logFC
genes = genes[!is.na(genes)]
genesDF = data.frame(row.names = genes)
for (d in 1:length(dataset_files)) {
  dat = dataset_files[[d]]
  dat = dat[match(genes, dat$gene),]
  genesDF$placeholderName = dat$log2FoldChange
  names(genesDF)[names(genesDF)=="placeholderName"] = dataset_names[d]
}
genesDF[is.na(genesDF)] = 0

#genes you want to highlight - pick the interesting genes to you!
highlight = c("genes", "of", "interest")

highlight = highlight[highlight %in% genes]

##color

col = magma(256, direction = -1)
col = viridis(256, direction = -1)

#make heatmap

heatmap <- Heatmap(genesDF,
                   column_names_rot = 45, 
                   col = col,
                   heatmap_legend_param = list(title = "Avg LogFC")) +
                   rowAnnotation(link = anno_mark(at = match(highlight, rownames(genesDF)),
                                 labels = highlight,
                                 labels_gp = gpar(fontsize = 11)))



pdf('filename.pdf', width = 8, height = 16)
print(heatmap)
dev.off()

tiff('filename.tiff', res = 300, width = 1400, height=2500) 
print(heatmap)
dev.off()

