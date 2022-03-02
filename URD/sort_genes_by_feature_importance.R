library(stringr)


workdir = "/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/classifier"
setwd(workdir)

## filter.heatmap.genes function
# Removes undesired (mitochondrial, ribosomal, tandem duplicated genes) from heatmaps for presentation purposes.
filter.heatmap.genes <- function(genes #(Character vector) genes to check
) {
  hb.genes <- grep("^HBB|HBA", ignore.case = T, genes, value = T)
  mt.genes <- grep("^MT-", ignore.case = T, genes, value = T)
  ribo.genes <- grep("^RPL|^RPS", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(hb.genes, mt.genes, ribo.genes)))
}


res.files = list.files(pattern=".csv$")

sorted.genes = list()

for (rf in res.files) {
fi.scores = read.table(rf, sep=',', header=F)

bp = str_split(rf, "_", simplify=T)[,1]
seg = gsub(".csv", "", str_split(rf, "_", simplify=T)[,2])
print(paste0(bp, '_', seg))

colnames(fi.scores) <- readRDS(paste0("../", bp, "_genes.rds"))

fi = colSums(fi.scores)/nrow(fi.scores)

fi = fi[filter.heatmap.genes(names(fi))]

# select all
fi = sort(fi, decreasing=T)
fi = fi[fi>0]
# select top 20 genes/TFs
#fi = fi[1:20]

# add to list
sorted.genes[[paste0(bp, "_", seg)]] = fi
}

saveRDS(sorted.genes, "all_genes_in_bp_seg_clean.rds")
gene.names = lapply(sorted.genes, function(x) names(x))
write.table(plyr::ldply(gene.names, rbind), "sorted_genes_by_feature_importance_per_bp_clean.tsv", sep='\t', quote=F, row.names=T, col.names=F, na="")

tf.ref = read.table("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/pseudotime_URD/merge_all/final/23days_aRG/downstream/dataCSV2748.csv", sep=',', header=T)
tfs = as.character(tf.ref$Name)

my.list = lapply(sorted.genes, function(x) x[intersect(names(x), tfs)])

saveRDS(my.list, "all_TFs_in_bp_seg.rds")
tf.names = lapply(my.list, function(x) names(x))
write.table(plyr::ldply(tf.names, rbind), "sorted_TFs_by_feature_importance_per_bp_clean.tsv", sep='\t', quote=F, row.names=T, col.names=F, na="")

