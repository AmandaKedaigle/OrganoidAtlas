.libPaths(.libPaths()[2])

library(Seurat)
library(ggplot2)
library(stringr)
library(reshape2)


setwd("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/classifier")

# Dot plot
#out.name = "share-seq_RNA_top20_genes"
out.name = "share-seq_RNA_top20_TFs"
#sorted.genes = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/classifier/all_genes_in_bp_seg_clean.rds")
sorted.genes = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/classifier/all_TFs_in_bp_seg.rds")
top20 = lapply(sorted.genes, function(x) {x[1:20]})
seg.ord = c('bp1_segment1','bp1_segment3','bp1_segment2')
top20 = top20[seg.ord]
arr = unlist(top20, use.names=F)
gl = lapply(top20, function(x) {names(x)})
names(arr) = unlist(gl, use.names=F)
arr = arr[!duplicated(names(arr))]

seur <- readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/data/seur_RNA_matching_barcodes.rds")
tree <- readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/res_URD/trees/dm-100-20/tm-40-80/dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05/obj_urd_tree_dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05.rds")
meta.df = cbind(tree@group.ids, tree@pseudotime)
seur = AddMetaData(seur, meta.df)
Idents(seur) <- "segment"

# only the cells that are visited during random walk are included in the tree
sseur <- subset(seur, cells=setdiff(colnames(seur), names(seur$segment[which(is.na(seur$segment))])))
VariableFeatures(sseur) <- names(arr)
sseur <- ScaleData(sseur)
# average exp
avg.exp = AverageExpression(sseur, features=names(arr), slot='scale.data')[[1]]
avg.exp = avg.exp[, as.character(c(3,1,2))]
colnames(avg.exp) <- paste0("Segment", colnames(avg.exp))
avg.exp = as.data.frame(avg.exp)

df.list = list()
for (i in c(1)) {
  bp = paste0("bp", i)
  seg.idx = seq(1,3) + 3*(i-1)
  segs = gsub("s", "S", str_split(seg.ord[seg.idx], '_', simplify=T)[,2])
  cur.gl = gl[grep(bp, names(gl), value=T)]
  df = as.data.frame(matrix(, nrow=0, ncol=3, dimnames=list(NULL,segs)))
  for (g in unlist(cur.gl, use.names=F)) {
    dup = length(grep(g, rownames(df)))
    if (dup==0) {
      df = rbind(df, avg.exp[which(rownames(avg.exp)==g), segs])
    } else {
      new.row = avg.exp[which(rownames(avg.exp)==g), segs]
      rownames(new.row) = paste0(rownames(new.row), ".", dup)
      df = rbind(df, new.row)
    }
  }
  df.list[[bp]] = df
}


# add importance score
#all.genes = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/pseudotime_URD/merge_all/final/23days_aRG/downstream/sigma30/branchpoint_DE/classifier/all_genes_in_bp_seg_clean.rds")
#all.genes = all.genes[which(!(names(all.genes) %in% c('bp2_segment6','bp3_segment5')))]
all.genes = sorted.genes[seg.ord]

x.cols = c('darkblue', 'black', 'darkred')
y.cols = rev(c(rep('black', 20), rep('darkblue', 20), rep('darkred', 20)))

plist = list()

for (bp in paste0("bp", c(1))) {
  df = df.list[[bp]]
  print(dim(df))
  plot.data = melt(as.matrix(df), varnames=c('Gene', 'Segment'), value.name="Exp")
  plot.data[["Score"]] = rep(NA, nrow(plot.data))
  print(dim(plot.data))

  cur.gl = gl[grep(bp, names(gl), value=T)]
  gene.set = unlist(cur.gl, use.names=F)

  temp = all.genes[grep(bp, names(all.genes), value=T)]
  sc = c()
  for (li in temp) {
    sc = c(sc, li[gene.set])
  }
  plot.data$Score = sc
  plot.data$sqrt.Score = sqrt(plot.data$Score)
  saveRDS(plot.data, paste0("dotplot_", out.name, "_", bp, "_data.rds"))
  plot.data$Gene = factor(plot.data$Gene, levels=c(levels(plot.data$Gene)[c(21:40, 1:20, 41:60)]))
  plot.data$sqrt.Score = ifelse(plot.data$sqrt.Score==0, NA, plot.data$sqrt.Score)

  plist[[bp]] = ggplot(plot.data) + 
  	geom_point(aes(x=Segment, y=Gene, size=sqrt.Score, colour=Exp)) +
        theme_classic() + 
	theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, colour=x.cols),
	  axis.text.y=element_text(colour=y.cols),
	  axis.text=element_text(size=10),
	  axis.title.x=element_blank(),
	  axis.title.y=element_blank()) +
	scale_y_discrete(limits = rev(levels(plot.data$Gene))) +
        #scale_colour_gradient(low='#19D719', high='#F02D7D') +
	scale_colour_viridis_c() +
	scale_x_discrete(position='top') +
	geom_hline(yintercept=c(20.5, 40.5)) +
	ggtitle(gsub("bp", "Branching Point ", bp))
}

pdf(paste0("dotplot_split_", out.name, "_by_importance_0_removed.pdf"), height=10, width=17/6)
#gridExtra::grid.arrange(grobs=plist, ncol=length(plist), as.table=F)
patchwork::wrap_plots(plist, nrow=1)
dev.off()



redraw_bp_dotplot = function(out.name="top20_TFs", width=10) {
  require(ggplot2)
  require(patchwork)
  x.cols = c('darkblue', 'black', 'darkred')
  y.cols = rev(c(rep('black', 20), rep('darkblue', 20), rep('darkred', 20)))

  plist = list()

  for (bp in paste0("bp", c(1,2,3))) {
    plot.data = readRDS(paste0("dotplot_", out.name, "_", bp, "_data.rds"))
    plot.data$Gene = factor(plot.data$Gene, levels=c(levels(plot.data$Gene)[c(21:40, 1:20, 41:60)]))

    plist[[bp]] = ggplot(plot.data) + geom_point(aes(x=Segment, y=Gene, size=sqrt.Score, colour=Exp)) +
        theme_classic() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, colour=x.cols),
          axis.text.y=element_text(colour=y.cols),
          axis.text=element_text(size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
        scale_y_discrete(limits = rev(levels(plot.data$Gene))) +
        #scale_colour_gradient(low='#19D719', high='#F02D7D') +
        scale_colour_viridis_c() +
        scale_x_discrete(position='top') +
        geom_hline(yintercept=c(20.5, 40.5)) +
        ggtitle(gsub("bp", "Branching Point ", bp))
  }

  pdf(paste0("dotplot_split_", out.name, "_by_importance_reordered.pdf"), height=10, width=width)
  print(wrap_plots(plist, nrow=1))
  dev.off()

}

