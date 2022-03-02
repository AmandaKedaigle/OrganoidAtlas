library(Seurat)
library(reshape2)
library(gplots)
library(ggplot2)
library(viridis)
library(dplyr)
source("~/kwanho/src/cNMF.R")
source("~/kwanho/src/seurat_tools.R")

seur <- readRDS("/stanley/levin_dr/akedaigle/AtlasAna/CPNs/FromHodge_CPNBerg_ADULT.rds")
seur$cluster = droplevels(seur$cluster)

cpn_cols = c("#3c0d03", "#8d1c06", "#e67424", "#ed9b49", "#f5c34d")
names(cpn_cols) = levels(seur$cluster)
pdf("adult_CPN_subtypes.pdf")
DimPlot(seur, group.by='cluster', cols=cpn_cols)
dev.off()

print("Run cNMF!")
res <- Repro_NMF(seur , k=5)

print("Svaing!")
saveRDS(res, "res_cNMF.rds")

print("cNMF finished!")

#res = readRDS("res_cNMF.rds")

k=5
topn = 100

H = res$H
colnames(H)=sub("^","cNMF_",1:k)
rownames(H)=rownames(seur@meta.data)
seur[["cnmf"]] <- CreateDimReducObject(embeddings = H, key = "cNMF_", assay = DefaultAssay(seur))

W = res$W
colnames(W)=sub("^", "cNMF_", 1:k)
rownames(W) <- rownames(seur)
seur@reductions$cnmf@feature.loadings = W


# feature plot cell loading
seur <- AddMetaData(seur, as.data.frame.matrix(H))
plot_feature2(seur, features=colnames(H), nc=5, size=5, filename="feature_plot_cNMF_cell_loading.png", cols=magma(50), dev='png')

# heatmap cell loading
meta.df = seur@meta.data[, c('cluster', colnames(H))]
df = as.data.frame.matrix(meta.df)
df$cluster = factor(df$cluster, levels=levels(seur$cluster))
df = df[order(df$cluster), ]
my_sample_col <- data.frame(Cluster = rep(names(table(df$cluster)), times=table(df$cluster)))
rownames(my_sample_col) <- rownames(df)
my_sample_col[,1] = factor(my_sample_col[,1], levels=levels(seur$cluster))
df = df[,-1]
df = t(df)
library(pheatmap)
anno_cols = list(Cluster=cpn_cols)
pdf("heatmap_cNMF_cell_loading.pdf", height=10, width=10)
pheatmap(df, annotation_col=my_sample_col, cluster_rows=F, cluster_cols=F, show_colnames=F, scale='column', annotation_colors = anno_cols)
dev.off()

# violin cell loading
df = as.data.frame.matrix(meta.df)
scaled.df = cbind.data.frame(df[,1], scale(df[,2:ncol(df)]))
colnames(scaled.df)[1] = 'cluster'
mat = melt(df)
p <- ggplot(mat, aes(x=cluster, y=value, fill=cluster)) +
        geom_violin(trim=F, scale='width') +
        stat_summary(aes(y = value), fun='median', geom='crossbar') +
        facet_grid(cols=vars(variable)) +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='top')
pdf("violin_cNMF_cell_loading.pdf", height=8, width=12)
p + scale_fill_manual(values=cpn_cols)
dev.off()


top.genes = list()
for (i in 1:k) {
nam.col = paste0("cNMF_", i)
colnames(W)[i] = "col"
as.data.frame(W) %>% top_n(n=topn, wt=col) -> genes
print(quantile(genes$col))
top.genes[[nam.col]] = rownames(genes)
colnames(W)[i] = nam.col
}
names(top.genes) <- paste0(names(top.genes), "_top", topn, "_genes")

# save top genes in each module
saveRDS(top.genes, paste0("cNMF_top", topn, "_modules.rds"))
write.table(plyr::ldply(top.genes, rbind), paste0("table_cNMF_top", topn, "_genes.tsv"), sep='\t', quote=F, row.names=T, col.names=F)

outdir = paste0("top_", topn)
if (!dir.exists(outdir)) {
dir.create(outdir)
}
setwd(outdir)

source("~/kwanho/src/seurat_tools.R")
seur <- MyModuleScore(seur, gene.list=top.genes, save=T, filename="metadata_cNMF_module_scores.rds")

# feature plot module score
plot_feature2(seur, features=names(top.genes), title="cNMF Modules", size=3, filename="featureplot_gene_modules.pdf")

# violin plot module score
pl = list()
for (x in names(top.genes)) {
pl[[x]] <- VlnPlot(seur, features=x, pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
}
pdf("violin_gene_modules.pdf", height=20, width=20)
patchwork::wrap_plots(pl, ncol=3)
dev.off()

dir.create("individual_gene_plots")

# individual gene plots
ngenes = 20
for (nam in names(top.genes)) {
print(nam)
gtp = top.genes[[nam]]
gtp = intersect(gtp, rownames(seur))
gtp = gtp[1:ngenes]
pdf(paste0("violin_stacked_", gsub(topn, ngenes, nam), ".pdf"), height=25, width=10)
print(myStackedVlnPlot(seur, gtp))
dev.off()
FeaturePlotSingleLegend(seur, features=gtp, title=nam, size=3, filename=paste0("individual_gene_plots/featureplot_", nam, ".pdf"))
}

