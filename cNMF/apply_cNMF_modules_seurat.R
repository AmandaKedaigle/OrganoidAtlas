library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
source("~/kwanho/src/seurat_tools.R")


gene.list <- readRDS("cNMF_top100_modules.rds")
print(lapply(gene.list, function(x) grep("LAMP5|LTK|LINC00507|FREM3|GLP2R|RORB|CARM1P1|COL22A1", x, value=T)))
#mod.keep = c(1,2,3,4,5,12,15)
#gene.list = gene.list[mod.keep]


datadir = "/stanley/levin_dr/akedaigle/AtlasAna/CPNs"

obj.list = list.files(datadir, pattern="*.rds", full.names=T)
obj.list = setdiff(obj.list, grep("FromHodge_CPNBerg_ADULT.rds", obj.list, value=T))
obj.names = gsub(".rds", "", basename(obj.list))

print(obj.names)

outdir = "top_100/exp_on_other_datasets"
if (!dir.exists(outdir)) {
dir.create(outdir, recursive=T)
}

for (i in 1:length(obj.list)) {
nam = obj.names[i]
print(nam)

res.seur.file = paste0(outdir, "/seur_", nam, ".rds")
if (!file.exists(res.seur.file)) {
seur <- readRDS(obj.list[i])

print(table(Idents(seur)))

matching.genes = lapply(gene.list, function(x) intersect(x, rownames(seur)))
names(matching.genes) <- paste0(names(matching.genes), '_', sapply(matching.genes, length), "_match")
print(names(matching.genes))
df=plyr::ldply(matching.genes, rbind)
write.table(df, paste0(outdir, "/matching_genes_", nam, ".tsv"), sep='\t', row.names=T, col.names=F, na="")

# module score file name
msfn = paste0(outdir, "/metadata_", nam, "_cNMF_module_score.rds")
if (file.exists(msfn)) {
seur <- AddMetaData(seur, readRDS(msfn))
} else {
seur <- MyModuleScore(seur, gene.list=matching.genes, save=T, filename=msfn)
}

plot_feature2(seur, features=names(matching.genes), title=paste0("In vivo cNMF module exp on ", nam), filename=paste0(outdir, "/feature_plot_", nam, "_cNMF_module_scores.pdf"), size=4)

# process using module genes as variable genes and identify clusters
var.genes = unique(unname(unlist(matching.genes)))
VariableFeatures(seur) <- var.genes

seur <- seur %>%
	ScaleData() %>%
	RunPCA()

pdf(paste0(outdir, "/pca_elbow_plot_", nam, ".pdf"))
print(ElbowPlot(seur))
dev.off()

seur <- seur %>%
	FindNeighbors(dims=1:10) %>%
	FindClusters(resolution=0.3) %>%
	RunUMAP(dims=1:10)

print(table(Idents(seur)))
saveRDS(seur, res.seur.file)
} else {
print("processed seur found!")
seur <- readRDS(res.seur.file)
matching.genes = lapply(gene.list, function(x) intersect(x, rownames(seur)))
names(matching.genes) <- paste0(names(matching.genes), '_', sapply(matching.genes, length), "_match")
}

glo.min=min(apply(seur@meta.data[,grep("cNMF_", colnames(seur@meta.data), value=T)], 2, min))
glo.max=max(apply(seur@meta.data[,grep("cNMF_", colnames(seur@meta.data), value=T)], 2, max))

cols = magma(50)
vpl = list()
for (fea in names(matching.genes)) {
vpl[[fea]] = VlnPlot(seur, features=fea, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend() + ylim(round(glo.min), round(glo.max))
}

# module score on processed seur
#pdf(paste0(outdir, "/module_score_plots_", nam, ".pdf"), height=8, width=8)
#print(DimPlot(seur, label=T) + ggtitle(paste0(nam, ": variable genes = CPN module genes")))
#plot_feature(seur, features=names(matching.genes), title=paste0("In vivo cNMF module exp on ", nam), cols=cols)
#print(patchwork::wrap_plots(vpl, ncol=3))
#dev.off()

png(paste0(outdir, "/umap_", nam, ".png"), height=5, width=5, unit='in', res=300)
print(DimPlot(seur, label=T) + ggtitle(paste0(nam, ": variable genes = CPN module genes")))
dev.off()

plot_feature2(seur, features=names(matching.genes), title=paste0("In vivo cNMF module exp on ", nam), cols=cols, nc=5, dev='png', size=5, filename=paste0(outdir, "/feature_plot_", nam, "_new_umap.png"), res=300)

png(paste0(outdir, "/violin_", nam, "_new_clusters.png"), height=5, width=20, unit='in', res=300)
print(patchwork::wrap_plots(vpl, ncol=5))
dev.off()
}


