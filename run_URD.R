library(Seurat)
library(URD)
library(useful)


args = commandArgs(trailingOnly=T)
knn = as.numeric(args[1])
sigma = as.numeric(args[2])
rc = as.numeric(args[3])  # 0: all aRG, 1: 23days aRG as root
root.cells = readRDS("root_cells_all_aRG.rds")
rc.nam = "all_aRG"
if (rc == 1) {
print("Using 23days aRG as the root cells")
root.cells = readRDS("root_cells_23d_aRG.rds")
rc.nam = "23days_aRG"
} else {
print("Using all aRG as the root cells")
}

source("tip_clusters.R")

print("Input arguments:")
print(paste0("knn=", knn))
print(paste0("sigma=", sigma))
print(paste0("root.cells=", rc.nam))


if (!file.exists("obj_urd_initial.rds")) {
print("loading Seurat!")
seur <- readRDS("seur_combined_final_downsampled.rds")
print(table(seur$FinalName))

seur$age_FinalName = paste(seur$age, seur$FinalName, sep='_')
cluster_ord = c()
for (tp in c('^23days', '^1m', '^1.5m','^2m','^3m','^4m','^5m','^6m')) {
cluster_ord = c(cluster_ord, grep(tp, names(table(seur$age_FinalName)), value=T))
}
seur$age_FinalName = factor(seur$age_FinalName, levels=cluster_ord)

# Create empty URD object
obj.urd <- new("URD")

# transfer metadata
obj.urd@logupx.data = seur@assays$RNA@data
obj.urd@meta <- seur@meta.data[,c('nCount_RNA','nFeature_RNA','percent.mito')]
obj.urd@group.ids <- seur@meta.data[,c('age_FinalName', 'FinalName', 'age', 'orig.ident')]

# transfer var.genes
obj.urd@var.genes <- seur@assays$RNA@var.features

# transfer UMAP embedding into tSNE slot of URD object
obj.urd@tsne.y <- as.data.frame(seur@reductions$umap@cell.embeddings)
colnames(obj.urd@tsne.y) = c("tSNE1", "tSNE2")

# transfer harmony results into PCA slot of URD object
obj.urd@pca.load <- as.data.frame(Loadings(seur, reduction = "pca"))
obj.urd@pca.scores <- as.data.frame(seur@reductions$pca@cell.embeddings)
obj.urd@pca.sdev <- as.numeric(apply(seur@reductions$pca@cell.embeddings, 2, stats::sd))
obj.urd@pca.sig <- c(rep(TRUE, 30), rep(FALSE, 20))

saveRDS(obj.urd, "obj_urd_initial.rds")

if (!dir.exists("figures")) {
dir.create("figures")
}

pdf("figures/urd_umap_check.pdf")
plotDim(obj.urd, "FinalName", plot.title="FinalName") + NoAxes()
plotDim(obj.urd, "age", plot.title = "age") + NoAxes()
dev.off()

rm(list=c('seur'))

} else {
print("initial URD object found!")
obj.urd <- readRDS("obj_urd_initial.rds")
}

print("DONE URD object creation")


###########################################

if (!dir.exists(rc.nam)) {
dir.create(rc.nam)
}
setwd(rc.nam)

###########################################
#obj.urd = readRDS("obj_urd_initial.rds")

if (!dir.exists("dm")) {
dir.create("dm")
}

# Calc diffusion map
dm.outfile = paste0("dm/dm_", knn, "_", sigma, ".rds")
if (!file.exists(dm.outfile)) {
print("Compute diffusion map!")
obj.urd <- calcDM(obj.urd, knn=knn, sigma.use=sigma)  # sigma.use
saveRDS(obj.urd@dm, dm.outfile)

pdf(paste0("dm/diffusion_map_knn=",knn,"_sigma=",sigma,".pdf"))
print(plotDimArray(obj.urd, reduction.use = "dm", dims.to.plot = 1:18, outer.title = paste0("Diffusion Map - knn=", knn, ", sigma=", sigma), label="FinalName", plot.title="", legend=F))
print(plotDim(obj.urd, "FinalName", transitions.plot = 10000, plot.title="Clusters with transitions"))
dev.off()
} else {
print("Computed diffusion map found!")
obj.urd@dm <- readRDS(dm.outfile)
}

if (!dir.exists("floods")) {
dir.create("floods")
}

flood.outfile = paste0("floods/flood_dm_", knn, "_", sigma, ".rds")
if (!file.exists(flood.outfile)) {
print("Do the flood simulations!")
flood.res = floodPseudotime(obj.urd, root.cells=root.cells, verbose=T)
saveRDS(flood.res, flood.outfile)
} else {
print("Flood simulation results found!")
flood.res <- readRDS(flood.outfile)
}

if (!dir.exists("pseudotime")) {
dir.create("pseudotime")
}

ptime.outfile = paste0("pseudotime/urd_obj_dm_", knn, "_", sigma, ".rds")
if (!file.exists(ptime.outfile)) {
print("Process simulations to pseudotime!")
obj.urd = floodPseudotimeProcess(obj.urd, flood.res, floods.name='pseudotime')
saveRDS(obj.urd, ptime.outfile)

pdf(paste0("pseudotime/pseudotime_stability_dm_", knn, "_", sigma, ".pdf"))
pseudotimePlotStabilityOverall(obj.urd)
dev.off()
} else {
print("Object with pseudotime found!")
obj.urd <- readRDS(ptime.outfile)
}

print("DONE pseudotime")


#########################################
# build trees

library(dplyr)
library(stringr)
library(cowplot)
library(grid)

walks.to.do = 10000
ocf = 40  # optimal.cells.forward
mcb = 80  # max.cells.back

#source("tip_clusters.R")

# Vis tips
print("Applying tips chosen!")
print(tipClusters)

Tips = list()
tipCells = c()
for (tc in tipClusters) {
curTipCells = cellsInCluster(obj.urd, "FinalName", tc)
tipCells = c(tipCells, curTipCells)
Tips[[tc]] <- curTipCells
}

obj.urd@group.ids[tipCells, "tip.clusters"] <- obj.urd@group.ids[tipCells, "FinalName"]
print(table(as.factor(obj.urd@group.ids$tip.clusters)))
tipName = names(table(as.factor(obj.urd@group.ids$tip.clusters)))
obj.urd@group.ids$tip.clusters <- as.factor(as.numeric(as.factor(obj.urd@group.ids$tip.clusters)))
print(table(obj.urd@group.ids$tip.clusters))
tipNum = names(table(obj.urd@group.ids$tip.clusters))

print("tips to be used:")
print(tipNum)
print(tipName)

outdir = paste0("dm-", knn, "-", sigma, "/tm-", ocf, "-", mcb)
if (!file.exists(paste0("walks/", outdir, "/obj_urd_walk_complete.rds"))) {

if (!dir.exists("vis_tips_on_DCs")) {
dir.create("vis_tips_on_DCs")
}

print("Visualizing tips on DCs!")
pdf(paste0("vis_tips_on_DCs/plot_tips_on_DCs_", knn, "_", sigma, ".pdf"))
print(plotDimArray(obj.urd, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Tip Clusters", label="tip.clusters", plot.title="", legend=F))
dev.off()


# biased random walk
if (!dir.exists("walks")) {
dir.create("walks")
}
setwd("walks")
dir.create(outdir, recursive=T)
setwd(outdir)
print(getwd())

if (!file.exists("diffusion_logistic.rds")) {
print("Find logistic function params!")
pdf("random_walk_logistics.pdf")
diffusion.logistic <- pseudotimeDetermineLogistic(obj.urd, "pseudotime", optimal.cells.forward=ocf,
	max.cells.back=mcb, pseudotime.direction="<", do.plot=T, print.values=T)
dev.off()
saveRDS(diffusion.logistic, "diffusion_logistic.rds")
} else {
print("Reading logistic function params!")
diffusion.logistic <- readRDS("diffusion_logistic.rds")
}

if (!file.exists("biased_transition_matrix.rds")) {
print("Compute biased TM!")
biased.tm <- pseudotimeWeightTransitionMatrix(obj.urd, pseudotime='pseudotime', 
	logistic.params=diffusion.logistic, pseudotime.direction="<")
saveRDS(biased.tm, "biased_transition_matrix.rds")
} else {
print("Reading biased TM!")
biased.tm <- readRDS("biased_transition_matrix.rds")
}

print("Biased TM dimension:")
print(dim(biased.tm))

obj.urd <- urdSubset(obj.urd, cells.keep=rownames(biased.tm))

if (!file.exists("random_walks.rds")) {
print("Simulate random walk!")
walks = simulateRandomWalksFromTips(obj.urd, tip.group.id="tip.clusters", root.cells=root.cells, 
	transition.matrix = biased.tm, n.per.tip=walks.to.do, root.visits=1, 
	verbose=T, max.steps=5000)
saveRDS(walks, "random_walks.rds")
} else {
print("Reading random walk!")
walks = readRDS("random_walks.rds")
}

print("Processing random walks!")
obj.urd <- processRandomWalksFromTips(obj.urd, walks, verbose = T)
saveRDS(obj.urd, "obj_urd_walk_complete.rds")

} else {
setwd(paste0("walks/", outdir))
print("Reading URD object with random walk complete")
obj.urd <- readRDS("obj_urd_walk_complete.rds")
}

# build tree
print("Building trees!")
setwd("../../../")

if (!dir.exists("trees")) {
dir.create("trees")
}
setwd("trees")

dir.create(outdir, recursive=T)
setwd(outdir)
print(getwd())

markers = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/pseudotime_URD/merge_all/subset/all_aRG_as_root/DEG_FinalName/markers_MAST_merged_FinalName.rds")
markers %>% group_by(cluster) %>% top_n(n=6, wt=avg_logFC) -> top.genes

myTree_unprocessed <- loadTipCells(obj.urd, "tip.clusters")

divergence.methods <- c("preference")
p.threshs = c(0.05)
cells.per.pseudotime.bins <- c(80, 50, 100)  # 80
bins.per.pseudotime.windows <- c(5, 3, 8)  # 5
visit.thresholds <- c(0.7)  # 0.7
minimum.visitss <- c(10)  # 10

for (divergence.method in divergence.methods ){
  for (visit.threshold in visit.thresholds ){
    for (minimum.visits in minimum.visitss){
      for (cells.per.pseudotime.bin in cells.per.pseudotime.bins ){
        for (bins.per.pseudotime.window in bins.per.pseudotime.windows ){
          for (p.thresh in p.threshs ){

            myTree = myTree_unprocessed

            paramSet = paste0("dm", divergence.method, "_vt", visit.threshold, "_mv", minimum.visits, "_cppb", cells.per.pseudotime.bin, "_bppw", bins.per.pseudotime.window, "_pt", p.thresh)
            print(paramSet)

            dir.create(paramSet)
            setwd(paramSet)

            treeName = paste0("obj_urd_tree_", paramSet, ".rds")

            tryCatch(
                {
                  myTree <- buildTree(myTree,
                        pseudotime = "pseudotime",
                        tips.use=tipNum,
                        divergence.method = divergence.method,
                        visit.threshold = visit.threshold,
                        minimum.visits = minimum.visits,
                        cells.per.pseudotime.bin = cells.per.pseudotime.bin,
                        bins.per.pseudotime.window = bins.per.pseudotime.window,
                        save.all.breakpoint.info = T,
                        save.breakpoint.plots = "breakpoint_plot",
                        min.cells.per.segment = 1,  # 1
                        p.thresh = p.thresh)

                  myTree <- nameSegments(myTree, segments=tipNum, segment.names = tipName, short.names = tipName)
                  segs <- myTree@tree$segments

                  for (tip in tipClusters) {
                    value <- list()
                    for (i in segs) {
                      prop = 100*sum(Tips[[tip]] %in% myTree@tree$cells.in.segment[[i]], na.rm = TRUE)/length(Tips[[tip]])
                      if(prop>0){
                        line = paste0(tip,"   ",i,"   ",prop)
                        #print(line)
                        cat(line, file="tip_prop.txt", sep='\n', append=T)
                      }
                    }
                  }


                  print("saving tree!")

                  saveRDS(myTree, treeName)

                  print("plotting!")

                  pdf("tree_and_segments.pdf")
                  print(plotTree(myTree, "age"))
                  print(plotTree(myTree, "FinalName"))
                  p2 <- plotTree(myTree, "segment",label.segments=T)
                  p3 <- plotDim(myTree, "segment",label.clusters = T, legend = F, plot.title="tree segment")
                  print(plot_grid(p2,p3, ncol = 2))
                  dev.off()

		  for (cur.clust in c('CPN','CFuPN','Astroglia','Glial precursors')) {
		    print(cur.clust)
		    cur.genes = top.genes %>% filter(cluster==cur.clust) %>% pull(gene)
		    cur.clust = gsub(" |/", "_", cur.clust)
                    
		    plist = list()
		    for (cg in cur.genes) {
		      plist[[cg]] = plotTree(myTree, cg, tree.size = .5, cell.alpha = 0.5, cell.size = 0.25)
		    }

		    pdf(paste0("tree_gene_", cur.clust, ".pdf"), height=8, width=9)
                    print(gridExtra::grid.arrange(grobs=plist, nrow=2,
		    			top = textGrob(cur.clust, gp=gpar(fontsize=20,font=3))))
                    dev.off()
		  }

                },
                error = function(cond) {
                  message(paste("paramSet didn't work:", paramSet))
                  message("Here's the original error message:")
                  message(cond)
                },
                finally={
                setwd('../')
                }
            )
          }
        }
      }
    }
  }
}

print("DONE")

