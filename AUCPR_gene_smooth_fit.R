library(URD)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)


##############
## Functions
##############

## threshold.tree.markers function
# Function to threshold markers from a markersAUCPRAlongTree test with additional criteria
threshold.tree.markers <- function(markers, #list of results from markersAUCPRAlongTree tests (gene.markers)
                                   tip, #which tip (or element of the list to pursue)
                                   global.fc = 0.1, #fold.change that gene must have along the trajectory pursued vs. rest of the data
                                   branch.fc = 0.4, #fold.change that gene must have (in best case) vs. the opposing branch at any branchpoint along the trajectory.
                                   aucpr.ratio.all = 1.03 #classifier score that gene must exhibit along trajectory test vs. rest of the data
) #Returns markers with only a subset of rows retained.
{
  m <- markers[[tip]][[1]]
  print(dim(m))
  # First off -- lose global FC < x
  bye.globalfc <- rownames(m)[m$expfc.all < global.fc]
  # Second -- get rid of branch FC < x
  bye.branchfc <- rownames(m)[m$expfc.maxBranch < branch.fc]
  # Third -- get rid of stuff essentially worse than random classification on
  # global level
  bye.badglobalaucpr <- rownames(m)[m$AUCPR.ratio.all < aucpr.ratio.all]
  bye.all <- unique(c(bye.globalfc, bye.branchfc, bye.badglobalaucpr))
  print(length(bye.all))
  m.return <- m[setdiff(rownames(m), bye.all), ]
  print(dim(m.return))
  return(m.return)
}




# function to determine timing
determine.timing_DDB <- function(s, #result from geneSmoothFit
                                 genes = rownames(s$mean.expression) #genes to order; default is all genes that were fit.
) #Returns s but with an additional list entry ($timing) of the order to plot genes
{
  s$timing <- as.data.frame(do.call("rbind", lapply(genes, function(g) {
    sv <- as.numeric(s$scaled.smooth[g, ])
    pt <- as.numeric(colnames(s$scaled.smooth))
    # Figure out baseline expression & threshold for finding peaks
    min.val <- max(min(sv), 0)
    peak.val <- ((1 - min.val)/2) + min.val
    exp.val <- ((1 - min.val)/5) + min.val
    # Run-length encoding of above/below the peak-threshold
    peak.rle <- rle(sv >= peak.val)
    peak.rle <- data.frame(lengths = peak.rle$lengths, values = peak.rle$values)
    peak.rle$end <- cumsum(peak.rle$lengths)
    peak.rle$start <- head(c(0, peak.rle$end) + 1, -1)
    # Run-length encoding of above/below the expressed-threshold
    exp.rle <- rle(sv >= exp.val)
    exp.rle <- data.frame(lengths = exp.rle$lengths, values = exp.rle$values)
    exp.rle$end <- cumsum(exp.rle$lengths)
    exp.rle$start <- head(c(0, exp.rle$end) + 1, -1)
    # Take top-two longest peak RLE & select later one. Find stretches that are
    # above peak value
    peak <- which(peak.rle$values)
    
    ## DDB - nov 2020: im adding a modifiation here.
    # Check if there is a peaak. If there is none, remove gene
    if (length(peak) == 0){
      v <- c(NA, NA, NA, NA, NA)
      names(v) <- c("pt.start", "pt.peak", "pt.end", "exp.start", "exp.end")
      return(v)
    } else {
      # Order by length and take 1 or 2 longest ones
      peak <- peak[order(peak.rle[peak, "lengths"], decreasing = T)][1:min(2, length(peak))]
      # Order by start and take latest one.
      peak <- peak[order(peak.rle[peak, "start"], decreasing = T)][1]
      # Identify the actual peak value within that stretch
      peak <- which.max(sv[peak.rle[peak, "start"]:peak.rle[peak, "end"]]) + peak.rle[peak,
                                                                                      "start"] - 1
      # Identify the start and stop of the expressed stretch that contains the peak
      exp.start <- exp.rle[which(exp.rle$end >= peak & exp.rle$start <= peak),
                           "start"]
      exp.end <- exp.rle[which(exp.rle$end >= peak & exp.rle$start <= peak), "end"]
      # Identify values of expression at start and stop
      smooth.start <- sv[exp.start]
      smooth.end <- sv[exp.end]
      # Convert to pseudotime?
      exp.start <- pt[exp.start]
      exp.end <- pt[exp.end]
      peak <- pt[peak]
      # Return a vector
      v <- c(exp.start, peak, exp.end, smooth.start, smooth.end)
      names(v) <- c("pt.start", "pt.peak", "pt.end", "exp.start", "exp.end")
      return(v)
    }})))
  rownames(s$timing) <- genes
  # Remove rows that contain NA - genes where no peak was found
  s$timing = s$timing[complete.cases(s$timing),]
  # Decide on ordering of genes
  s$gene.order <- rownames(s$timing)[order(s$timing$pt.peak, s$timing$pt.start,
                                           s$timing$pt.end, s$timing$exp.end, decreasing = c(F, F, F, T), method = "radix")]
  return(s)
}  


## filter.heatmap.genes function
# Removes undesired (mitochondrial, ribosomal, tandem duplicated genes) from heatmaps for presentation purposes.
filter.heatmap.genes <- function(genes #(Character vector) genes to check
) {
  hb.genes <- grep("^HBB|HBA", ignore.case = T, genes, value = T)
  mt.genes <- grep("^MT-", ignore.case = T, genes, value = T)
  ribo.genes <- grep("^RPL|^RPS", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(hb.genes, mt.genes, ribo.genes)))
}



###########
## Start!
###########
print("loading data!")
#runs = list.dirs("trees", recursive=F)

#treepath = file.path("trees/dm-100-30", "tm-40-80", "dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05")  # default param
#treefile = list.files(path=treepath, pattern="^obj")
#treename = file.path(treepath, treefile)
#print(treename)

tree = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/res_URD/trees/dm-100-20/tm-40-80/dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05/obj_urd_tree_dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05.rds")

source("../tip_clusters.R")  # loads tipClusters

## aucpr
if (!dir.exists("aucpr")) {
dir.create("aucpr")
}
setwd("aucpr")

if (!file.exists("list_aucpr_res.rds")) {

# Used in aucprTestAlongTree. Errors out if omitted
tree@meta$n.Genes = tree@meta$nFeature_RNA
tree@meta$n.Trans = tree@meta$nCount_RNA

res <- list()
for (tip in tipClusters) {
print(paste0(Sys.time(), ": ", tip))
markers = aucprTestAlongTree(tree, pseudotime="pseudotime", tips=tip, report.debug=T)
tipname <- gsub(" ", "_", tip)
saveRDS(markers, paste0("marker_genes_", tipname, ".rds"))
res[[tip]] <- markers
}
saveRDS(res, "list_aucpr_res.rds")
markers.to.use=res
} else {
print("gene expression cascade results found! skipping")
markers.to.use <- readRDS("list_aucpr_res.rds")
}

setwd('../')


tips=names(markers.to.use)
seg.tip = list(CFuPN=as.character(c(3,1)), CPN=as.character(c(3,2)))

if (!dir.exists("cascade")) {
dir.create("cascade")
}
setwd("cascade")

## filtering param
filt_fc = 0.1
filt_AUCPR = 1.05

## for saving:
cascades_spline = list()

for (tip in tips){
  print(tip)
  # Get markers from aucprTestAlongTree and filter genes based on fc and aupcr ratio
  markers <- threshold.tree.markers(markers.to.use,
                                    tip = tip,
                                    global.fc = filt_fc,
                                    aucpr.ratio.all = filt_AUCPR)
  print(dim(markers))
  # Calculate spline curves Using segments for each traj
  spline <- geneSmoothFit(tree, pseudotime = "pseudotime",
                          cells = cellsInCluster(tree,"segment", seg.tip[[tip]]),
                          genes = rownames(markers),
                          method = "spline",
                          moving.window = 5,
                          cells.per.window = 25, 
			  #pseudotime.per.window = 0.005,
                          spar = 0.5, verbose = T)
  saveRDS(spline, paste0("spline_fit_", tip, ".rds"))
  # save spline (most time consuming step)
  cascades_spline[[tip]] = spline
}

print("saving cascades")
saveRDS(cascades_spline, "list_cascades_spline.rds")



# determine timing
print("determine timing")
cascades_spline_timing = list()
cascades_order = list()
for (i in 1:length(cascades_spline)){
  spline = cascades_spline[[i]]
  tip = tips[i]
  # Calculate gene expression timing for ordering rows
  spline <- determine.timing_DDB(s = spline)
  order <- filter.heatmap.genes(spline$gene.order)
  markers <- markers.to.use[[i]][["diff.exp"]]
  # Output gene table
  table.save <- data.frame(gene = order, stringsAsFactors = F)
  table.save$rgc.AUCPR.ratio.all <- markers[table.save$gene, "AUCPR.ratio.all"]
  table.save$rgc.AUCPR.ratio.maxBranch <- markers[table.save$gene, "AUCPR.ratio.maxBranch"]
  table.save$rgc.exp.fc.all <- markers[table.save$gene, "expfc.all"]
  table.save$rgc.exp.fc.best <- markers[table.save$gene, "expfc.maxBranch"]
  write.csv(table.save, quote = F, file = paste0("data_timing_", tip, ".csv"))

  # save spline with timing
  #cascades_spline_timing[[i]]= spline
  cascades_spline_timing[[i]]= spline
  # save order
  cascades_order[[i]] = order
}
names(cascades_spline_timing) = names(markers.to.use)
names(cascades_order) = names(markers.to.use)
saveRDS(cascades_spline_timing, "cascades_spline_timing.rds")
saveRDS(cascades_order, "cascades_gene_order.rds")


## GENERATE HEATMAP: ALL GENES
# colors
library(scales)
cols <- scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd"))(seq(0,1,length.out = 50))
library(viridis)
cols <- viridis(50)

# Make sure any values <0 in the spline curves get set to 0 so that the heatmap scale doesn't get messed up.
for (i in 1:length(tips)){
  spline = cascades_spline[[i]]
  order = cascades_order[[i]]
  tip = tips[i]
  spline$scaled.smooth[spline$scaled.smooth < 0] <- 0
  pdf(paste0("heatmap_", tip, "_viridis.pdf"), width=5, height=10)
  gplots::heatmap.2(x = as.matrix(spline$scaled.smooth[order, ]),
                    Rowv = F, Colv = F,
                    dendrogram = "none",
                    col = cols,
                    trace = "none",
                    density.info = "none", key = F,
                    cexCol = 0.8, cexRow = 0.15, margins = c(4, 4),
                    lwid = c(0.3, 4), lhei = c(0.3,4), labCol = NA)
  title(main = tip)
  dev.off()

  #tiff(paste0("Heatmaps_", tip, "-viridis.tiff"), width=3, height=10, units = 'in', res = 300)
  #gplots::heatmap.2(x = as.matrix(spline$scaled.smooth[order, ]),
  #                  Rowv = F, Colv = F,
  #                  dendrogram = "none",
  #                  col = cols,
  #                  trace = "none",
  #                  density.info = "none", key = F,
  #                  cexCol = 0.8, cexRow = 0.3,
  #                  lwid = c(0.3, 4), lhei = c(0.3,4), labCol = NA)
  #title(main = titles[i])
  #dev.off()
}



