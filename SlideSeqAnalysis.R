library(Seurat)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(gplots)
library(RCTD)
library(Matrix)
library(xlsx)

#by Xin Jin & from RTCD code

### RTCD need R-3.5

reference <- readRDS("SCRef.RDS")
puck <- read.SpatialRNA(paste0(filename, "/RTCD")) # read in the SpatialRNA object

# subset from a puck
#given a puck object, returns a puck with counts filtered based on UMI threshold and gene list
restrict_counts <- function(puck, gene_list, UMI_thresh = 1, UMI_max = 20000) {
  keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
  puck@coords = puck@coords[keep_loc,]
  puck@counts = puck@counts[gene_list,keep_loc]
  if(length(puck@cell_labels) > 0) #check cell_labels non null
    puck@cell_labels = puck@cell_labels[keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}
# crop 
restrict_y <- function(puck, gene_list, ymin = 1, ymax = 20000) {
  keep_loc = (puck@coords$y > ymin) & (puck@coords$y < ymax)
  puck@counts = puck@counts[gene_list,keep_loc]
  puck@coords = puck@coords[keep_loc,]
  if(length(puck@cell_labels) > 0) #check cell_labels non null
    puck@cell_labels = puck@cell_labels[keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}

# crop 
restrict_x <- function(puck, gene_list, xmin = 1, xmax = 20000) {
  keep_loc = (puck@coords$x > xmin) & (puck@coords$x < xmax)
  puck@counts = puck@counts[gene_list,keep_loc]
  puck@coords = puck@coords[keep_loc,]
  if(length(puck@cell_labels) > 0) #check cell_labels non null
    puck@cell_labels = puck@cell_labels[keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}

subpuck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 200) #21896 dots
subpuck <- restrict_y(subpuck, rownames(puck@counts), ymin =0, ymax=3000)  # y was 531-5697
subpuck <- restrict_x(subpuck, rownames(puck@counts), xmin =1500, xmax=4200)  # y was 531-5697

# barcodes <- colnames(puck@counts) #pixels to be used (a list of barcode names). 
plot_puck_continuous(puck=subpuck, barcode=colnames(subpuck@counts)[subpuck@nUMI >= 200], plot_val=puck@nUMI, ylimit = c(0,900), 
                     title ='subset plot of nUMI')
ggsave("201210_half_puck_above200.pdf")

plot_puck_continuous(puck=puck, barcode=colnames(puck@counts)[puck@nUMI >= 200], plot_val=puck@nUMI, ylimit = c(0,900), 
                     title ='subset plot of nUMI')
ggsave("201210_full_puck_above200.pdf")

# create RTCD
myRCTD <- create.RCTD(subpuck, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = TRUE)

save(myRCTD, file = "201210_myRCTD.Robj")

## plot
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- '201210_RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)

# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
ggsave("201210_Plot1.pdf")
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
ggsave("201210_Plot2.pdf")

# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
ggsave("201210_Plot3.pdf")

# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
ggsave("201210_Plot4.pdf")



# find marker genes
library(ggpubr)

load("201210_myRCTD.Robj")
puck <- myRCTD@spatialRNA

cur_cell_types <- c("CFuPN","Newborn CPN","Newborn CFuPN","IP",
                    "oRG","aRG","PN","Cortical hem")
cell_type_info_restr = myRCTD@cell_type_info$info
cell_type_info_restr[[1]] = cell_type_info_restr[[1]][,cur_cell_types]
cell_type_info_restr[[2]] = cur_cell_types; cell_type_info_restr[[3]] = length(cur_cell_types)

get_marker_data <- function(cell_type_names, cell_type_means, gene_list) {
  marker_means = cell_type_means[gene_list,]
  marker_norm = marker_means / rowSums(marker_means)
  marker_data = as.data.frame(cell_type_names[max.col(marker_means)])
  marker_data$max_epr <- apply(cell_type_means[gene_list,],1,max)
  colnames(marker_data) = c("cell_type",'max_epr')
  rownames(marker_data) = gene_list
  marker_data$log_fc <- 0
  epsilon <- 1e-9
  for(cell_type in unique(marker_data$cell_type)) {
    cur_genes <- gene_list[marker_data$cell_type == cell_type]
    other_mean = rowMeans(cell_type_means[cur_genes,cell_type_names != cell_type])
    marker_data$log_fc[marker_data$cell_type == cell_type] <- log(epsilon + cell_type_means[cur_genes,cell_type]) - log(epsilon + other_mean)
  }
  return(marker_data)
}

de_genes <- get_de_genes(cell_type_info_restr, puck, fc_thresh = 3, expr_thresh = .0001, MIN_OBS = 3)
marker_data_de = get_marker_data(cell_type_info_restr[[2]], cell_type_info_restr[[1]], de_genes)
saveRDS(marker_data_de, '201210_marker_data_de_standard.RDS')
write.csv(marker_data_de, "201210_markergene.csv")

results_df <- results$results_df
weights_doublet <- results$weights_doublet
marker_data_de <- readRDS('201210_marker_data_de_standard.RDS')
barcodes = rownames(results_df[results_df$spot_class != "reject" & puck@nUMI >= 100,])
my_table = puck@coords[barcodes,]
my_table$class = results_df[barcodes,]$first_type
n_levels = cell_type_info_restr[[3]]
my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
names(my_pal) = cell_type_info_restr[[2]]
my_pal_curr <- my_pal
pres = unique(as.integer(my_table$class))
pres = pres[order(pres)]
p1 <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + 
  ggplot2::geom_point(ggplot2::aes(size = .5, shape=19,color=class)) + 
  ggplot2::scale_color_manual("",values = my_pal_curr[pres], breaks = cur_cell_types, labels = cur_cell_types)+ 
  ggplot2::scale_shape_identity() + 
  ggplot2::theme_classic() + ggplot2::scale_size_identity() + 
  coord_fixed() + theme(legend.position="top")+ 
  guides(colour = guide_legend(override.aes = list(size=2)))+ 
  scale_x_continuous(breaks = c(1000,3000,5000), limits = c(900,5600)) + 
  scale_y_continuous(breaks = c(500,2000,3500), limits = c(500,3400))+ 
  geom_segment(aes(x = 1300, y = 700, xend = 1684.6, yend = 700), color = "black")+  
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggarrange(p1)
ggsave("201210_plot.pdf")

# plot doublet efficiency
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",]
doub_occur <- table(doublets$second_type, doublets$first_type)
rest_cell <- cur_cell_types
doub_chart <- doub_occur[rest_cell, rest_cell]
diag(doub_chart) <- 0.5 * table(results_df[results_df$spot_class == "singlet","first_type"])[rest_cell]
doub_chart[doub_chart > log(100,2)] <- log(100,2)
library(reshape2)
data <- melt(as.matrix(doub_chart))
colnames(data) = c('Prediction','Reference','value')
p4 <- ggplot(data, aes(Reference, Prediction, fill= value)) +  geom_tile() +theme_classic() +scale_fill_gradientn(colors = pals::brewer.blues(20)[2:20],name = "Log Doublet Count", labels = c(1,6),breaks = c(1,6))+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Cell Type 1')+ ylab('Cell Type 2')+ theme(legend.position="top")
ggarrange(p4)
ggsave("201210_double.pdf")

## cluster and re do DE analysis
load("201210_myRCTD.Robj")
puck <- myRCTD@spatialRNA

target <- CreateSeuratObject(puck@counts, project = "Slideseq",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
target[['image']] <- new(Class = 'SlideSeq', assay = 'Spatial', coordinates = puck@coords)

target <- NormalizeData(target, verbose = FALSE)
target <- FindVariableFeatures(target, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
target <- ScaleData(target)
target <- RunPCA(target, assay = "RNA", verbose = FALSE)
target <- RunUMAP(target, dims = 1:30)
target <- FindNeighbors(target, dims = 1:30)
target <- FindClusters(target, resolution = 1.2, verbose = FALSE)
puck@cell_labels <- target@active.ident
n_clusters = length(levels(puck@cell_labels))

p1 <- DimPlot(target,order = TRUE, label=T,pt.size = 0.05)
p2 <- SpatialDimPlot(target, stroke = NA, alpha = 1, pt.size = 1.5)
p <- grid.arrange(p1, p2, nrow = 1)
ggsave(filename = "201210_tsne_and_spatial.pdf", plot = p, width = 15)

save(target, file = "201210_target.Robj")

FeaturePlot(target, features = c("SOX2","MAP2"))
ggsave(filename = "201210_seurat_marker.pdf", width = 15)

m = FindAllMarkers(object = target, only.pos = TRUE)
write.csv(x = m, file = "201210_seuratmarkers.csv")

#Plot scRNA signature genes
puck = myRCTD@spatialRNA

file.names = list.files(path = "3mscRNAsignatures", recursive = T, full.names = T, pattern = "*.xlsx")
names = sapply(strsplit(sapply(strsplit(file.names, "/"),"[",2), ".", fixed=T), "[", 1)

for (i in 1:length(file.names)) {
  markers = read_xlsx(file.names[[i]])
  markerGenes = markers[1:50, "gene"]$gene
  markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
  print(length(markerGenes))
  
  vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
  plot_puck_continuous(puck, colnames(puck@counts), vals, ylimit = c(0,quantile(vals, probs = 0.95))) + #, xlim = c(2100, 4200)) +
    theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=quantile(vals, probs = 0.3))
  ggsave(paste0("markers_",names[[i]],".pdf"), width=3.5, height=3)
}

#plot RTCD weights
weights = myRCTD@results$weights
pdf("cellTypeWeights.pdf", width=3.5, height=3)
for (ct in colnames(weights)) {
  print(plot_puck_continuous(puck, colnames(puck@counts), weights[,ct], ylimit = c(0,quantile(weights[,ct], probs = 0.95)), xlim = c(2100, 4200)) +
          theme_void() + scale_color_gradient2(low="white", mid="powderblue", high="darkblue", midpoint=quantile(weights[,ct], probs = 0.3)) + ggtitle(ct))
}
dev.off()

