library(RCTD)
library(Seurat)
library(readxl)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggridges)
library(forcats)

setwd("~/Documents/SlideSeq/Exp1/1m_X8/sub1/")

load("../myRCTD_sub1.Robj")
myRCTD = myRCTD1
puck = myRCTD@spatialRNA

xmin = 1300 #1-2000 2-500  3-1300 4-1100 5-2000 6-1250 7-350  8-400  9-2200 10-1900 11-2600 14-1350 X8_1-1300 X8_2-3600 X8_3-2200 X3-1800 X4-2250
xmax = 2650 #1-3900 2-2600 3-3800 4-3600 5-4800 6-4000 7-3600 8-3800 9-3000 10-2750 11-4550 14-3800 X8_1-2650 X8_2-4650 X8_3-3500 X3-4200 X4-4200
ymin = 2900 #1-1500 2-2100 3-1300 4-1500 5-1300 6-1500 7-2400 8-1700 9-1430 10-1000 11-2300 14-700  X8_1-2900 X8_2-2100 X8-3-1900 X3-500  X4-2400
ymax = 4100 #1-3250 2-3900 3-3800 5-4000 5-3800 6-4800 7-5000 8-5000 9-2000 10-1900 11-4500 14-2850 X8_1-4100 X8_2-2900 X8-3-2900 X3-2800 X4-4200
plot_puck_continuous(puck, colnames(puck@counts), puck@nUMI, ylimit=c(0,max(puck@nUMI)/2), xlim = c(xmin, xmax), ylim=c(ymin,ymax))

pos = puck@coords
barcodes = rownames(pos)[(pos$x>xmin & pos$x<xmax) & (pos$y>ymin & pos$y<ymax)]
mean(puck@nUMI[names(puck@nUMI) %in% barcodes])


file.names = list.files(path = "~/Documents/AtlasAna/CompareSignatures/2mo_Signatures/", recursive = T, full.names = T, pattern = "*.xlsx")
names = sapply(strsplit(sapply(strsplit(file.names, "/"),"[",9), ".", fixed=T), "[", 1)

for (i in 1:length(file.names)) {
  markers = read_xlsx(file.names[[i]])
  markerGenes = markers[1:50, "gene"]$gene
  markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
  print(length(markerGenes))

  vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
  plot_puck_continuous(puck, colnames(puck@counts), vals, ylimit = c(0,quantile(vals, probs = 0.95)), xlim = c(xmin, xmax), ylim=c(ymin,ymax))+
      theme_void()+ scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=quantile(vals, probs = 0.3))
  ggsave(paste0("markers_",names[[i]],".pdf"), width=3.5, height=3)
}


areas = read_xlsx("../../Areas/Areas_Bhaduri_2021.xlsx")
names = unique(areas$area)
for (i in names) {
  markerGenes = areas[areas$area==i,"gene"]$gene
  markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
  print(length(markerGenes))
  
  vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
  plot_puck_continuous(puck, colnames(puck@counts), vals, ylimit = c(0,quantile(vals, probs = 0.95)), xlim = c(xmin, xmax), ylim=c(ymin,ymax))+
    theme_void()+ scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=quantile(vals, probs = 0.3))
  ggsave(paste0("markers_area_",i,".pdf"), width=3.5, height=3)
}

file.names = list.files(path = "~/Documents/HumanData/Geschwind2019/RandomForest/txtLists/", recursive = T, full.names = T, pattern = "*.txt")
names = sapply(strsplit(sapply(strsplit(file.names, "/"),"[",10), ".", fixed=T), "[", 1)

for (i in 1:length(file.names)) {
  markerGenes = read.table(file.names[[i]])$V1
  markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
  print(length(markerGenes))
  
  vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
  plot_puck_continuous(puck, colnames(puck@counts), vals, ylimit = c(0,quantile(vals, probs = 0.95)), xlim=c(xmin, xmax), ylim=c(ymin,ymax)) +
    theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkred", midpoint=quantile(vals, probs = 0.3))
  ggsave(paste0("markers_",names[[i]],".pdf"), width=3.5, height=3)
}


weights = myRCTD@results$weights
pdf("cellTypeWeights.pdf", width=3.5, height=3)
for (ct in colnames(weights)) {
   print(plot_puck_continuous(puck, colnames(puck@counts), weights[,ct], ylimit = c(0,quantile(weights[,ct], probs = 0.95)),xlim=c(xmin, xmax), ylim=c(ymin,ymax)) +
    theme_void() + scale_color_gradient2(low="white", mid="powderblue", high="darkblue", midpoint=quantile(weights[,ct], probs = 0.3)) + ggtitle(ct))
}
dev.off()

#Cell types colors & weights bar chart
results_df <- myRCTD@results$results_df
puck = myRCTD@spatialRNA
weights = myRCTD@results$weights
barcodes = rownames(results_df[results_df$spot_class != "reject" & puck@nUMI >= 100,])
pos = puck@coords[barcodes,]
barcodes = barcodes[(pos$x>xmin & pos$x<xmax) & (pos$y>ymin & pos$y<ymax)]
pos = pos[barcodes,]
weights = weights[barcodes,]
assigns = results_df[barcodes,]$first_type
pos$assign = assigns
celltypes = colnames(weights)
cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9', "#1D91C0", "#EC7014")
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown', "Newborn CFuPN", "Newborn CPN")
cols = cols[match(celltypes,levels)]
ggplot(pos, aes(x=x, y=y, color=assign)) +
  geom_point(size=1) + scale_color_manual(values=cols) + theme_void()
ggsave("AssignedCelltypes.png", height=2.5, width=4)
p = list()
for (ct in celltypes) {
  weightsub = as.matrix(weights[assigns==ct,])
  if (ncol(weightsub)>3) {
  summary = data.frame("CellType" = colnames(weightsub), "mean" = colMeans(weightsub), sd = apply(weightsub, 2, sd))
  p[[ct]] = ggplot(summary, aes(x=CellType, y = mean, fill=CellType)) +
    geom_bar(stat="identity", color="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
    scale_fill_manual(values=cols) +
    theme_classic() + NoLegend() + labs(x="", y="", title=paste0("Annotated ", ct))
  }
}
pdf("weightsDist.pdf", height=3*3, width=3*3.5)
wrap_plots(p)
dev.off()

#Metabolism -----
results_df <- myRCTD@results$results_df
weights = myRCTD@results$weights
barcodes = rownames(results_df[results_df$spot_class != "reject" & puck@nUMI >= 100,])
pos = puck@coords[barcodes,]
barcodes = barcodes[(pos$x>xmin & pos$x<xmax) & (pos$y>ymin & pos$y<ymax)]
pos = pos[barcodes,]
weights = weights[barcodes,]
assigns = results_df[barcodes,]$first_type
pos$assign = assigns
cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9', "#1D91C0", "#EC7014")
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown', "Newborn CFuPN", "Newborn CPN")
pos$cell = as.factor(rownames(pos))

#distance from edge defined as shortest of the distances to each line of convex hull
pos$distFromEdge = NA
hpts = chull(pos[,c('x','y')]) #gives the indices of points lying on the "convex hull" of the points
hpts = c(hpts, hpts[1])
plot(pos[,c('x','y')])
lines(pos[hpts,c('x','y')])
for (point in 1:nrow(pos)) {
  dist = Inf
  for (i in 1:(length(hpts)-1)) {
    pt1 = c(pos[hpts[i],'x'], pos[hpts[i],'y'])
    pt2 = c(pos[hpts[i+1],'x'], pos[hpts[i+1],'y'])
    poi = c(pos[point, 'x'], pos[point, 'y'])
    #distance between point poi and line defined by pt1 and pt2
    dist1 = abs(((pt2[1]-pt1[1])*(pt1[2]-poi[2])) - ((pt1[1]-poi[1])*(pt2[2]-pt1[2]))) / sqrt(((pt2[1]-pt1[1])^2) + ((pt2[2] - pt1[2])^2))
    if (dist1 < dist) { dist = dist1 }
  }
  pos[point, "distFromEdge"] = dist
}
ggplot(pos, aes(x=x, y=y, col=distFromEdge)) + geom_point()

#pos$x = NULL
#pos$y = NULL

file.names = list.files(path = "~/Documents/HumanData/Geschwind2019/RandomForest/txtLists/", recursive = T, full.names = T, pattern = "*.txt")
names = sapply(strsplit(sapply(strsplit(file.names, "/"),"[",10), ".", fixed=T), "[", 1)

for (i in 1:length(file.names)) {
  markerGenes = read.table(file.names[[i]])$V1
  markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
  print(length(markerGenes))
  
  vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
  pos[,names[i]] = scale(vals[rownames(pos)])
}

m = melt(pos[,3:ncol(pos)], id =c("cell","assign","distFromEdge"))
m$variable = factor(m$variable, levels = names[c(2:4,8:15,17,19:29,31,32,34:40,1,7,18,5,6,16,30,33)])
pathcols = c(rep("lightgray",32),'#A6CEE3','#33A02C','#B2DF8A','#1F78B4','#E78AC3','#E31A1C','#FDBF6F','#FF7F00')
pdf("Metabolism-allPathways-smooth.pdf",height=6, width=18)
ggplot(m, aes(x=distFromEdge, y=value, color=variable, fill=variable)) +
  geom_smooth(size=1, method = "loess") + theme_classic() + scale_color_manual(values=pathcols) + scale_fill_manual(values=pathcols) +
  labs(x = "Distance From Edge of Organoid", y="Scaled Normalized Gene Set Expression")
dev.off()
m = m[m$assign %in% names(table(m$assign)[table(m$assign)>400]),]
colsh = cols[match(levels(factor(m$assign)),levels)]
p1=ggplot(m[m$variable=="HALLMARK_HYPOXIA",], aes(x=distFromEdge, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance From Edge of Organoid", y="Scaled Normalized Hypoxia Expression")
p2=ggplot(m[m$variable=="HALLMARK_GLYCOLYSIS",], aes(x=distFromEdge, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance From Edge of Organoid", y="Scaled Normalized Glycolysis Expression")
p3=ggplot(m[m$variable=="HALLMARK_APOPTOSIS",], aes(x=distFromEdge, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance From Edge of Organoid", y="Scaled Normalized Apoptosis Expression")
pdf("Metabolism-3pathways-smooth.pdf", height=9, width=6.5)
p1/p2/p3
dev.off()

#distribution of cells
pos2 = pos[pos$assign %in% names(table(pos$assign)[table(pos$assign)>3]),]
pos2$assign = fct_reorder(pos2$assign, pos2$distFromEdge, median)
colsd = cols[match(levels(factor(pos2$assign)),levels)]
ggplot(pos2, aes(x=distFromEdge, y=assign, fill=assign)) +
  geom_density_ridges(size=0.3) + scale_fill_manual(values=colsd) + theme_classic() +
  labs(x = "Distance From Edge of Organoid", y="Density of Cells per Cell Type")
ggsave("DistributionOfCellTypes.pdf", height=5, width=6)

#Distribution from center of rosettes
#manually define centers based on IP and aRG signals

#1m X8_1:
markers = read_xlsx("~/Documents/AtlasAna/CompareSignatures/DESeq-organoids-1mo/IP.vs.Other.xlsx")
markerGenes = markers[1:50, "gene"]$gene
markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
names(vals) = colnames(puck@counts)
pos$IP_sig = vals[rownames(pos)]
pos$IP_sig[pos$IP_sig>1.7] = 1.7
markers = read_xlsx("~/Documents/AtlasAna/CompareSignatures/DESeq-organoids-1mo/aRG.vs.Other.xlsx")
markerGenes = markers[1:50, "gene"]$gene
markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
names(vals) = colnames(puck@counts)
pos$aRG_sig = vals[rownames(pos)]
pos$aRG_sig[pos$aRG_sig>1.7] = 1.7

centers = data.frame(x=c(1750, 2175, 2550, 2525, 2185, 1850, 2300, 2050), 
                     y=c(3400, 3275, 3175, 3500, 3600, 3700, 3850, 3950))

pos2 = pos[(pos$x>1600 & pos$x<1900 & pos$y>3200 & pos$y<3550) |
           (pos$x>1970 & pos$x<2350 & pos$y>3100 & pos$y<3450) |
           (pos$x>2400 & pos$y>3050 & pos$y<3300) |
           (pos$x>2350 & pos$y>3350 & pos$y<3670) |
           (pos$x>2030 & pos$x<2350 & pos$y>3300 & pos$y<3700) |
           (pos$x>1700 & pos$x<2000 & pos$y>3570 & pos$y<3850) |
           (pos$x>1900 & pos$x<2185 & pos$y>3750) |
           (pos$x>2185 & pos$x<2450 & pos$y>3750 & pos$y<4050),]
pos = pos2

ggplot(pos, aes(x=x, y=y, color=IP_sig)) +
  geom_point(size=0.5) + theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=0.5) +
  geom_point(data = centers, aes(x=x, y=y), color="black", size=2)
ggsave("RosetteCenters-IP-boxes.tiff")
ggplot(pos, aes(x=x, y=y, color=aRG_sig)) +
  geom_point(size=0.5) + theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=0.5) +
  geom_point(data = centers, aes(x=x, y=y), color="black", size=2)
ggsave("RosetteCenters-aRG-boxes.tiff")

#2m X3
markers = read_xlsx("~/Documents/AtlasAna/CompareSignatures/2mo_Signatures/IP.vs.Other.DEGs.xlsx")
markerGenes = markers[1:50, "gene"]$gene
markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
names(vals) = colnames(puck@counts)
pos$IP_sig = vals[rownames(pos)]
pos$IP_sig[pos$IP_sig>4] = 4
markers = read_xlsx("~/Documents/AtlasAna/CompareSignatures/2mo_Signatures/aRG.vs.Other.DEGs.xlsx")
markerGenes = markers[1:50, "gene"]$gene
markerGenes = markerGenes[markerGenes %in% rownames(puck@counts)]
vals = 100*colSums(as.matrix(puck@counts[markerGenes,]))/puck@nUMI
names(vals) = colnames(puck@counts)
pos$aRG_sig = vals[rownames(pos)]
pos$aRG_sig[pos$aRG_sig>4] = 4

centers = data.frame(x=c(3140, 2580, 2360, 2400, 2400, 2550, 2900, 3300, 3700, 4000, 4020, 4000, 3900, 3650), 
                     y=c(800, 900, 1050, 1400, 1750, 2100, 2400, 2450, 2350, 2100, 1750, 1350, 1000,850))

ggplot(pos, aes(x=x, y=y, color=aRG_sig)) +
  geom_point(size=0.5) + theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=1.75) +
  geom_point(data = centers, aes(x=x, y=y), color="black", size=2)
ggsave("RosetteCenters-aRG.tiff")
ggplot(pos, aes(x=x, y=y, color=IP_sig)) +
  geom_point(size=0.5) + theme_void() + scale_color_gradient2(low="white", mid="lightblue", high="darkmagenta", midpoint=1) +
  geom_point(data = centers, aes(x=x, y=y), color="black", size=2)
ggsave("RosetteCenters-IP.tiff")

#plots!
pos$distFromCenter = NA
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
for (point in 1:nrow(pos)) {
  dist = Inf
  for (i in 1:nrow(centers)) {
    dist1 <- euc.dist(pos[point,c('x','y')],centers[i,])
    if (dist1 < dist) { dist = dist1 }
  }
  pos[point, "distFromCenter"] = dist
}


m = melt(pos[,c(3:4,6:48)], id =c("cell","assign","distFromCenter"))
m = m[m$variable != "aRG_sig" & m$variable != "IP_sig",]
#m = m[m$assign!="Cortical hem" & m$assign != "Subcortical neurons" & m$assign != "Subcortical progenitors" & m$assign != "Cajal Retzius",]
m$variable = factor(m$variable, levels = names[c(2:4,8:15,17,19:29,31,32,34:40,1,7,18,5,6,16,30,33)])
pathcols = c(rep("lightgray",32),'#A6CEE3','#33A02C','#B2DF8A','#1F78B4','#E78AC3','#E31A1C','#FDBF6F','#FF7F00')
pdf("RosetteCenter-Metabolism-allPathways-smooth-boxes.pdf",height=6, width=18)
ggplot(m, aes(x=distFromCenter, y=value, color=variable, fill=variable)) +
  geom_smooth(size=1, method = "loess") + theme_classic() + scale_color_manual(values=pathcols) + scale_fill_manual(values=pathcols) +
  labs(x = "Distance To Closest Rosette Center", y="Scaled Normalized Gene Set Expression")
dev.off()

m = melt(pos[,c(3:4,6:45,48)], id =c("cell","assign","distFromCenter"))
m = m[m$assign %in% names(table(m$assign)[table(m$assign)>400]),]
cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9', "#1D91C0", "#EC7014")
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown', "Newborn CFuPN", "Newborn CPN")
colsh = cols[match(levels(factor(m$assign)),levels)]
p1=ggplot(m[m$variable=="HALLMARK_HYPOXIA",], aes(x=distFromCenter, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance To Closest Rosette Center", y="Scaled Normalized Hypoxia Expression")
p2=ggplot(m[m$variable=="HALLMARK_GLYCOLYSIS",], aes(x=distFromCenter, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance To Closest Rosette Center", y="Scaled Normalized Glycolysis Expression")
p3=ggplot(m[m$variable=="HALLMARK_APOPTOSIS",], aes(x=distFromCenter, y=value, color=assign, fill=assign)) +
  geom_smooth(size=1, method="loess") + scale_color_manual(values=colsh) + theme_classic() + scale_fill_manual(values=colsh) +
  labs(x = "Distance To Closest Rosette Center", y="Scaled Normalized Apoptosis Expression")
pdf("RosetteCenter-Metabolism-3pathways-smooth.pdf", height=9, width=6.5)
p1/p2/p3
dev.off()

pos2 = pos[pos$assign %in% names(table(pos$assign)[table(pos$assign)>3]),]
pos2$assign = fct_reorder(pos2$assign, pos2$distFromCenter, median)
colsd = cols[match(levels(factor(pos2$assign)),levels)]
ggplot(pos2, aes(x=distFromEdge, y=assign, fill=assign)) +
  geom_density_ridges(size=0.3) + scale_fill_manual(values=colsd) + theme_classic() +
  labs(x = "Distance From Closest Rosette Center", y="Density of Cells per Cell Type")
ggsave("RosetteCenter-DistributionOfCellTypes.pdf", height=5, width=6)
