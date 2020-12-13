##Plots

##DimPlot individual organoids

Idents(obj) <- "org"
numOrgs = 21
l = levels(factor(obj$org))
cols= c('colors', 'of', 'interest', 'up', 'to', 'numOrgs')
for (i in 1:numOrgs) {
  tiff(paste0('file',l[i], '.tiff'), res = 300,  width = 1600, height=1600)
  print(DimPlot(obj, cells = WhichCells(obj, idents=l[i]), pt.size = 0.05, cols = cols[i]) + NoAxes() + NoLegend())
  dev.off()
}

##FeaturePlot

genes <- c("genes", "of", "interest")

fp <- FeaturePlot(obj, features = genes, combine = F, pt.size = 0.05, cols = c("color1", "color2")) 

for (i in 1:length(fp)){
  fp[[i]] <- fp[[i]] + NoLegend() + NoAxes() + ggtitle(element_blank())
  tiff(paste0('GENES', i, ".tiff"), res = 300, width = 1600, height=1600)
  print(fp[[i]])
  dev.off()
}

##alternatively 

fp <- FeaturePlot(obj, features = genes, combine = F, pt.size = 0.05) 

for (i in 1:length(fp)){
  fp[[i]] <- fp[[i]] + NoAxes() + ggtitle(element_blank()) + scale_color_gradient2(low="color1",mid="color2",high="color3", midpoint = 0)
  
  tiff(paste0('GENES', i, ".tiff"), res = 300, width = 1600, height=1600)
  print(fp[[i]])
  dev.off()
}

##VlnPlot

p0 <- VlnPlot(obj, 
              
              features = c("genes", "of", "interest"), 
              idents = c("cells", "of", "interest"),
              ncol = 3, pt.size = 0, combine = F, cols = c("colors", "of", "interest"))

for(j in 1:length(p0)) { 
  p0[[j]] <- p0[[j]] + 
    NoLegend() + 
    coord_flip() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          plot.title = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          axis.line.y = element_blank())
  
}

tiff("filename.tiff", res = 300, width = 2200, height = 1500)
print(finalplot)
dev.off()

##alternatively
genes <- c("genes", "of", "interest",)
p0 <- VlnPlot(obj,  features = genes, pt.size = 0, cols = c('color1','color2','color3'))

for (i in 1:length(p0)){
  p0[[i]]$layers[[1]]$aes_params$size <- 0.2
  p0[[i]] <- p0[[i]]+ ggtitle(element_blank()) +
    NoLegend() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=5, angle = 45), 
          axis.ticks.y = element_line(size = 0.2), 
          axis.ticks.x = element_line(size = 0.2), 
          axis.line.x = element_line(size = 0.2),
          axis.line.y = element_line(size = 0.2))
  
  pdf(paste0('filename', i, ".pdf"), width = 1.5, height = 1.5)
  print(p0[[i]])
  dev.off()   
}


