##Barchart
library(Seurat)
library(ggplot2)
library(reshape2)

setwd("~/directory/")

#Sets working directory. Change to a folder on your computer
obj= readRDS("obj.rds") #Reads in a rds (R dataset) object. Change filename to the object you want

counts = as.matrix(table(obj$CellTypes, obj$orig.ident)) #set so that first column is clusters/celltypes, second is organoids, or however you want to split up the x-axis
counts = t(t(counts)/colSums(counts)) #transform from raw counts to percentages, to normalize for library size
counts = melt(counts,varnames = c('CellType','Org')) #you need the reshape2 library(from above) for this - formats for ggplot
counts$Cluster = factor(as.character(counts$Cluster), levels=as.character(0:100)) #If your clusters are numbers (as opposed to words like Cell Types) run this line

#plot!
p1 <- ggplot(counts, aes(fill= CellType , y=value, x=Org)) + 
  geom_bar(stat="identity", position="fill", show.legend=T)+
  scale_fill_manual(values= c('chose_your_colors')) 

p1 <- p1 + theme_minimal(base_line_size = 0) + NoLegend() + NoAxes()

tiff('filename.tiff', res = 300, width = 1100, height=1000)
print(p1) 
dev.off()
