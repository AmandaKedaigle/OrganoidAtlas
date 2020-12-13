##Module score

gene.set <- read.table("~/directory/genelist.txt")$V1

#alternatively, create a vector with your genes of interest 
gene.set = c("SOX2", "EMX1") 

#module score function
obj = AddModuleScore(obj, features=list(gene.set))

#plot the results!
FeaturePlot(obj, "Cluster1")
VlnPlot(obj, "Cluster1")