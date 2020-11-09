#Harmony 

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

obj1 = readRDS("obj1.rds")
obj1$dataset = 1 #Add a metadata column to each object so that at the end you can keep track of which dataset each cell is from
obj2 = readRDS("obj2.rds")
obj2$dataset = 2
obj3 = readRDS("obj3.rds")
obj3$dataset = 3
obj4 = readRDS("obj4.rds")
obj4$dataset = 4
obj5 = readRDS("obj5.rds")
obj5$dataset = 5
obj6 = readRDS("obj6.rds")
obj6$dataset = 6
obj7 = readRDS("obj7.rds")
obj7$dataset = 7

combined = merge(obj1, list(obj2, obj3, obj4, obj5, obj6, obj7))

#Begin the analysis of "combined" object

combined = NormalizeData(combined, normalization.method = "LogNormalize", scale.factor=1000000)
combined = FindVariableFeatures(combined, selection.method='mean.var.plot')
combined <- ScaleData(combined, features=VariableFeatures(combined),vars.to.regress=c("nCount_RNA", "CC.Difference"))
combined <- RunPCA(combined)

#Use the Seurat wrapper function "RunHarmony" to do batch correction!
#This will take a few minutes and output several progress bars (usually around 3-10 steps are needed for convergence)
#"theta" is the diversity penalty parameter. Higher theta is more aggresive batch correction, theta=0 is basically no batch correction. Default is 2.
combined = RunHarmony(combined, group.by.vars = "dataset", theta=2)

combined = RunUMAP(combined, reduction="harmony", dims=1:30) #Just RunUMAP or RunTSNE using reduction="harmony" to see the overlapped!

#Save your harmonized object :)
saveRDS(combined, "harmonizedObj.rds")

#You can now make other umaps/tsnes and feature plots just as you would normally!
DimPlot(combined, group.by = 'dataset') #i.e., color cells by which dataset they came from

