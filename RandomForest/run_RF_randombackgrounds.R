library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dplyr)
library(randomForest)
library(ROCR)
library(reshape2)

args = commandArgs(trailingOnly = T) #prefix, reference, job number

source("../RF_utils.R")

# Import datasets
#prefix = "3m_bg_HALLMARK_HYPOXIA"
prefix = args[1]
ref = readRDS(args[2])
#ref = readRDS("3mo_downsampled.rds")
metlist = as.character(read.table(args[3])$V1)
#metlist = as.character(read.table("../txtLists/HALLMARK_HYPOXIA.txt")$V1)
num = args[4] #go up to 5-10,000 background runs per list

dir.create(prefix)
setwd(prefix)

query <- readRDS("../../../Geschwin_clusteredSeur.rds")
var.genes = VariableFeatures(ref)
var.genes <- intersect(var.genes, rownames(query@assays$RNA@counts)) #792 genes

query$Cluster = as.character(query$Cluster)
'%ni%' = Negate('%in%')
query = subset(query, subset=Cluster %ni% c("End","Mic","Per"))
query$Cluster[query$Cluster=="ExDp1"] = "ExDp"
query$Cluster[query$Cluster=="ExDp2"] = "ExDp"

# Set the cell types in the active.idents slot of the reference seurat object
Idents(ref) = "FinalName"
Idents(query) = "Cluster"

#select genes to remove
#Expression matched to each list, method as in ModuleScores
data.avg <- Matrix::rowMeans(x = GetAssayData(ref))
data.avg <- data.avg[order(data.avg)]
data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = 24, labels = FALSE, right = FALSE) #24 default nbin from AddModuleScore
names(x = data.cut) <- names(x = data.avg)
features.use <- metlist
features.use <- features.use[features.use %in% names(data.cut)]
ctrl.use = character()
for (j in 1:length(features.use)) {
  newgene = names(sample(data.cut[which(x = data.cut == data.cut[features.use[j]])],size = 1))
  k = 1
  while (k<5) {
    if (newgene %in% ctrl.use) { newgene = names(sample(data.cut[which(x = data.cut == data.cut[features.use[j]])],size = 1)) }
    k=k+1
  }
  ctrl.use[[j]] =  newgene
}

var.genes = var.genes[!var.genes %in% ctrl.use]
saveRDS(ctrl.use, paste0(num, "_removed.genes.rds"))

# Prepare datasets to be used in training/testing
dat1 <- data.frame(t(as.matrix(ref@assays$RNA@data[var.genes,])))
dat1["CellType"]=ref@active.ident  # add label
test.dat <- data.frame(t(as.matrix(query@assays$RNA@data[var.genes,])))

# Split training data for validation
train.size <- floor(0.8*nrow(dat1))
sampled.rows <- sample(seq_len(nrow(dat1)), size=train.size)
train.dat <- dat1[sampled.rows,]
val.dat <- dat1[-sampled.rows,]
print("Training dataset size")
print(dim(train.dat))
# Training
print("Train")
model <- fitMethod(train.dat)
saveRDS(model, paste0(num, "_RF_model.rds"))

# Test
print("Test")
resRM <- testMethod(model, test.dat)
outRM <- resRM[[1]]
probRM <- resRM[[2]]

query$RF_prediction <- outRM

t = as.data.frame(table(query$Cluster,query$RF_prediction)/colSums(table(query$RF_prediction, query$Cluster)))
t = melt(t)
saveRDS(t, paste0(num, "_tBG.rds"))

