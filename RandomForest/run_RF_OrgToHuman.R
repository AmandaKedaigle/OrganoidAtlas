library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dplyr)
library(randomForest)
library(ROCR)
library(reshape2)

args = commandArgs(trailingOnly = T) #prefix, query, gene list

source("../RF_utils.R")

# Import datasets
prefix = args[1] #i.e. "1m"
query = readRDS(args[2]) #i.e. "1mo_111520_final.rds"

system(paste0("mkdir -p ", prefix))
setwd(prefix)

ref <- readRDS("../../../Geschwin_clusteredSeur.rds")
var.genes = VariableFeatures(ref)
var.genes <- intersect(var.genes, rownames(query@assays$RNA@counts))
query@assays$RNA@var.features <- var.genes
ref@assays$RNA@var.features <- var.genes

# Set the cell types in the active.idents slot of the reference seurat object
ref = subset(ref, cells=rownames(ref@meta.data)[!is.na(ref$Cluster)]) #removing 10 cells
ref = subset(ref, cells=rownames(ref@meta.data)[!ref$Cluster %in% c("Mic","Per","End")])
ref$Cluster = as.character(ref$Cluster)
ref$Cluster[ref$Cluster=="ExDp1"] = "ExDp"
ref$Cluster[ref$Cluster=="ExDp2"] = "ExDp"
Idents(ref) <- ref$Cluster

cat("Reference data cell types:\n")
print(table(Idents(ref)))

# downsample the reference dataset to have the equal number of cells in each cell type
num <- min(table(Idents(ref)))
cat(paste0("Number of cells to be sampled from each cluster: ", num, "\n"))
ref<-subset(ref,downsample=num)

# save randomly sampled cells (just in case, since I don't set random seed here)
write.table(colnames(ref), file="sampled_ref_cells.tsv", sep='\t', quote=F, row.names=F, col.names=F)
print(table(Idents(ref)))


# Prepare datasets to be used in training/testing
dat1 <- data.frame(t(as.matrix(ref@assays$RNA@data[var.genes,])))
dat1["CellType"]=ref@active.ident  # add label
test.dat <- data.frame(t(as.matrix(query@assays$RNA@data[var.genes,])))

# Split training data for validation
train.size <- floor(0.8*nrow(dat1))
sampled.rows <- sample(seq_len(nrow(dat1)), size=train.size)
train.dat <- dat1[sampled.rows,]
saveRDS(train.dat, "trainingData.rds")
val.dat <- dat1[-sampled.rows,]
saveRDS(val.dat, "valData.rds")
print("Training dataset size")
print(dim(train.dat))
print("Validation dataset size")
print(dim(val.dat))
print("Test dataset size")
print(dim(test.dat))

# Also prepare a null training dataset by shuffling labels
null.train <- train.dat
labels <- levels(train.dat$CellType)
nLabel <- nlevels(train.dat$CellType)
sampled.labels <- sample(1:nLabel, nrow(train.dat), replace=T)
null.train$CellType <- as.factor(labels[sampled.labels])

# Training
print("Train")
model <- fitMethod(train.dat)
saveRDS(model, paste0(prefix, "_RF_model.rds"))

# Training null model
print("Train Null")
null_model <- fitMethod(null.train)
saveRDS(null_model, paste0(prefix, "_null_model.rds"))

# Validation
#human colors
names = c('vRG', 'IP', 'oRG', 'ExN', 'ExM', 'ExM-U', 'ExDp', 'OPC', 'InMGE','InCGE')
cols = c('#66C2A5', '#8DA0CB', '#FFD92F', '#FEE0D2','#FC9272','#DE2D26', '#1F78B4', '#016C59', '#FA9FB5','#F768A1') #'#980043','#969696', '#737373' Mic Per End
names(cols) = names
#colors <- readRDS("colors.rds")  # character array of colors, and the names must match the cell types
valMethod(model, val.dat, model_name="RF", cols = cols)
valMethod(null_model, val.dat, model_name="null", cols = cols)


# Test
print("Test")
res <- testMethod(model, test.dat)
out <- res[[1]]
prob <- res[[2]]
saveRDS(out, paste0(prefix, "_prediction_output.rds"))
saveRDS(prob, paste0(prefix, "_prediction_probs.rds"))

query$RF_prediction <- out
p1 = DimPlot(query, reduction='umap', pt.size=.1, group.by="RF_prediction", label=TRUE) + NoAxes() + ggtitle("Prediction")
p2 = DimPlot(query, reduction='umap', pt.size=.1, group.by="FinalName", label=TRUE) + NoAxes() + ggtitle("Original")
png(paste0(prefix, "_prediction_result.png"), res=100, width = 1000, height=500)
p1 + p2
dev.off()

# Plot predicted labels vs probability scores
plotRes(out, prob, prefix)

#Plot agreement between our cell types and predicted labels
orgnames = query$FinalName
t = as.data.frame(table(orgnames,out)/colSums(table(out, orgnames)))
t = t[t$orgnames != "Unknown",]
t$orgnames = factor(t$orgnames, levels = c("aRG","oRG", "oRG II", "oRG/Astroglia", "Cortical hem", "Subcortical neurons",
                                           "Subcortical progenitors","IP", "Newborn PN", "Newborn DL PN",
                                           "PN", "Cajal Retzius", "CPN","CFuPN", "IN progenitors", "Immature IN", 
                                           "Glial precursors","Astroglia"))
t$out = factor(t$out, levels = c("PgG2M", "PgS","vRG", "oRG", "IP", "ExN", "ExM", "ExM-U","ExDp", "InCGE", "InMGE", "OPC"))

pdf(paste0(prefix, "_compare_labels_pretty.pdf"), width=6.5, height=4.5, useDingbats = F)
ggplot(t, aes(orgnames, out)) + 
  geom_point(aes(size =  Freq, color = Freq)) + theme_bw() +
  ylab("Human Cell Types") + 
  xlab("Organoid Cell Types") + 
  scale_color_gradient(low="lightgray", high="darkmagenta", name="Fraction of Organoid Cells") +
  scale_size_continuous(name="Fraction of Organoid Cells") +
  guides(color=guide_legend(), size=guide_legend()) +
  theme(axis.text.x=element_text(angle=45, hjust=0.9))
dev.off()
