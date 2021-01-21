library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dplyr)
library(randomForest)
library(ROCR)
library(reshape2)

args = commandArgs(trailingOnly = T) #prefix, reference, gene list

source("../RF_utils.R")

#Reference datasets have been built like this: ----
#ref <- readRDS("../../../../wt_merge/FinalObjects/6mo_harmonizedObj_103020.rds")
#ref$FinalName = as.character(ref$FinalName)
#ref = subset(ref, subset=FinalName !="Unknown")
#ref = subset(ref, subset=FinalName !="Cajal Retzius")
#ref = subset(ref, subset=FinalName !="Cortical hem")
#ref = subset(ref, subset=FinalName !="Subcortical neurons")
#ref = subset(ref, subset=FinalName !="Subcortical progenitors")
#Idents(ref) = "FinalName"
# downsample the reference dataset to have the equal number of cells in each cell type
#num <- min(table(Idents(ref)))
#cat(paste0("Number of cells to be sampled from each cluster: ", num, "\n"))
#ref<-subset(ref,downsample=num)
#saveRDS(ref, "1mo_downsampled.rds")

# Import datasets ----
#prefix = "1m_minusGlyc"
prefix = args[1]
ref = readRDS(args[2])
#genefile = "../../txtLists/HALLMARK_GLYCOLYSIS.txt"
genefile = args[3]

system(paste0("mkdir -p ", prefix))
setwd(prefix)

query <- readRDS("../../../Geschwin_clusteredSeur.rds")
var.genes = VariableFeatures(ref)
var.genes <- intersect(var.genes, rownames(query@assays$RNA@counts))

query$Cluster = as.character(query$Cluster)
'%ni%' = Negate('%in%')
query = subset(query, subset=Cluster %ni% c("End","Mic","Per"))
query$Cluster[query$Cluster=="ExDp1"] = "ExDp"
query$Cluster[query$Cluster=="ExDp2"] = "ExDp"

# Set the cell types in the active.idents slot of the reference seurat object
Idents(ref) = "FinalName"
Idents(query) = "Cluster"

f = file(paste0(prefix, "_genesUsed.txt"), "w")
writeLines(c("Variable Genes in Reference Object", length(var.genes)), f)
if (!genefile=="None") {
  add_genes = as.character(read.table(genefile, skip=2)$V1)
  writeLines(c("Number of Genes in Additional Gene List", length(add_genes)), f)
  ol = add_genes[add_genes %in% var.genes]
  writeLines(c("Overlapping Genes", ol), f)
  var.genes = var.genes[!var.genes %in% add_genes]
  writeLines(c("Number of Final Genes Used", length(var.genes)), f)
}
close(f)

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
saveRDS(null_model, paste0(prefix, "_rm_null_model.rds"))

# Validation
cols = c('#8dd3c7','#bebada', '#fb9a99', '#80b1d3', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#ffffb3','#fee391','#a1d99b', '#67a9cf','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
names = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
          'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical', #neurons
          'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
          'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
          'Glial precursors', 'IN progenitors','Immature IN','Unknown')
names(cols) = names
cols = cols[names(cols) %in% levels(val.dat$CellType)]
valMethod(model, val.dat, model_name="RF", cols = cols)

# Test
print("Test")
resRM <- testMethod(model, test.dat)
outRM <- resRM[[1]]
probRM <- resRM[[2]]
saveRDS(outRM, paste0(prefix, "_prediction_output.rds"))
saveRDS(probRM, paste0(prefix, "_prediction_probs.rds"))

query$RF_prediction <- outRM

# Plot predicted labels vs probability scores
plotRes(outRM, probRM, prefix)

#Plot agreement between our cell types and predicted labels
tRM = as.data.frame(table(query$Cluster,query$RF_prediction)/colSums(table(query$RF_prediction, query$Cluster)))
tRM = melt(tRM)
saveRDS(tRM, "tRM.rds")
pdf(paste0(prefix, "_compare_labels_pretty.pdf"), width=6.5, height=3.5, useDingbats=F)
ggplot(tRM, aes(Var1, Var2)) + 
  geom_point(aes(size = value, color = value)) + theme_bw() +
  xlab("Organoid Cell Types") + 
  ylab("Human Cell Types") + 
  #geom_text(aes(label = ifelse(value>0.005, round(value,2), ''))) +
  scale_color_gradient(low="lightgray", high="darkmagenta", name="Fraction of Human Cells") +
  scale_size_continuous(name="Fraction of Human Cells") +
  guides(color=guide_legend(), size=guide_legend()) +
  theme(axis.text.x=element_text(angle=45, hjust=0.9))
dev.off()
