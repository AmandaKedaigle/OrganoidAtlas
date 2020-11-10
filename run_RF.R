library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dplyr)
library(randomForest)
library(ROCR)
library(reshape2)

setwd("~/Documents/HumanData/Geschwind2019/RandomForest/")

###Utility functions
fitMethod <- function(dat)
{
  y=factor(dat[,"CellType"])
  print(levels(y))
  print(table(y))
  x <- dat %>% select(-CellType)
  ret=tuneRF(x,y,doBest=T)
  dev.off()
  return(ret)
}


valMethod <- function(model, dat, cols=NA, model_name="RF")
{
  pdf(paste0(model_name, "_model_ROC.pdf"))
  pred_probs <- predict(model, dat, type='prob')
  classes <- levels(dat$CellType)
  auc_vals <- c()
  for (i in 1:length(classes))
  {
    true_val <- ifelse(dat$CellType==classes[i], 1, 0)
    pred <- prediction(pred_probs[,i], true_val)
    perf <- performance(pred, "tpr", "fpr")
    if(i==1)
    {
      if (is.na(cols)) { plot(perf,main="ROC Curve") }
      else { plot(perf,main="ROC Curve", col=cols[classes[i]]) }
    }
    else
    {
      if (is.na(cols)) { plot(perf,main="ROC Curve",add=TRUE) }
      else { plot(perf,main="ROC Curve", col=cols[classes[i]],add=TRUE) }
    }
    auc.perf <- performance(pred, measure = "auc")
    print(auc.perf@y.values)
    auc_vals <- c(auc_vals, auc.perf@y.values[[1]])
  }
  if (!is.na(cols)) { legend("bottomright", legend=names(cols), col=cols, lty=1) }
  dev.off()
  
  names(auc_vals) <- classes
  write.table(auc_vals, paste0(model_name, "_model_AUC.tsv"))
}


testMethod <- function(model, dat)
{
  y=predict(model,dat)
  y_prob = predict(model, dat, type='prob')
  
  return(list(y, y_prob))
}


plotRes <- function(out, prob, prefix, cols=NA)
{
  ord <- order(out)
  out <- out[ord]
  prob <- prob[ord,levels(out)]
  anno = as.character(out)
  names(anno) <- names(out)
  anno <- as.data.frame(anno)
  colnames(anno) <- "Predicted Label"
  if (!is.na(cols)) { anno_cols <- list("Predicted Label"=cols) }
  
  pdf(paste0(prefix, "_prob_heatmap_rf.pdf"), height=8, width=12)
  if (is.na(cols)) { pheatmap(t(prob), legend=T, scale='column', show_colnames=F, main="Prediction Probabilities - column normalized", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno) }
  else { pheatmap(t(prob), legend=T, scale='column', show_colnames=F, main="Prediction Probabilities - column normalized", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols) }
  
  if (is.na(cols)) { pheatmap(t(prob), legend=T, scale='none', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno) }
  else { pheatmap(t(prob), legend=T, scale='none', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols) }
  dev.off()
  
}



# Import datasets
print("Loading data")

query <- readRDS("../../../wt_merge/FinalObjects/3mo_harmonizedObj_102820.rds")
prefix = "3m_harmonized"

system(paste0("mkdir -p ", prefix))
setwd(prefix)

ref <- readRDS("../../Geschwin_clusteredSeur.rds")
# read in a list of curated genes to use as the input features for RF
#var.genes <- readRDS("var_genes_3311_without_IEGs.rds")
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
val.dat <- dat1[-sampled.rows,]
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
model <- fitMethod(train.dat)
saveRDS(model, paste0(prefix, "_RF_model.rds"))

# Training null model
null_model <- fitMethod(null.train)
saveRDS(null_model, paste0(prefix, "_null_model.rds"))

# Validation
#colors <- readRDS("colors.rds")  # character array of colors, and the names must match the cell types
valMethod(model, val.dat, model_name="RF")
valMethod(null_model, val.dat, model_name="null")


# Test
res <- testMethod(model, test.dat)
out <- res[[1]]
prob <- res[[2]]
saveRDS(out, paste0(prefix, "_prediction_output.rds"))
saveRDS(prob, paste0(prefix, "_prediction_probs.rds"))

query$RF_prediction <- out
p1 = DimPlot(query, reduction='umap', pt.size=.1, group.by="RF_prediction", label=TRUE) + NoAxes() + ggtitle("Prediction")
p2 = DimPlot(query, reduction='umap', pt.size=.1, label=TRUE) + NoAxes() + ggtitle("Original")
png(paste0(prefix, "_prediction_result.png"), res=100, width = 1000, height=500)
p1 + p2
dev.off()

# Plot predicted labels vs probability scores
plotRes(out, prob, prefix)

#Plot agreement between our cell types and predicted labels
t = as.data.frame(table(query$CellType, query$RF_prediction))
t = melt(t)
png(paste0(prefix, "_compare_labels.png"), res=100, width=800, height=700)
ggplot(t, aes(Var1, Var2)) + 
  geom_point(aes(size = value, color = value)) + theme_bw() +
  xlab("Human Cell Types") + 
  ylab("Organoid Cell Types") + 
  geom_text(aes(label = value)) +
  scale_color_gradient(low="lightgray", high="violet", name="Number of Organoid Cells") +
  scale_size_continuous(name="Number of Organoid Cells") +
  guides(color=guide_legend(), size=guide_legend()) +
  theme(axis.text.x=element_text(angle=45, hjust=0.9))
dev.off()
