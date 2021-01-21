library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(randomForest)
library(ROCR)
library(pheatmap)
library(patchwork)


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


valMethod <- function(model, dat, cols='', model_name="RF")
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
if (cols=='') { plot(perf,main="ROC Curve") }
else { plot(perf,main="ROC Curve", col=cols[classes[i]]) }
}
else
{
if (cols=='') { plot(perf,main="ROC Curve",add=TRUE) }
else { plot(perf,main="ROC Curve", col=cols[classes[i]],add=TRUE) }
}
auc.perf <- performance(pred, measure = "auc")
print(auc.perf@y.values)
auc_vals <- c(auc_vals, auc.perf@y.values[[1]])
}
if (cols!='') { legend("bottomright", legend=names(cols), col=cols, lty=1) }
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

#if (is.na(cols)) { p1 = pheatmap(t(prob), legend=T, scale='column', show_colnames=F, main="Prediction Probabilities - column normalized", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno) }
#else { p1 = pheatmap(t(prob), legend=T, scale='column', show_colnames=F, main="Prediction Probabilities - column normalized", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols) }

png(paste0(prefix, "_prob_heatmap_rf.png"), res=100, width=800, height=800)
if (is.na(cols)) { pheatmap(t(prob), legend=T, scale='none', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno) }
else { pheatmap(t(prob), legend=T, scale='none', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols) }
dev.off()

}
