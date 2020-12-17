library(Seurat)
library(WGCNA)

choosePower <- function(seur, powers=NULL, genesToUse=NULL, save=T, plot=T, lowCutoff=0.15, highCutoff=9)
{

print("gene expression cutoffs:")
print(paste0("lowCutoff=", lowCutoff))
print(paste0("highCutoff=", highCutoff))

dat=seur@assays$RNA@data

if (is.null(genesToUse)) {
dat<-dat[!(Matrix::rowMeans(dat, na.rm = TRUE) < lowCutoff),] #get rid of lowly expressed genes
dat<-dat[!(Matrix::rowMeans(dat, na.rm = TRUE) > highCutoff),] # get rid of highly expressed genes
dat <- dat[!(rownames(dat) %in% grep(pattern = "^RP", x = rownames(dat), value = TRUE)),] ##Remove ribosomal genes
dat <- dat[!(rownames(dat) %in% grep(pattern = "^MT-", x = rownames(dat), value = TRUE)),] ##Remove mito genes
} else {
dat <- dat[genesToUse, ]
}

dat = as.matrix(dat)

dat=t(dat)

if (is.null(powers)) {
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
}

# Call the network topology analysis function
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5)

if (save) {
saveRDS(sft, "pickSoftThreshold_res.rds")
}

if (plot) {
# Plot the results:
pdf("power_selection.pdf", width=9, height=5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
}

return(sft)
}
