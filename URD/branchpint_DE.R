.libPaths(.libPaths()[2])

library(URD)
library(dplyr)


# import tree
tree = readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/res_URD/trees/dm-100-20/tm-40-80/dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05/obj_urd_tree_dmpreference_vt0.7_mv10_cppb80_bppw5_pt0.05.rds")

if (!dir.exists("branchpoint_DE")) {
dir.create("branchpoint_DE")
}
setwd("branchpoint_DE")

# function to select cells at each branching point
cellsNearBranchpoint = function(urd.obj, parent.seg, child.seg1, child.seg2, spread=0.04) {
p.cells = cellsInCluster(urd.obj, "segment", parent.seg)
c1.cells = cellsInCluster(urd.obj, "segment", child.seg1)
c2.cells = cellsInCluster(urd.obj, "segment", child.seg2)
p.obj = urdSubset(urd.obj, p.cells)
c1.obj = urdSubset(urd.obj, c1.cells)
c2.obj = urdSubset(urd.obj, c2.cells)

ptime = max(p.obj@pseudotime$pseudotime)
print("pseudotime of this branching point:")
print(ptime)

min.ptime = ptime - spread
max.ptime = ptime + spread

cells = c()
for (obj in c(p.obj, c1.obj, c2.obj)) {
cur.cells = rownames(obj@pseudotime)[which(obj@pseudotime$pseudotime > min.ptime & obj@pseudotime$pseudotime < max.ptime)]
cells = c(cells, cur.cells)
}

return(cells)
#return(urdSubset(tree, cells))
}

#bp1 = subsetURDBranchpoint(tree, 7, 5, 6)
#bp2 = subsetURDBranchpoint(tree, 6, 1, 4)
#bp3 = subsetURDBranchpoint(tree, 5, 2, 3)

library(Seurat)
# import Seurat
seur <- readRDS("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/data/seur_RNA_matching_barcodes.rds")

# my branching points
p.arr = c(3)
c1.arr = c(1)
c2.arr = c(2)

library(glmnet)
library(openxlsx)

# find DEGs for each segment connected to a branching point
for (i in 1:length(p.arr)) {
print(paste0("working on branching point ", i))

# select cells near bp
bp.cells = cellsNearBranchpoint(tree, p.arr[i], c1.arr[i], c2.arr[i])

# subset data
sseur <- subset(seur, cells=bp.cells)
bp.tree <- urdSubset(tree, bp.cells)

# add segment info onto Seurat metadata
meta.df = cbind(bp.tree@group.ids, bp.tree@pseudotime)
sseur = AddMetaData(sseur, meta.df)
Idents(sseur) <- "segment"

# find variable features
sseur <- FindVariableFeatures(sseur)

# lasso regression to select genes vary with pseudotime
print("Lasso regression on variable genes")
y = sseur$pseudotime
x = data.matrix(t(sseur@assays$RNA@data[VariableFeatures(sseur),]))
cv.model = cv.glmnet(x, y, alpha=1)
best.lambda <- cv.model$lambda.min  # lambda that minimizes MSE
print(best.lambda)
pdf(paste0("plot_MSE_by_lambda_value_bp", i, ".pdf"))
print(plot(cv.model))
dev.off()

best.model = glmnet(x, y, alpha=1, lambda=best.lambda)
saveRDS(best.model, paste0("model_lasso_bp", i, ".rds"))
print(dim(coef(best.model)))

var.feas = coef(best.model)
select.feas = setdiff(rownames(var.feas)[which(var.feas[,1] != 0)], "(Intercept)")

# compute DEGs for each segment related to bp
markers <- FindAllMarkers(sseur, test.use="wilcox", logfc.threshold=0.25, only.pos=T)
saveRDS(markers, paste0("markers_bp", i, "_wilcox.rds"))
marker.feas <- markers %>% subset(p_val_adj < 0.05) %>% pull(gene)

# selected variable features + DEGs will be used to fit gradient boosting classifier
select.feas = unique(c(select.feas, marker.feas))
saveRDS(select.feas, paste0("bp", i, "_genes.rds"))
dat <- as.data.frame(t(sseur@assays$RNA@data[select.feas,]))
dat[["class"]] <- sseur$segment
print("saving feature matrix!")
write.table(dat, paste0("data_varGenesLasso+DEG_bp", i, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)

# export to excel
outfname = paste0("markers_bp", i, "_wilcox_sig_only.xlsx")
wb = createWorkbook()
for (cls in levels(markers$cluster)) {
m <- markers %>% subset(cluster==cls) %>% subset(p_val_adj < 0.05)
if (nrow(m)>0) {
sheetName=paste0("segment",cls)
addWorksheet(wb, sheetName)
freezePane(wb, sheetName, firstRow=T)
writeDataTable(wb, sheetName, m[order(-m$avg_log2FC),])
saveWorkbook(wb, outfname, overwrite=T)
}
}
}

print("DONE!")
