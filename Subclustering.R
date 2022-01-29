## Subset and scale new object
Idents(obj) <- "Ident_needed"
DimPlot(obj, label = T)
obj1 <- subset(obj, idents = c("cells.needed")) 
DimPlot(obj1, label = T)

## Normalize and scale
obj1 = NormalizeData(obj1, normalization.method = "LogNormalize", scale.factor=1000000) ## normalize the data
obj1 = FindVariableFeatures(obj1, selection.method='mean.var.plot', x.low.cutoff=1)
obj1 <- ScaleData(obj1, vars.to.regress = c("nCount_RNA","CC.Difference"), features = VariableFeatures(obj1)) ## scale the data

print(ElbowPlot(obj1))

obj1 = RunPCA(obj1,features = VariableFeatures(obj1))
print(ElbowPlot(obj1))
obj1 = FindNeighbors(obj1, dims=c(1:20)) 
obj1 = FindClusters(obj1, resolution=c(1,0.5, 0.2, 0.1))
obj1 <- RunTSNE(object = obj1, dims = c(1:20)) 
obj1 <- RunUMAP(obj1, dims = 1:20)

## Make some plots to decide resolution
DimPlot(obj1, group.by = c("RNA_snn_res.1", "RNA_snn_res.0.5", "RNA_snn_res.0.2", "RNA_snn_res.0.1"), label = T)
Idents(obj1) = 'RNA_snn_res.0.2' ## define resolution
DimPlot(obj1)

## FindMarkers
markers_obj1 <- FindAllMarkers(obj1, min.pct = 0.25, only.pos = T, rownames = F) ## only positive 
write.csv(markers_obj1, file = "findmarkersobj1_filename.csv")

##Give names
Idents(obj1) <- "RNA_snn_res.0.2" ## desired resolution needs to be active identity
new.cluster.ids <- c("cell1", "cell2") ## create a vector with new cell types in the same order i.e. cluster 1, position 1 
names(new.cluster.ids) <- levels(obj1) ## levels of obj1: active identities ## with this function you assign the names to the clusters using the vector generated before 
obj1 <- RenameIdents(obj1, new.cluster.ids) ## for each cell, changing the name of the level assign to the cells to the new names
obj1$new_name_column <- Idents(obj1) ## in metadata

saveRDS(obj1, file = "obj1_name") ## save obj1

## bringing identities from subclustering (obj1) to the original obj

Idents(obj) <- "Ident_needed"

## Generate a new column in the metadata
obj$NewColumn <- as.character(Idents(obj))

## Change the information of cells containing sub-cluster information in the new metadata column
obj$NewColumn[Cells(obj1)] <- paste(Idents(obj1))
DimPlot(obj, group.by = "NewColumn") #plot to verify the identities have been changed

## Subclustering information Kedaigle, Uzquiano et al.

#For the Mito 210 c1 b1 23 days data set, clusters 16, 6, and 4 were further subclustered (resolution=0.5), resulting in Preplate/Subplate cells and FOXG1- EMX1- neurons. 
#For the GM08330 b2 1 month data set, cluster 10 was further subclustered (resolution=0.5), resulting in Subcortical neurons and Cajal Retzius cells. 
#For the Mito 210 c1 b3 1 month data set, no subclustering. 
#For the Mito 210 c1 b4 1 month data set, cluster 10 was further subclustered (resolution=0.1), resulting in Newborn PN and Subcortical neurons. 
#For the Mito 210 c1 b5 1 month data set, cluster 6 was further subclustered (resolution=0.2), resulting in Subcortical neurons and Cajal Retzius cells. 
#For the Mito 210 c1 b1 2 months data set, cluster 1 was further subclustered (resolution=0.5), resulting in aRG and oRG. Cluster 5 was further subclustered (resolution=0.2), resulting in Newborn CPN and CFuPN. Cluster 6 was further subclustered (resolution=0.2), resulting in aRG and unspecified PN. Cluster 16 was further subclustered (resolution=0.1), resulting in oRG and IP. 
#For the GM08330 b2 3 months data set, cluster 6 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 8 was further subclustered (resolution=0.2), resulting in aRG, oRG and unspecified PN. #For the Mito 210 c1 b6 3 months data set, cluster 14 was further subclustered (resolution=0.2), resulting in aRG and oRG.
#For the Mito 210 c2 b5 3 months data set, cluster 6 was further subclustered (resolution=0.2), resulting in aRG and oRG. Cluster 9 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 16 was further subclustered (resolution=0.5), resulting in aRG and oRG. 
#For the HUES66 b7 3 months data set, cluster 7 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 12 was further subclustered (resolution=0.2), resulting in aRG and oRG. Cluster 16 was further subclustered (resolution=0.1), resulting in aRG and oRG. 
#For the HUES66 b8 3 months data set, cluster 9 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 16 was further subclustered (resolution=0.2), resulting in aRG and oRG. 
#For the PGP1 b9 3 months data set, cluster 10 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 13 was further subclustered (resolution=0.5), resulting in aRG and oRG.
#For the PGP1 b10 3 months data set, cluster 8 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 11 was further subclustered (resolution=0.2), resulting in aRG and oRG. 
#For the Mito 210 c1 b1 4 months data set, cluster 3 was further subclustered (resolution=0.5), resulting in Unspecified PN and CPN. Cluster 12 was further subclustered (resolution=0.1), resulting in IP and oRG. Cluster 13 was further subclustered (resolution=0.5), resulting in IP and oRG.
#For the Mito 210 c1 b1 5 months data set, cluster 5 was further subclustered (resolution=0.5), resulting in IP, oRG and IN progenitors. Cluster 7 was further subclustered (resolution=0.5), resulting in IP, oRG, Immature IN and Glial precursors. Cluster 18 was further subclustered (resolution=0.5), resulting in IP and oRG.
#For the 11a b11 6 months data set, cluster 3 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 11 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors.
#For the GM08330 b12 6 months data set, cluster 6 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 9 was further subclustered (resolution=0.2), resulting oRG and IN progenitors. Cluster 11 was further subclustered (resolution=0.5), resulting in Glial precursors and Immature IN. 
#For the Mito 210 c1 b3 6 months data set, cluster 5 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. 
#For the Mito 210 c2 b13 6 months data set, cluster 1 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 4 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 13 was further subclustered (resolution=0.5), resulting in IN and CPN. 
#For the HUES66 b8 6 months data set, cluster 4 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 8 was further subclustered (resolution=0.2), resulting in aRG and Immature IN. Cluster 11 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. 
#For the PGP1 b9 6 months data set, cluster 13 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 14 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. 
#For the PGP1 b14 6 months dataset, cluster 1 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 5 was further subclustered (resolution=0.2), resulting in Glial precursors and Immature IN. 

#For the SHARE-seq RNA data, each time point was isolated and subclustered with resolution=1.0, except day_59, where resolution=1.8
#For day_90, clusters 7, 8, 9, were isolated and further subclustered (resolution=1.0), resulting in aRG and IP.