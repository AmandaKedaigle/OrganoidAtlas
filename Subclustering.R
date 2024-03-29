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
obj1 = FindNeighbors(obj1, dims=c(1:20)) #15 or 20
obj1 = FindClusters(obj1, resolution=c(1,0.5, 0.2, 0.1))
obj1 <- RunTSNE(object = obj1, dims = c(1:20)) #15 or 20
obj1 <- RunUMAP(obj1, dims = 1:20) #15 or 20

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

#Dorsal organoids
#For the Mito 210 c1 b1 23 days data set, clusters 16, 6, and 4 were further subclustered (resolution=0.5), resulting in Preplate/Subplate cells and FOXG1- EMX1- neurons. 
#For the PGP1 c1 b2 23 days data set, cluster 10 was further subclustered (resolution=0.2), resulting in Preplate/Subplate and FOXG1-EMX1- neurons. Cluster 11 was further subclustered (resolution=0.5), resulting in Subcortical neuronal precursors, Subcortical neurons, Subcortical interneurons. 
#For the GM08330 b3 1 month data set, cluster 10 was further subclustered (resolution=0.5), resulting in Subcortical neurons and Cajal Retzius cells. 
#For the Mito 210 c1 b5 1 month data set, cluster 10 was further subclustered (resolution=0.1), resulting in Newborn DL PN and Subcortical neurons. 
#For the Mito 210 c2 b6 1 month data set, cluster 6 was further subclustered (resolution=0.2), resulting in Subcortical neurons and Cajal Retzius cells. 
#For the Mito 210 c1 b7 1.5 months data set, cluster 10 was further subclustered (resolution=0.2), resulting in Newborn DL PN and Immature DL PN. Cluster 12 was further subclustered (resolution=0.5) resulting in Subcortical neurons and Cajal Retzius. 
#For the PGP1 c1 b2 1.5 months data set, cluster 7 was further subclustered (resolution=0.2), resulting in Newborn DL PN and Immature DL PN. Cluster 8 was further subclustered (resolution=0.2), resulting in IP and Newborn DL PN.
#For the Mito 210 c1 b1 2 months data set, cluster 1 was further subclustered (resolution=0.5), resulting in aRG and oRG. Cluster 5 was further subclustered (resolution=0.2), resulting in Newborn CPN and CFuPN. Cluster 6 was further subclustered (resolution=0.2), resulting in aRG and unspecified PN. Cluster 16 was further subclustered (resolution=0.1), resulting in oRG and IP. 
#For the PGP1 c1 b2 2 months data set, cluster 0 was further subclustered (resolution=0.5), resulting in aRG and oRG. Cluster 1 was further subclustered (resolution=0.5), resulting in CFuPN and Newborn CPN. Cluster 4 was further subclustered (resolution=0.2), resulting in CFuPN and unspecified PN. Cluster 5 was further subclustered (resolution=0.1), resulting in Cortical hem and aRG. Cluster 6 was further subclustered (resolution=0.5), resulting in IP and unspecified PN. Cluster 12 was further subclustered (resolution=0.1), resulting in IP and oRG. 
#For the GM08330 b3 3 months data set, cluster 6 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 8 was further subclustered (resolution=0.2), resulting in aRG, oRG and unspecified PN. #For the Mito 210 c1 b6 3 months data set, cluster 14 was further subclustered (resolution=0.2), resulting in aRG and oRG.
#For the Mito 210 c2 b6 3 months data set, cluster 6 was further subclustered (resolution=0.2), resulting in aRG and oRG. Cluster 9 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 16 was further subclustered (resolution=0.5), resulting in aRG and oRG. 
#For the HUES66 b9 3 months data set, cluster 7 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 12 was further subclustered (resolution=0.2), resulting in aRG and oRG. Cluster 16 was further subclustered (resolution=0.1), resulting in aRG and oRG. #For the HUES66 b10 3 months data set, cluster 9 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 16 was further subclustered (resolution=0.2), resulting in aRG and oRG. 
#For the PGP1 c2 b11 3 months data set, cluster 10 was further subclustered (resolution=0.1), resulting in oRG and IP. Cluster 13 was further subclustered (resolution=0.5), resulting in aRG and oRG.
#For the PGP1 c2 b12 3 months data set, cluster 8 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 11 was further subclustered (resolution=0.2), resulting in aRG and oRG. 
#For the Mito 210 c1 b1 4 months data set, cluster 3 was further subclustered (resolution=0.5), resulting in Unspecified PN and CPN. Cluster 12 was further subclustered (resolution=0.1), resulting in IP and oRG. Cluster 13 was further subclustered (resolution=0.5), resulting in IP and oRG.
#For the PGP1 c1 b13 4 months data set, cluster 3 was further subclustered (resolution=0.5), resulting in unspecified PN and CPN. Cluster 12 and 17 were further subclustered (resolution=0.2), resulting in IP and oRG.
#For the Mito 210 c1 b1 5 months data set, cluster 5 was further subclustered (resolution=0.5), resulting in IP, oRG and IN progenitors. Cluster 7 was further subclustered (resolution=0.5), resulting in IP, oRG, Immature IN and Glial precursors. Cluster 18 was further subclustered (resolution=0.5), resulting in IP and oRG.
#For the PGP1 c1 b13 5 months data set, cluster 6 was further subclusted (resolution=0.2), resulting in IP, oRG and Glial precursors. Cluster 10 was further subclustered (resolution=0.5), resulting in oRG, Glial precursors and Immature IN. Cluster 12 was further subclustered (resolution=0.2), resulting in aRG and oRG. Cluster 15 was further subclustered (resolution=0.2),resulting in oRG/Astroglia and unspecified PN.
#For the 11a b14 6 months data set, cluster 3 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 11 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors.
#For the GM08330 b15 6 months data set, cluster 6 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 9 was further subclustered (resolution=0.2), resulting oRG and IN progenitors. Cluster 11 was further subclustered (resolution=0.5), resulting in Glial precursors and Immature IN. 
#For the Mito 210 c1 b4 6 months data set, cluster 5 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. 
#For the Mito 210 c2 b16 6 months data set, cluster 1 was further subclustered (resolution=0.2), resulting in oRG and IP. Cluster 4 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 13 was further subclustered (resolution=0.5), resulting in IN and CPN. 
#For the HUES66 b10 6 months data set, cluster 4 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 8 was further subclustered (resolution=0.2), resulting in aRG and Immature IN. Cluster 11 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. 
#For the PGP1 c2 b11 6 months data set, cluster 13 was further subclustered (resolution=0.1), resulting in oRG and IN progenitors. Cluster 14 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. 
#For the PGP1 c2 b17 6 months dataset, cluster 1 was further subclustered (resolution=0.2), resulting in oRG and IN progenitors. Cluster 5 was further subclustered (resolution=0.2), resulting in Glial precursors and Immature IN. 

#Whole brain organoids
#For the GM08330 3 months data set, cluster 7 was further subclustered (resolution=0.1), resulting in RG and IP.
#For the HUES66 3 months data set, cluster 7 was further subclustered (resolution=0.2), resulting in RG and IP. Cluster 10 was further subclustered (resolution=0.1), resulting in retinal ganglion cells and amacrine cells. 

#Fetal data
#Cluster 0 was further subclustered (resolution=0.2), resulting in IP and CPN. Cluster 22 was further subclustered (resolution=0.1), resulting in IN progenitors and Migrating IN. Clusters 11, 13 and 30 were subclustered together (resolution=0.2), resulting in aRG and oRG.

#For the SHARE-seq RNA data, each time point was isolated and subclustered with resolution=1.0, except day_59, where resolution=1.8
#For day_90, clusters 7, 8, 9, were isolated and further subclustered (resolution=1.0), resulting in aRG and IP.