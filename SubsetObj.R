#subsets the object in neurons vs progenitors
#creates metadata column for DESeq2 analysis, to compare each population to the other. E.g. CFuPN to other excitatory neurons (PN, CPN)

Idents(obj) <- "FinalName"
DimPlot(obj)

neurons <- subset(obj, idents = c("PN", "CFuPN", "CPN"))
Idents(neurons) <- "FinalName"
DimPlot(neurons)

saveRDS(neurons, file = "neurons.rds")

#To Compare CPN to the rest, ie make CPN and other neurons. Name all the other neurons 'other'
#To save column, make a new column
neurons@meta.data$CPNAnalysis <- neurons@meta.data$FinalName # create a column like the one already generated without modifying the original
neurons$CPNAnalysis <- as.vector(neurons$CPNAnalysis)
neurons$CPNAnalysis [neurons$CPNAnalysis=="PN"] = "Other"
neurons$CPNAnalysis [neurons$CPNAnalysis=="CFuPN"] = "Other"
Idents(neurons) <- "CPNAnalysis"
DimPlot(neurons)

## find markers
markers_neurons_CPN <- FindMarkers(neurons, ident.1 = 'CPN', ident.2 = 'Other', min.pct = 0.25, only.pos = F) 
write.csv(markers_neurons_CPN, file = "markers_neurons_CPN.csv")

##Compare CFuPN to the rest

neurons@meta.data$CFuPNAnalysis <- neurons@meta.data$FinalName # create a column like the one already generated without modifying the original
neurons$CFuPNAnalysis<- as.vector(neurons$CFuPNAnalysis)
neurons$CFuPNAnalysis [neurons$CFuPNAnalysis=="CPN"] = "Other"
neurons$CFuPNAnalysis [neurons$CFuPNAnalysis=="PN"] = "Other"
Idents(neurons) <- "CFuPNAnalysis"
DimPlot(neurons)

## find markers
markers_neurons_CFuPN <- FindMarkers(neurons, ident.1 = "CFuPN", ident.2 = "Other", min.pct = 0.25, only.pos = F) 
write.csv(markers_neurons_CFuPN, file = "markers_neurons_CFuPN.csv")

##Compare PN to the rest 

neurons@meta.data$PNAnalysis <- neurons@meta.data$FinalName # create a column like the one already generated without modifying the original
neurons$PNAnalysis<- as.vector(neurons$PNAnalysis)
neurons$PNAnalysis [neurons$PNAnalysis=="CPN"] = "Other"
neurons$PNAnalysis [neurons$PNAnalysis=="CFuPN"] = "Other"
Idents(neurons) <- "PNAnalysis"
DimPlot(neurons)

## find markers
markers_neurons_PN <- FindMarkers(neurons, ident.1 = "PN", ident.2 = "Other", min.pct = 0.25, only.pos = F) 
write.csv(markers_neurons_PN, file = "markers_neurons_PN.csv")

#save the neurons object with the new columns
saveRDS(neurons, file = "neurons.rds")

## for progenitors

Idents(obj) <- "FinalName"
DimPlot(obj)
progs <- subset(obj, idents = c("aRG", "IP", "oRG"))
Idents(progs) <- "FinalName"
DimPlot(progs)

saveRDS(progs, file = "progs.rds")

#Compare aRG to the rest, ie make aRG and other progs. Name all the other progs 'other'
#To save column, make a new column
progs@meta.data$aRGAnalysis <- progs@meta.data$FinalName # create a column like the one already generated without modifying the original
progs$aRGAnalysis <- as.vector(progs$aRGAnalysis)
progs$aRGAnalysis [progs$aRGAnalysis=="IP"] = "Other"
progs$aRGAnalysis [progs$aRGAnalysis=="oRG"] = "Other"
Idents(progs) <- "aRGAnalysis"
DimPlot(progs)

## find markers
markers_progs_aRG <- FindMarkers(progs, ident.1 = 'aRG', ident.2 = 'Other', min.pct = 0.25, only.pos = F) 
write.csv(markers_progs_aRG, file = "markers_progs_aRG.csv")

#Compare IP to the rest

progs@meta.data$IPAnalysis <- progs@meta.data$FinalName # create a column like the one already generated without modifying the original
progs$IPAnalysis <- as.vector(progs$IPAnalysis)
progs$IPAnalysis [progs$IPAnalysis=="aRG"] = "Other"
progs$IPAnalysis [progs$IPAnalysis=="oRG"] = "Other"
Idents(progs) <- "IPAnalysis"
DimPlot(progs)

## find markers
markers_progs_IP <- FindMarkers(progs, ident.1 = 'IP', ident.2 = 'Other', min.pct = 0.25, only.pos = F) 
write.csv(markers_progs_IP, file = "markers_progs_IP.csv")

#Compare oRG to the rest

progs@meta.data$oRGAnalysis <- progs@meta.data$FinalName # create a column like the one already generated without modifying the original
progs$oRGAnalysis <- as.vector(progs$oRGAnalysis)
progs$oRGAnalysis [progs$oRGAnalysis=="aRG"] = "Other"
progs$oRGAnalysis [progs$oRGAnalysis=="IP"] = "Other"
Idents(progs) <- "oRGAnalysis"
DimPlot(progs)

## find markers
markers_progs_oRG <- FindMarkers(progs, ident.1 = 'oRG', ident.2 = 'Other', min.pct = 0.25, only.pos = F) 
write.csv(markers_progs_oRG, file = "markers_progs_oRG.csv")

#save the progenitors object with the new columns
saveRDS(progs, file = "progs.rds")




