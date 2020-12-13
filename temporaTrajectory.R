library(Tempora)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(igraph)

#Load and cluster graph from FLOWMAP
g = read_graph("../2020-12-12_Dorsal Organoid scRNA_FLOW-MAP_xy_orig_time_22.06.07.graphml", format = "graphml")
c = cluster_louvain(g) # If the graph has a weight edge attribute, then this is used by default
g = set_vertex_attr(g, "community", value=membership(c))
write.table(cbind(V(g)$id, V(g)$community), "vertex_communities.txt", sep="\t",row.names = F, quote=F)
saveRDS(g, "igraphObj.rds")

attrs = as_data_frame(g, what="vertices")
meta = attrs[,c(paste0("PC_", 1:10), "celltypes", "Timepoint", "percent.total", "name", "size", "x","y","id", "community")]
exprMat = as.matrix(attrs[,!colnames(attrs) %in% colnames(meta)])
colnames(meta)[colnames(meta)=="Timepoint"] = "Timepoints"
colnames(meta)[colnames(meta)=="community"] = "Clusters"
meta$Timepoints[meta$Timepoints==1] = "23d"
meta$Timepoints[meta$Timepoints==2] = "1m"
meta$Timepoints[meta$Timepoints==3] = "2m"
meta$Timepoints[meta$Timepoints==4] = "3m"

#Calling clusters based on majority cell type
t = table(meta$Clusters, meta$celltypes)
calls = as.numeric(colnames(t)[apply(t, which.max, MARGIN = 1)])
typesMap = readRDS("../celltypesMap.rds")
types = typesMap[calls]

#tempura ----

cortex_tempora <- CreateTemporaObject(t(exprMat), meta,
                                      timepoint_order = c("23d", "1m", "2m", "3m"),
                                      cluster_labels=types)

#Estimate pathway enrichment profiles of clusters
cortex_tempora <- CalculatePWProfiles(cortex_tempora, 
                                      gmt_path = "Human_GOBP_AllPathways_no_GO_iea_November_17_2020_symbol.gmt",
                                      method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)

#Build trajectory with 6 PCs 
cortex_tempora <- BuildTrajectory(cortex_tempora, n_pcs = 10, difference_threshold = 0.01)

#Fit GAMs on pathway enrichment profile
cortex_temporaP <- IdentifyVaryingPWs(cortex_tempora, pval_threshold = 0.1)

#Plot expression trends of significant time-varying pathways
pdf("Pathways.pdf")
PlotVaryingPWs(cortex_temporaP)
dev.off()

#Visualize the trajectory
edge_graph <- graph_from_data_frame(d=cortex_tempora@trajectory, vertices = cortex_tempora@cluster.metadata, directed = T)
write_graph(edge_graph, "temporaGraph.graphml", format="graphml")
#loading that into cytoscape