library(FLOWMAPR)
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
library(RColorBrewer)

seuset = readRDS("../combined.rds") #merged (non-batch-corrected) object of dorsal cells from 23days-3months. PCA has been run

#Take PCA score, first 10 dimensions based on elbow plot
pca.values <- FetchData(seuset, vars = paste0("PC_", 1:10), slot="data")

## Get metadata
celltypes = as.integer(factor(seuset$FinalName))
celltypesMap = levels(factor(seuset$FinalName)) #to go back to cell type labels from integers for labeling
times = as.integer(factor(seuset$time, levels = c("23d", "1m","2m","3m")))

#Concatenate variables for visualization on FLOWMAP
df <- as.data.frame(cbind(t(GetAssayData(seuset, slot = "scale.data")),celltypes,pca.values,times))
dfl = split(df, df$times)

## FLOWMAP
#Parameters
mode <- "single"
project.name <- "Dorsal Organoid scRNA"
save.folder <- "."
clustering.var <- colnames(pca.values) #clustering based on PCA
distance.metric <- "euclidean"
minimum <- 2
maximum <- 5
seed.X <- 1
set.seed(seed.X)
name.sort <- FALSE
savePDFs <- FALSE
which.palette <- "bluered"
time.col.label <- "times"
condition.col.label <- NULL
clustering <- TRUE
cluster.numbers = 3000 #Clustering each time point to 3k representative nodes
cluster.mode = "kmeans"

#Run FLOWMAP
FLOWMAPR::FLOWMAPfromDF(mode = mode, df = dfl, project.name = project.name,
                        time.col.label = time.col.label, condition.col.label = condition.col.label,
                        clustering.var = clustering.var, distance.metric = distance.metric,
                        minimum = minimum, maximum = maximum,
                        save.folder = save.folder,
                        name.sort = name.sort, clustering = clustering, 
                        cluster.numbers = cluster.numbers, cluster.mode = cluster.mode,
                        seed.X = seed.X, savePDFs = savePDFs, which.palette = which.palette)
