#My token for FLOWMAPR: fc520579ab3596e6d11bfed367d6504349347fa5
library(FLOWMAPR)
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
library(RColorBrewer)

setwd("~/Documents/wt_merge/FlowMAP/allTimepts-224/")

seur23 = readRDS("../../FinalObjects/Dorsal/dorsal_23d_111520.rds")
seur23$time = "23d"
seur1 = readRDS("../../FinalObjects/Dorsal/dorsal_1moharm_111520_final.rds")
seur1$time = "1m"
seur2 = readRDS("../../FinalObjects/Dorsal/dorsal_2mo_111520_final.rds")
seur2$time = "2m"
seur3 = readRDS("../../FinalObjects/Dorsal/dorsal_3mo_111520_final.rds")
seur3$time = "3m"
seur4 = readRDS("../../FinalObjects/Mito210c1_4mo_022321.rds")
seur4$time = "4m"
seur5 = readRDS("../../FinalObjects/Mito210c1_5mo_022321.rds")
seur5$time = "5m"
seur6 = readRDS("../../FinalObjects/6mo_harmonizedObj_103020.rds")
seur6$time = "6m"

seuset = merge(seur23, c(seur1, seur2, seur3, seur4, seur5, seur6), add.cell.ids = c("23d", 1:6), merge.data=T)
seuset = FindVariableFeatures(seuset, selection.method='mean.var.plot')
seuset<-ScaleData(seuset,features=VariableFeatures(seuset),vars.to.regress=c("nCount_RNA","CC.Difference"))
seuset = RunPCA(seuset)
Idents(seuset) = "time"
saveRDS(seuset, "combined.rds")
#seuset = subset(seuset, downsample=3000) #initial trials with downsampled data, hopefully not necessary

#Take PCA score, first 5 dimensions based on elbow plot
pca.values <- FetchData(seuset, vars = paste0("PC_", 1:15), slot="data")

## Get metadata
celltypes = as.integer(factor(seuset$FinalName))
celltypesMap = levels(factor(seuset$FinalName))
names(celltypesMap) = 1:length(celltypesMap)
saveRDS(celltypesMap, "celltypeMap.rds")
times = as.integer(factor(seuset$time, levels = c("23d", "1m","2m","3m", "4m", "5m", "6m")))

#Concatenate variables for visualization on FLOWMAP
df <- as.data.frame(cbind(t(GetAssayData(seuset, slot = "scale.data")),celltypes,pca.values,times))
dfl = split(df, df$times)

#rm(seuset)

## FLOWMAP
#Parameters
mode <- "single"
project.name <- "Dorsal Organoid scRNA"
save.folder <- "."
clustering.var <- colnames(pca.values)
distance.metric <- "euclidean"
#per <- 1
minimum <- 2
maximum <- 5
seed.X <- 1
set.seed(seed.X)
name.sort <- FALSE
savePDFs <- TRUE
which.palette <- "bluered"
time.col.label <- "times"
condition.col.label <- NULL
#clustering <- FALSE #can get away with doing indivudal cells up to ~12k cells
clustering <- TRUE
cluster.numbers = 3000
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
