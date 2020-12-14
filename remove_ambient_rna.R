library(foreach)
library(doParallel)
library(Seurat)
library(stringr)
options(stringsAsFactors=F)
source("RunDecontX.R")


# input file
fname = commandArgs(trailingOnly=T)
seur_path = "final_objects/"
prefix = str_split(fname, "_", simplify=T)[,1]
cat("object path:\n")
cat(paste0(seur_path, fname, '\n'))

obj.big <- readRDS(paste0(seur_path, fname))
orgs <- names(table(obj.big$org))

# parallelize
numCores <- detectCores()
cat("numCores\n")
print(numCores)
registerDoParallel(numCores)

nOrgs = length(table(obj.big$org))

foreach (i = seq(1,nOrgs,3)) %dopar% {
# file too large to convert sparse matrix directly to matrix
if (nOrgs==20 && i==19) {
cur.prefix = paste0(prefix, "_org", i, "-", i+1)
print(cur.prefix)
obj <- subset(obj.big, org %in% orgs[i:(i+1)])
} else {
cur.prefix = paste0(prefix, "_org", i, "-", i+2)
print(cur.prefix)
obj <- subset(obj.big, org %in% orgs[i:(i+2)])
}
RunDecontX(obj, cur.prefix)
}

print("DONE")
