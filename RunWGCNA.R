library(WGCNA)
library(Seurat)





RunWGCNA<-function(seur,output,myPower,genesToUse=NULL, lowCutoff=.15, highCutoff=9)
{
    print("gene expression cutoffs:")
    print(paste0("lowCutoff=", lowCutoff))
    print(paste0("highCutoff=", highCutoff))

    print("Filter Data!")
    dat=seur@assays$RNA@data

    if (is.null(genesToUse)) {
    dat<-dat[!(Matrix::rowMeans(dat, na.rm = TRUE) < lowCutoff),] #get rid of lowly expressed genes
    dat<-dat[!(Matrix::rowMeans(dat, na.rm = TRUE) > highCutoff),] # get rid of highly expressed genes
    dat <- dat[!(rownames(dat) %in% grep(pattern = "^RP", x = rownames(dat), value = TRUE)),] ##Remove ribosomal genes
    dat <- dat[!(rownames(dat) %in% grep(pattern = "^MT-", x = rownames(dat), value = TRUE)),] ##Remove mito genes
    } else {
    sprintf("Using custom DEGs of length: %d", length(genesToUse))
    dat <- dat[genesToUse,]
    }

    dat = as.matrix(dat)
    
    print("Remaining data:")
    print(dim(dat))
    print("Columns")
    print(head(colnames(dat)))
    print("Rows")
    print(head(rownames(dat)))
    
    if (!(file.exists("datExpr.rds"))) {
    print("Saving expr data matrix")
    saveRDS(dat, "datExpr.rds")
    }
    
    print("Transpose expr matrix")
    dat=t(dat)
    
    # Module parameters
    networkType = "signed"
    minModuleSize = 4  # 7
    mergeCutHeight = 0.15

    print("Get modules!")
    net = blockwiseModules(
    # Input data
    dat,
    
    # Data checking options
    checkMissingData = TRUE,
    
    # Options for splitting data into blocks
    blocks = NULL,
    maxBlockSize = 5000,  
    randomSeed = 59069,
    
    # load TOM from previously saved file?
    loadTOM = FALSE,
    
    # Network construction arguments: correlation options
    corType = "bicor", # more robust for non-normal? Song 2012
    maxPOutliers = 0.1,  
    quickCor = 0,
    pearsonFallback = "individual",
    cosineCorrelation = FALSE,
    
    # Adjacency function options
    # this is where power gets input
    power = myPower,
    networkType = networkType,
    
    # Topological overlap options
    TOMType = "signed",
    TOMDenom = "min",
    
    # Saving or returning TOM
    getTOMs = NULL,
    saveTOMs = TRUE,
    saveTOMFileBase = paste0("blockwiseTOM_power", myPower),
    
    # Basic tree cut options
    deepSplit = 3,  
    detectCutHeight = 0.995,
    minModuleSize = minModuleSize,
    
    # Advanced tree cut options
    maxCoreScatter = NULL, minGap = NULL,
    maxAbsCoreScatter = NULL, minAbsGap = NULL,
    minSplitHeight = NULL, minAbsSplitHeight = NULL,
    useBranchEigennodeDissim = FALSE,
    minBranchEigennodeDissim = mergeCutHeight,
    pamStage = TRUE, pamRespectsDendro = TRUE,
    
    # Gene reassignment, module trimming, and module "significance" criteria
    reassignThreshold = 1e-6,
    minCoreKME = 0.5,
    minCoreKMESize = minModuleSize/3,
    minKMEtoStay = 0.3,
    
    # Module merging options
    mergeCutHeight = 0.15,
    impute = TRUE,
    trapErrors = TRUE,
    
    # Output options
    numericLabels = FALSE,
    
    # Options controlling behaviour
    nThreads = 6,
    verbose = 3, 
    indent = 0)
    
    mods=NULL
    
    print("Save modules!")
    saveRDS(net,file=output)
    
    print("Done!")
    
    return(net)
    
}


