#Takes multiple Seurat Objects and performs differential expression analysis 
#Between two clusters/celltypes
#By splitting the objects into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise
#Runs DESeq2 using a design formula which pairs samples based on which sample (i.e. organoid) the cells came from

library(Seurat)
library(DESeq2)
library(writexl)


#Load this function as-is for use later
combineDEpair<-function(seurList,id,condition="CellType",base="Cycling",batch="dataset",combineOn="orig.ident", #These are defaults, but will be overrided by what you put below
                        #Below are default settings you can change if you know what you're doing
                        minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                        minBatches=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                        minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                        genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                        form="" #design formula for DESeq2 if you want to include more variables in addition to "~ + sample + condition"
)
{
  data.all = list()
  for (i in 1:length(seurList)) { 
    Idents(seurList[[i]]) = condition 
    print(paste("Subsample",i))
    seurList[[i]] = subset(seurList[[i]],idents=c(id))
  
    #in this setting, we have cells in 1 organoid that fall into several conditions (cell types)
    #we'll need a column in the metadata to distinguish this
    seurList[[i]]$org.condition = paste(seurList[[i]]@meta.data[,combineOn], seurList[[i]]@meta.data[,condition], sep=".")
  
    Idents(seurList[[i]])="org.condition"
  
    genes.use=rownames(GetAssayData(seurList[[i]],slot="counts"))
    if(length(genes)>0){genes.use=genes}
  
    print(paste("Combine data for sample",i))
    data.all[[i]]=data.frame(row.names = genes.use)
    for(t in levels(Idents(seurList[[i]]))) {
      temp.cells=WhichCells(seurList[[i]],ident=t)
      if (length(temp.cells)==1) data.temp=(GetAssayData(seurList[[i]],slot="counts")[genes.use,temp.cells])
      if (length(temp.cells)>1) data.temp=apply(GetAssayData(seurList[[i]],slot="counts")[genes.use,temp.cells],1,sum)
      data.all[[i]]=cbind(data.all[[i]],data.temp)
      colnames(data.all[[i]])[ncol(data.all[[i]])]=t
    }
  }
  print("Combine data together")
  data.all = do.call(cbind, data.all)
  
  print("Filter samples for minimum cells")
  summary = do.call(c,lapply(seurList, function(x) summary(Idents(x))))
  keepOrgs=names(summary)[summary>minCells]
  numOrg=length(keepOrgs)
  print(paste("Keeping", numOrg, "samples"))
  data.all=data.all[,keepOrgs]
  
  extraColumns<-strsplit(form,"+",fixed=T)[[1]]
  val=do.call(rbind, lapply(seurList, function(x) x@meta.data[,c(condition,combineOn,"org.condition",batch,extraColumns)]))
  val=val[!duplicated(val[,3]),]
  rownames(val)=val[,3]
  keepBatch=as.character(val[keepOrgs,1])
  levels = levels(factor(keepBatch))
  if(length(levels)<2) {
    print("Not enough batches per treatment group with minimum # of cells!")
    return(NULL)
  }
  for (level in levels) {
    if(sum(keepBatch==level)<2) {
      print("Not enough batches per treatment group with minimum # of cells!")
      return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepBatch)
  if (base != "") { colDat = relevel(colDat, ref=base)}
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[keepOrgs,"sample"] = factor(val[keepOrgs,2])
  colDat[keepOrgs,"batch"] = factor(val[keepOrgs,4])
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Run DESeq2")
  design= ~ sample + condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form,"sample + condition",sep=""))
  }
  print(design)
  
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  dds <- DESeq(dds)
  out=data.frame(results(dds))
  out=out[order(out$pvalue),]
  return(out)
}

#Load Seurat Objects - if you have multiple individual objects
seur = readRDS("v3wt.DorsalKadoshima.HUES66_3mon.seur.rds")
seur2 = readRDS("v3wt.DorsalKadoshima.GM_3mon.seur.rds")
seurList = c(seur, seur2)

#Alternative to above - if you have one object and want to split it on a metadata column into a seurList
#change "dataset" in two lines below to the name of the column you want to split on!
seur = readRDS("3mo_harmonizedObj_102820.rds")
levels = levels(factor(seur$dataset)) #gets a list of all datasets in the object
seurList = lapply(levels, function(x) subset(seur, subset=dataset==x))

#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
#Make sure all of the seurat objects you are using have metadata columns using the same name for this metadata column and values
#i.e. if in two seurat objects, CPNs are called cluster "1" in one object and cluster "7" in another, this code won't combine them for you - you'd have to add new metadata to those objects that matches
condition = "CellType"
clusterOfInterest = "CPNs"
base = "CFuPNs"

#Set batch to the metadata column that distinguishes the seurat objects
batch = "dataset"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "org"

#This will run DE analysis and save a .xlsx file!
id = c(clusterOfInterest,base)
degs <- combineDEpair(seurList, id=id, condition=condition, base=base, batch=batch, combineOn=combineOn)
degs$gene = rownames(degs)
if (length(degs)>0) {
  write_xlsx(degs, path=paste0(clusterOfInterest,".vs.",base,".DEGs-H66andGM.xlsx"), format_headers = T)
}






