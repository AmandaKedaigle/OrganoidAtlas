#Takes multiple Seurat Objects and performs differential expression analysis 
#Between *the seurat objects*
#By splitting the objects into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise
#Runs DESeq2 using a design formula which batches samples based on which batch the cells came from

library(Seurat)
library(DESeq2)
library(writexl)

#Load this function as-is for use later
combineDEbatch<-function(seurList,condition="age",base="Other",batch="dataset",combineOn="organoid", #These are defaults, but will be overrided by what you put below
                        #Below are default settings you can change if you know what you're doing
                        minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                        minSamples=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                        minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                        genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                        form="" #design formula for DESeq2 if you want to include more variables in addition to "~ batch + condition"
)
{
  data.all = list()
  for (i in 1:length(seurList)) { 
   Idents(seurList[[i]])=combineOn
    
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
  val=do.call(rbind, lapply(seurList, function(x) x@meta.data[,c(condition,combineOn,batch,extraColumns)]))
  val=val[!duplicated(val[,2]),]
  rownames(val)=val[,2]
  keepGroup=as.character(val[keepOrgs,1]) #age for each org
  levels = levels(factor(keepGroup))
  if(length(levels)<2) {
    #less than 2 ages left after minCells
    print("Not enough batches per treatment group with minimum # of cells!")
    return(NULL)
  }
  for (level in levels) {
    if(sum(keepGroup==level)<minSamples) {
      #less than minSamples organoids left in this age
      print("Not enough batches per treatment group with minimum # of cells!")
      return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepGroup)
  if (base != "") { colDat = relevel(colDat, ref=base)}
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[keepOrgs,"sample"] = factor(val[keepOrgs,2])
  colDat[keepOrgs,"batch"] = factor(val[keepOrgs,3])
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Run DESeq2")
  design= ~  condition #This is different! I am saying batch instead of sample - sample was for pairing CellTypes within the same organoid, here we group by batches
  if(nchar(form)>0){
    design=as.formula(paste("~",form,"batch + condition",sep=""))
  }
  print(design)
  
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  dds <- DESeq(dds)
  saveRDS(dds, "dds.rds")
  out=data.frame(results(dds))
  out=out[order(out$pvalue),]
  saveRDS(out, "res.rds")
  return(out)
}

#Load Seurat Objects
seur1 = readRDS("obj1.rds")
seur2 = readRDS ("obj2.rds")
seur3 = readRDS("obj3.rds")
seurList = c(seur1, seur2, seur3)

#Each seur should have a metadata column for organoid that seperates organoids from each dataset (so like GM_1m_org1 is different from Mito_1m_org1)
#and a metadata column that just has dataset (so GM_1m vs Mito_1m)
#and a metadata column that has age that you want calculate for (so like "1m" and the others say "Other")

#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
#Make sure all of the seurat objects you are using have metadata columns using the same name for this metadata column and values
#i.e. if in two seurat objects, CPNs are called cluster "1" in one object and cluster "7" in another, this code won't combine them for you - you'd have to add new metadata to those objects that matches
condition = "age" 
clusterOfInterest = "1m"
base = "Other"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "org"

#This will run DE analysis and save a .xlsx file!
degs <- combineDEbatch(seurList, condition=condition, base=base, combineOn=combineOn)
degs$gene = rownames(degs)
if (length(degs)>0) {
  write_xlsx(degs, path=paste0(clusterOfInterest,".vs.",base,".DEGs.xlsx"), format_headers = T)
}
