#Takes a Seurat Object and processes the same way we did for differential expression
#by summing across each organoid
#and makes plots for a list of example genes

library(Seurat)
library(DESeq2)
library(ggplot2)

#Load this function as-is for use later
DEfigure <-function(seur,id,condition="treat",base="wt",combineOn="orig.ident", #These are defaults, but will be overrided by what you put below
                    colorBy = "",
                    #Below are default settings you can change if you know what you're doing
                    minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                    minBatches=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                    minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                    genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                    form="" #design formula for DESeq2 if you want to include more variables in addition to "~ sample + condition"
)
{
  print("Subsample")
  seur<-subset(seur,idents=c(id))
  
  #in this setting, we have cells in 1 organoid that fall into several conditions (cell types)
  #we'll need a column in the metadata to distinguish this
  seur$org.condition = paste(seur@meta.data[,combineOn], seur@meta.data[,condition], sep=".")
  
  Idents(seur)="org.condition"
  genes.use=rownames(GetAssayData(seur,slot="counts"))
  if(length(genes)>0){genes.use=genes}
  
  print("Combine data per sample")
  data.all=data.frame(row.names = genes.use)
  for(i in levels(Idents(seur))) {
    temp.cells=WhichCells(seur,ident=i)
    if (length(temp.cells)==1) data.temp=(GetAssayData(seur,slot="counts")[genes.use,temp.cells])
    if (length(temp.cells)>1) data.temp=apply(GetAssayData(seur,slot="counts")[genes.use,temp.cells],1,sum)
    data.all=cbind(data.all,data.temp)
    colnames(data.all)[ncol(data.all)]=i
  }
  
  print("Filter samples for minimum cells")
  keepOrgs=names(summary(Idents(seur)))[summary(Idents(seur))>minCells]
  numOrg=length(keepOrgs)
  print(paste("Keeping", numOrg, "samples"))
  data.all=data.all[,keepOrgs]
  
  extraColumns<-strsplit(form,"+",fixed=T)[[1]]
  extraColumns = append(extraColumns, colorBy)
  val=seur@meta.data[,c(condition,combineOn,"org.condition",extraColumns)]
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
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  colnames(colDat)[ncol(colDat)] = "colorBy"
  
  print("Set Up DESeq2")
  design= ~ sample + condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form," + sample + condition",sep=""))
  }
  print(design)
  
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  
  return(dds)
}

#Load Seurat Object
seur = readRDS("obj.rds")
#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
condition = "CPNAnalysis"
clusterOfInterest = "CPN"
base = "Other"

#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "org"
#set colorBy to the metadata column you want each sample to be colored by in the final figure
colorBy = "dataset"
seur@meta.data[,colorBy] = factor(seur@meta.data[,colorBy]) #make it a factor so it is not 1-4 gradient

#This will run DE analysis and save a .xlsx file!
Idents(seur) = condition
id = c(clusterOfInterest,base)
dds <- DEfigure(seur, id=id, condition=condition, base=base, combineOn=combineOn, colorBy=colorBy)

#If you might want to run this again in the future (i.e. on different genes),
#you can save the dds object instead of re-running DEfigure
saveRDS(dds, "3m-harmonized-CPNvsOther-dds.rds")

#list of genes you want to plot

genes = c("CACNA2D1","DPYSL3", "PTPRK")

p=data.frame()
for (g in genes) { 
  p1 = plotCounts(dds, g, intgroup = c("condition", "sample", "colorBy"), returnData = T)
  p1$gene = g
  p = rbind(p, p1)
}

##plot!
p$condition = factor(p$condition, levels = c("CPN", "Other"))

for (g in genes){
  pdf(paste0('genes',g, ".pdf"), height = 6, width = 6)
  psub = p[p$gene == g,]
  print(ggplot(psub, aes(x=condition, y=count)) +
          geom_line(aes(group=sample), color="lightgray", size = 0.3) +
          geom_point(aes(color=colorBy))  +
          labs(y = "",x = "") +
          scale_color_manual(name="dataset", values=c('chose_your_colors')) +
          theme(panel.background = element_blank(), 
                axis.title.x = element_blank(), ## title x axis ##blank for it to disappear
                axis.title.y = element_blank(), ## title y axis, disappear
                axis.text.y = element_text(size = 5), 
                axis.text.x = element_text(size = 5), 
                plot.title = element_text(size = 7), 
                axis.ticks.y = element_line(size = 0.2), 
                axis.ticks.x = element_blank(), 
                axis.line.x = element_line(size = 0.2),
                axis.line.y = element_line(size =  0.2)) + NoLegend())
  dev.off()
}


