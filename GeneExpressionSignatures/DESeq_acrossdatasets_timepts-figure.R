#Takes multiple Seurat Objects and performs differential expression analysis 
#Between *the seurat objects*
#By splitting the objects into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise

library(Seurat)
library(DESeq2)
library(ggplot2)
library(plyr)
library(dplyr)

#Load this function as-is for use later
combineDEfigure<-function(seurList,condition="age",batch="dataset",combineOn="organoid", #These are defaults, but will be overrided by what you put below
                        #Below are default settings you can change if you know what you're doing
                        minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                        minSamples=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                        minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                        genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                        form="" #design formula for DESeq2 if you want to include more variables in addition to "~ condition"
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
  summary = do.call(c,lapply(seurList, function(x){table(Idents(x))}))
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
      #return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepGroup)
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[,"sample"] = factor(val[keepOrgs,2])
  colDat[,"batch"] = factor(val[keepOrgs,3])
  colDat[,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Set Up DESeq2")
  design= ~  condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form," condition",sep=""))
  }
  print(design)
  
  dds <- DESeqDataSetFromMatrix(countData = round(data.all),colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]

  return(dds)
}
#rm(seur1)
#Load Seurat Objects
seur23 <- readRDS("23d_ambientRemovedFinal.rds")
seur1 = readRDS("1mo_ambientRemovedFinal.rds")
seur15 = readRDS("1.5mo_ambientRemovedFinal.rds")
seur2 = readRDS("2mo_ambientRemovedFinal.rds")
seur3 = readRDS("3mo_ambientRemovedFinal.rds")
seur4 = readRDS("4mo_ambientRemovedFinal.rds")
seur5 = readRDS("5mo_ambientRemovedFinal.rds")
seur6 = readRDS("6mo_ambientRemovedFinalB.rds")

seur23$CellType_Org <- paste(seur23$FinalName, seur23$Organoid, sep="_")
seur1$CellType_Org <- paste(seur1$FinalName, seur1$Organoid, sep="_")
seur15$CellType_Org <- paste(seur15$FinalName, seur15$Organoid, sep="_")
seur2$CellType_Org <- paste(seur2$FinalName, seur2$Organoid, sep="_")
seur3$CellType_Org <- paste(seur3$FinalName, seur3$Organoid, sep="_")
seur4$CellType_Org <- paste(seur4$FinalName, seur4$Organoid, sep="_")
seur5$CellType_Org <- paste(seur5$FinalName, seur5$Organoid, sep="_")
seur6$CellType_Org <- paste(seur6$FinalName, seur6$Organoid, sep="_")

seur23$celltype <- seur23$FinalName
seur1$celltype <- seur1$FinalName
seur15$celltype <- seur15$FinalName
seur2$celltype <- seur2$FinalName
seur3$celltype <- seur3$FinalName
seur4$celltype <- seur4$FinalName
seur5$celltype <- seur5$FinalName
seur6$celltype <- seur6$FinalName

seurList <- c(seur23, seur1, seur15, seur2, seur3, seur4, seur5, seur6)
rm(seur23, seur1, seur15, seur2, seur3, seur4, seur5, seur6) #to free up memory

#Set condition to the metadata column you want DEGs between
#Make sure all of the seurat objects you are using have metadata columns using the same name
condition = "Age"
#clustersOfInterest = c("23d", "1mo", "1.5mo", "2mo", "3mo", "4mo", "5mo", "6mo")
batch = "Batch" #for coloring figure
#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "CellType_Org"

#run function
dds <- combineDEfigure(seurList, condition=condition, combineOn=combineOn, batch=batch, form="celltype")

#If you might want to run this again in the future (i.e. on different genes),
#you can save the dds object instead of re-running DEfigure
saveRDS(dds, "Expression_across_time_dds2.rds")
#dds = readRDS("Expression_across_time_dds.rds") #Reads in a rds (R dataset) object. Change filename to the object you want

#dds <- readRDS("Expression_across_time_dds2.rds")

mySigGeneFiles <- list.files("DDS", pattern=paste0(celltype,"_DEGs"))
mySigGenes <- data.frame()
for(genFile in mySigGeneFiles){
  tmp <- readRDS(genFile)
  tmp$gene <- rownames(tmp)
  mySigGenes <- rbind(mySigGenes, tmp)
}

mySigGenes <- mySigGenes[order(mySigGenes$padj),]
topGenes <- unique(mySigGenes$gene[1:20])

#list of genes you want to plot
genes = c("PLK1", "MCM2", "MCM7", "CDC20", "NDE1", "TPX2", "EML1", "EML4", "NEUROD6", "NEUROD2", "NELL2", "INHBA", "SYT4", "CUX2", "SLA", "GAD2", "OLIG1", "OLIG2", "GAD1", "CDK4", "TPX2", "MEIS2", "ETV5", "TBR1")
genes = c(genes, "CDK4", "TPX2", "MEIS2", "ETV5", "TBR1")
genes = c(genes,"DLX5", "APOE", "MAGOH", "SOX3", "ID4", "FEZF1", "EBF2", "LHX5",  "FEZF2", "NEUROG2", "NEUROD4", "SATB2", "LMO7", "GRIN2B", "LPL", "BCL11B", "AUTS2", "HEPACAM", "PDGFRA", "LHX1", "LHX9", "FOXP2", "EBF3")
genes = c(genes,"CTGF", "EOMES", "FOXP2", "RELN", "CALB2", "NHLH1", "NEUROD1", "LMO3", "DLG4", "DISC1", "GRIK2", "GRIK5", "GRIN2B", "GRIN3B", "SOX5", "NHLH2", "GPR88", "GRIK1", "PLXNA4", "SEMA3A")
genes = c(genes,"BCL11B", "FEZF2", "SOX5", "PLXNA4", "SEMA3A", "SATB2", "LMO7", "LPL", "BCL11B", "AUTS2", "HEPACAM", "PDGFRA", "LHX1", "LHX9","FOXP2", "NEUROD1", "LMO3")

genes = c(genes,"TPX2", "CDK4", "TBR1", "BCL11B", "CUX2", "APOE", "PDGFRA", "GAD2")
genes <- unique(genes)
#genes <- factor(levels = c("TPX2", "CDK4", "TBR1", "BCL11B", "CUX2", "APOE", "PDGFRA", "GAD2"))
dds <- estimateSizeFactors(dds)
tmp <- counts(dds, normalized=T)[genes,]
myMelt <- melt(tmp)
colnames(myMelt) <- c("gene","OrgType","value")
celltype <- sapply(as.character(myMelt$OrgType), function(x){strsplit(x,"_")[[1]][1]})
Organoid <- sapply(as.character(myMelt$OrgType), function(x){paste(strsplit(x,"_")[[1]][-1], collapse="_")})
Batch <- sapply(as.character(Organoid), function(x){paste(head(strsplit(x,"_")[[1]],-1),collapse="_")})
Age <- sapply(as.character(myMelt$OrgType), function(x){tail(strsplit(x,"_")[[1]],2)[1]})
myMelt$Age <- Age
myMelt$Organoid <- Organoid
myMelt$Batch <- Batch
myMelt$celltype <- celltype
myMelt$AgeType <- paste(myMelt$celltype, myMelt$Age, sep="_")

tmp <- myMelt[myMelt$gene==genes[1],] %>% group_by(AgeType) %>% summarise(mean=mean(value), sd = sd(value))
tmp$gene <- genes[1]
for(myGene in genes[-1]){
  tmp2 <- myMelt[myMelt$gene==myGene,] %>% group_by(AgeType) %>% summarise(mean=mean(value), sd = sd(value))
  tmp2$gene <- myGene
  tmp <- rbind(tmp, tmp2)
}
tmp$Age <- sapply(tmp$AgeType, function(x){strsplit(x, "_")[[1]][2]})
tmp$celltype <- sapply(tmp$AgeType, function(x){strsplit(x, "_")[[1]][1]})

myMeltFull <- myMelt

plotGeneByCellType <- function(geneList=c("APOE","PDGFRA","GAD2"), cellTypes=c("aRG","IP","oRG"), title=""){
  tmpSub <- counts(dds, normalized=T)[geneList,dds$celltype %in% cellTypes]
  mySubMelt <- melt(tmpSub)
  colnames(mySubMelt) <- c("gene","OrgType","value")
  celltype <- sapply(as.character(mySubMelt$OrgType), function(x){strsplit(x,"_")[[1]][1]})
  Organoid <- sapply(as.character(mySubMelt$OrgType), function(x){paste(strsplit(x,"_")[[1]][-1], collapse="_")})
  Batch <- sapply(as.character(Organoid), function(x){paste(head(strsplit(x,"_")[[1]],-1),collapse="_")})
  Age <- sapply(as.character(mySubMelt$OrgType), function(x){tail(strsplit(x,"_")[[1]],2)[1]})
  mySubMelt$Age <- Age
  mySubMelt$Organoid <- Organoid
  mySubMelt$Batch <- Batch
  mySubMelt$celltype <- celltype
  mySubMelt$AgeType <- paste(mySubMelt$celltype, mySubMelt$Age, sep="_")
  tmpSub <- mySubMelt[mySubMelt$gene==geneList[1],] %>% group_by(AgeType) %>% summarise(mean=mean(value), sd = sd(value))
  tmpSub$gene <- geneList[1]
  for(myGene in geneList[-1]){
    tmp2 <- mySubMelt[mySubMelt$gene==myGene,] %>% group_by(AgeType) %>% summarise(mean=mean(value), sd = sd(value))
    tmp2$gene <- myGene
    tmpSub <- rbind(tmpSub, tmp2)
  }
  tmpSub$Age <- sapply(tmpSub$AgeType, function(x){strsplit(x, "_")[[1]][2]})
  tmpSub$celltype <- sapply(tmpSub$AgeType, function(x){strsplit(x, "_")[[1]][1]})
  tmpSub <- tmpSub[tmpSub$gene %in% geneList & tmpSub$celltype %in% cellTypes,]
  myExpan <- expand.grid(gene=levels(factor(tmpSub$gene)), Age=levels(factor(tmpSub$Age)), celltype=levels(factor(tmpSub$celltype)))
  myExpan <- cbind(AgeType=paste(myExpan$celltype, myExpan$Age, sep="_"), mean=NA, sd=NA, myExpan)
  isIn <- paste(myExpan$AgeType,myExpan$gene) %in% paste(tmpSub$AgeType, tmpSub$gene)
  tmpSubFull <- rbind(tmpSub, myExpan[!isIn,])
  p1 <- ggplot(tmpSubFull, aes(x=factor(Age,levels=c("d23","1mo","1.5mo","2mo","3mo","4mo","5mo","6mo")), y=mean, group=1)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), fill="gray90")+
    geom_line(color="black") +
    geom_point(data=mySubMelt, aes(x=Age, y=value,
                                   color=factor(Batch,levels=c("PGP1_d23","Mito210c1_d23",
                                                               "Mito210_PTEN_1mo","Mito210_d28_1mo","Mito_210_d35_1mo","Other_1mo",
                                                               "PGP1_1.5mo","Mito210c1_1.5mo",
                                                               "PGP1_2mo","Mito210c1_2mo",
                                                               "1_3mo","2_3mo","3_3mo","4_3mo","5_3mo","6_3mo","7_3mo",
                                                               "PGP1_4mo","Mito210c1_4mo",
                                                               "PGP1_5mo","Mito210c1_5mo",
                                                               "11a_6mo","GM_6mo","Mito210_6mo","Mito210_PTEN_6mo","HUES66_6mo","PGP1_b1_6mo","PGP1_b3_6mo")))) +
    scale_color_manual(values=c("#A63603","#A1D99B",
                                "#377EB8","#4DAF4A", "#A8DDB5","#984EA3",
                                "#A63603","#A1D99B",
                                "#A63603","#A1D99B",
                                "#984EA3","#4DAF4A", '#377EB8',"#AE017E", "#F768A1",'#FD8D3C', "#CC4C02",
                                "#A63603","#A1D99B",
                                "#A63603","#A1D99B",
                                "#E41A1C", "#984EA3", "#AE017E", "#4DAF4A", '#377EB8','#FD8D3C', "#FEC44F"))+
  
    facet_grid(rows=factor(gene,levels=geneList) ~ factor(celltype, levels=cellTypes), scales="free") +
    theme_minimal() + NoLegend() +
    xlab("Age") + ylab("Normalized Counts") +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.line=element_line(colour = "black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          strip.background=element_rect(fill="darkgray"))
  return(p1)
}

##
#c(1,2,3,4,5,6,7)
#c("GM","HUES66","HUES66_d106","Mito210","Mito210_PTEN","PGP1_b1","PGP1_b2")
#c("#984EA3","#4DAF4A", '#377EB8',"#AE017E", "#F768A1",'#FD8D3C', "#CC4C02")
#For 3 months: "#984EA3", "#AE017E", "#F768A1", "#4DAF4A", '#377EB8','#FD8D3C', "#CC4C02"
#               GM08330, Mito210 c1 bx, #Mito210 c2 bx, HUES66, HUES66, PGP1 c2 bx, PGP1 c2 bx.

plotGeneByCellType(geneList=c("TPX2","CDK4","TBR1","BCL11B","CUX2","APOE","PDGFRA","GAD2"), cellTypes=c("aRG","IP","oRG"))
plotGeneByCellType(geneList=c("APOE","PDGFRA","GAD2","EOMES"), cellTypes=c("aRG","IP","CPN"))

#  scale_color_manual(name="dataset", values=c("#e41a1c",'#984ea3',"#984ea3","#ae017e","#ae017e","#f768a1","#4daf4a", "#4daf4a" , "#4daf4a" ,'#377eb8' ,'#377eb8' ,'#fd8d3c', '#fd8d3c','#cc4c02','#cc4c02')) +
#  theme_minimal() + NoLegend() + ggtitle()

ddsARG <- dds[,dds$celltype=="aRG"]
ddsARG$Age <- as.numeric(as.character(mapvalues(ddsARG$condition, from=c("d23","1mo","1.5mo","2mo","3mo","4mo","5mo","6mo"), to=c(23,30,45,60,90,120,150,180))))
design(ddsARG) <- ~Age
ddsARG <- DESeq(ddsARG)
aRGout <- data.frame(results(ddsARG))
aRGout <- aRGout[order(aRGout$pvalue),]


RunDEGOverTimeCellType <- function(celltype){
  ddsTMP <- dds[,dds$celltype==celltype]
  ddsTMP$Age <- as.numeric(as.character(mapvalues(ddsTMP$condition, from=c("d23","1mo","1.5mo","2mo","3mo","4mo","5mo","6mo"), to=c(23,30,45,60,90,120,150,180))))
  design(ddsTMP) <- ~Age
  if(length(table(ddsTMP$Age))<2){
    return(NULL)
  }
  ddsTMP <- DESeq(ddsTMP)
  TMPout <- data.frame(results(ddsTMP))
  TMPout <- TMPout[order(TMPout$pvalue),]
  return(TMPout)
}
myDEGs <-list()
for(celltype in unique(myMelt$celltype)){
  if(!celltype %in% c("ThisisAPlaceholder")){
    print(paste("Processing", celltype))
    myDEGs[[celltype]] <- RunDEGOverTimeCellType(celltype)
  }
}
saveRDS(myDEGs, "DEGs_Celltypes_Over_Time.RDS")

for(i in names(myDEGs)){
  print(paste(i, "Sig:", sum(myDEGs[[i]]$padj[!is.na(myDEGs[[i]]$padj)]<0.05), "Insignificant:", sum(!myDEGs[[i]]$padj[!is.na(myDEGs[[i]]$padj)]<0.05)))
}

GOup <- list()
GOdown <- list()
for(celltype in names(myDEGs)){
  Degs <- myDEGs[[celltype]]
  Degs <- Degs[!is.na(Degs$padj),]
  Degup <- Degs[Degs$log2FoldChange > 0 & Degs$padj < 0.05,]
  print(paste(celltype,":", nrow(Degup),"sig up."))
  Degdown <- Degs[Degs$log2FoldChange < 0 & Degs$padj < 0.05,]
  print(paste(celltype,":", nrow(Degdown),"sig down."))
  GOup[[celltype]] <- rownames(Degup)[1:min(1000, nrow(Degup))]
  GOdown[[celltype]] <- rownames(Degdown)[1:min(1000, nrow(Degdown))]
}

library(clusterProfiler)
library(org.Hs.eg.db)
GOupRes <- list()
GOdownRes <- list()
for(celltype in names(GOup)){
  GOupRes[[celltype]] <- enrichGO(GOup[[celltype]], OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                                  minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
  GOdownRes[[celltype]] <- enrichGO(GOdown[[celltype]], OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                                  minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
}

simGOupRes <- list()
simGOdownRes <- list()
for(celltype in names(GOupRes)){
  simGOupRes[[celltype]] <- simplify(GOupRes[[celltype]], cutoff=0.4)
  simGOdownRes[[celltype]] <- simplify(GOdownRes[[celltype]], cutoff=0.4)
}

saveRDS(GOupRes, "GO_up_over_time.RDS")
saveRDS(GOdownRes, "GO_down_over_time.RDS")
saveRDS(simGOupRes, "Simple_GOup.RDS")
saveRDS(simGOdownRes, "Simple_GOdown.RDS")

aRG_topGO <- rbind(simGOupRes$aRG@result[1:10,], simGOdownRes$aRG@result[10:1,])
aRG_topGO$Direction <- c(rep("up",10),rep("down",10))
ggplot(aRG_topGO, aes(x=factor(Direction, levels=c("up","down")), y=factor(Description, levels=rev(Description)), size=Count, color=p.adjust)) +
  geom_point() +
  scale_color_gradient(low="orange", high="purple") +
  xlab("Direcion of change over time") + ylab("") +
  ggtitle("aRGs") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=16),
        axis.text=element_text(size=12),
        axis.title.x=element_text(face="bold",size=12))

PlotGOs <- function(celltype, title=""){
  if(title==""){
    title <- celltype
  }
  topGOs <- rbind(simGOupRes[[celltype]]@result[1:min(10,nrow(simGOupRes[[celltype]]@result)),], simGOdownRes[[celltype]]@result[min(10,nrow(simGOdownRes[[celltype]]@result)):1,])
  topGOs$Direction <- c(rep("up",min(10,nrow(simGOupRes[[celltype]]@result))),rep("down",min(10,nrow(simGOdownRes[[celltype]]@result))))
  p2 <- ggplot(topGOs, aes(x=factor(Direction, levels=c("up","down")), y=factor(Description, levels=rev(Description)), size=Count, color=p.adjust)) +
    geom_point() +
    scale_color_gradient(low="orange", high="purple") +
    xlab("Direcion of change over time") + ylab("") +
    ggtitle(title) +
    theme(plot.title=element_text(face="bold", hjust=0.5, size=16),
          axis.text=element_text(size=12),
          axis.title.x=element_text(face="bold",size=12))
  return(p2)
}

plotGeneByCellType(geneList = strsplit(simGOupRes$IP@result$geneID[1],"/")[[1]][1:14], "IP")
plotGeneByCellType(geneList = strsplit(simGOupRes$IP@result$geneID[1],"/")[[1]][15:28], "IP")
plotGeneByCellType(geneList = strsplit(simGOupRes$IP@result$geneID[1],"/")[[1]][29:42], "IP")


RunDEGEachTimeCellType <- function(celltype, age){
  ddsTMP <- dds[,dds$celltype==celltype]
  ddsTMP$Focus <- (ddsTMP$Age==age)
  design(ddsTMP) <- ~ Focus
  if(length(table(ddsTMP$Focus))<2){
    return(NULL)
  }
  ddsTMP <- DESeq(ddsTMP)
  TMPout <- data.frame(results(ddsTMP))
  TMPout <- TMPout[order(TMPout$pvalue),]
  return(TMPout)
}
oRG_DEGs <-list()
for(age in names(table(dds$Age[dds$celltype=="oRG"]))[table(dds$Age[dds$celltype=="oRG"])>0]){
  if(!celltype %in% c("ThisisAPlaceholder")){
    print(paste("Processing", age))
    oRG_DEGs[[age]] <- RunDEGEachTimeCellType(celltype = "oRG", age = age)
  }
}
aRG_DEGs <-list()
for(age in names(table(dds$Age[dds$celltype=="aRG"]))[table(dds$Age[dds$celltype=="aRG"])>0]){
  if(!celltype %in% c("ThisisAPlaceholder")){
    print(paste("Processing", age))
    aRG_DEGs[[age]] <- RunDEGEachTimeCellType(celltype = "aRG", age = age)
  }
}
IP_DEGs <-list()
for(age in names(table(dds$Age[dds$celltype=="IP"]))[table(dds$Age[dds$celltype=="IP"])>0]){
  if(!celltype %in% c("ThisisAPlaceholder")){
    print(paste("Processing", age))
    IP_DEGs[[age]] <- RunDEGEachTimeCellType(celltype = "IP", age = age)
  }
}
CPN_DEGs <-list()
for(age in names(table(dds$Age[dds$celltype=="CPN"]))[table(dds$Age[dds$celltype=="CPN"])>0]){
  if(!celltype %in% c("ThisisAPlaceholder")){
    print(paste("Processing", age))
    CPN_DEGs[[age]] <- RunDEGEachTimeCellType(celltype = "CPN", age = age)
  }
}

for(i in names(oRG_DEGs)){
  print(paste(i, "up:", paste(rownames(oRG_DEGs[[i]][oRG_DEGs[[i]]$log2FoldChange>0,])[1:50], collapse=", ")))
  print(paste(i, "down:", paste(rownames(oRG_DEGs[[i]][oRG_DEGs[[i]]$log2FoldChange<0,])[1:50], collapse=", ")))
}

aRGtop <- list(day23=(aRG_DEGs$d23[aRG_DEGs$d23$log2FoldChange>0 & aRG_DEGs$d23$padj<0.05 & !is.na(aRG_DEGs$d23$padj),]),
                     mo1=(aRG_DEGs$`1mo`[aRG_DEGs$`1mo`$log2FoldChange>0 & aRG_DEGs$`1mo`$padj<0.05 & !is.na(aRG_DEGs$`1mo`$padj),]),
                     mo1.5=(aRG_DEGs$`1.5mo`[aRG_DEGs$`1.5mo`$log2FoldChange>0 & aRG_DEGs$`1.5mo`$padj<0.05 & !is.na(aRG_DEGs$`1.5mo`$padj),]),
                     mo2=(aRG_DEGs$`2mo`[aRG_DEGs$`2mo`$log2FoldChange>0 & aRG_DEGs$`2mo`$padj<0.05 & !is.na(aRG_DEGs$`2mo`$padj),]),
                     mo3=(aRG_DEGs$`3mo`[aRG_DEGs$`3mo`$log2FoldChange>0 & aRG_DEGs$`3mo`$padj<0.05 & !is.na(aRG_DEGs$`3mo`$padj),]),
                     mo4=(aRG_DEGs$`4mo`[aRG_DEGs$`4mo`$log2FoldChange>0 & aRG_DEGs$`4mo`$padj<0.05 & !is.na(aRG_DEGs$`4mo`$padj),]),
                     mo5=(aRG_DEGs$`5mo`[aRG_DEGs$`5mo`$log2FoldChange>0 & aRG_DEGs$`5mo`$padj<0.05 & !is.na(aRG_DEGs$`5mo`$padj),]),
                     mo6=(aRG_DEGs$`6mo`[aRG_DEGs$`6mo`$log2FoldChange>0 & aRG_DEGs$`6mo`$padj<0.05 & !is.na(aRG_DEGs$`6mo`$padj),]))
for(i in names(aRGtop)){
  aRGtop[[i]]$gene  <- rownames(aRGtop[[i]])
}
write.xlsx(aRGtop, "aRG_DEGs_AllTimepoints.xlsx")
oRGtop <- list(mo2=(oRG_DEGs$`2mo`[oRG_DEGs$`2mo`$log2FoldChange>0 & oRG_DEGs$`2mo`$padj<0.05 & !is.na(oRG_DEGs$`2mo`$padj),]),
               mo3=(oRG_DEGs$`3mo`[oRG_DEGs$`3mo`$log2FoldChange>0 & oRG_DEGs$`3mo`$padj<0.05 & !is.na(oRG_DEGs$`3mo`$padj),]),
               mo4=(oRG_DEGs$`4mo`[oRG_DEGs$`4mo`$log2FoldChange>0 & oRG_DEGs$`4mo`$padj<0.05 & !is.na(oRG_DEGs$`4mo`$padj),]),
               mo5=(oRG_DEGs$`5mo`[oRG_DEGs$`5mo`$log2FoldChange>0 & oRG_DEGs$`5mo`$padj<0.05 & !is.na(oRG_DEGs$`5mo`$padj),]),
               mo6=(oRG_DEGs$`6mo`[oRG_DEGs$`6mo`$log2FoldChange>0 & oRG_DEGs$`6mo`$padj<0.05 & !is.na(oRG_DEGs$`6mo`$padj),]))
for(i in names(oRGtop)){
  oRGtop[[i]]$gene  <- rownames(oRGtop[[i]])
}
write.xlsx(oRGtop, "oRG_DEGs_AllTimepoints.xlsx")

IPtop <- list(day23=(IP_DEGs$d23[IP_DEGs$d23$log2FoldChange>0 & IP_DEGs$d23$padj<0.05 & !is.na(IP_DEGs$d23$padj),]),
               mo1=(IP_DEGs$`1mo`[IP_DEGs$`1mo`$log2FoldChange>0 & IP_DEGs$`1mo`$padj<0.05 & !is.na(IP_DEGs$`1mo`$padj),]),
               mo1.5=(IP_DEGs$`1.5mo`[IP_DEGs$`1.5mo`$log2FoldChange>0 & IP_DEGs$`1.5mo`$padj<0.05 & !is.na(IP_DEGs$`1.5mo`$padj),]),
               mo2=(IP_DEGs$`2mo`[IP_DEGs$`2mo`$log2FoldChange>0 & IP_DEGs$`2mo`$padj<0.05 & !is.na(IP_DEGs$`2mo`$padj),]),
               mo3=(IP_DEGs$`3mo`[IP_DEGs$`3mo`$log2FoldChange>0 & IP_DEGs$`3mo`$padj<0.05 & !is.na(IP_DEGs$`3mo`$padj),]),
               mo4=(IP_DEGs$`4mo`[IP_DEGs$`4mo`$log2FoldChange>0 & IP_DEGs$`4mo`$padj<0.05 & !is.na(IP_DEGs$`4mo`$padj),]),
               mo5=(IP_DEGs$`5mo`[IP_DEGs$`5mo`$log2FoldChange>0 & IP_DEGs$`5mo`$padj<0.05 & !is.na(IP_DEGs$`5mo`$padj),]),
               mo6=(IP_DEGs$`6mo`[IP_DEGs$`6mo`$log2FoldChange>0 & IP_DEGs$`6mo`$padj<0.05 & !is.na(IP_DEGs$`6mo`$padj),]))
for(i in names(IPtop)){
  IPtop[[i]]$gene  <- rownames(IPtop[[i]])
}
write.xlsx(IPtop, "IP_DEGs_AllTimepoints.xlsx")

aRGbestGO <- rbind(aRGsim$d23@result[1:5,],aRGsim$`1mo`@result[1:5,],aRGsim$`1.5mo`@result[1:5,],aRGsim$`3mo`@result[1:5,],aRGsim$`4mo`@result[1:2,],aRGsim$`5mo`@result[1:5,],aRGsim$`6mo`@result[1:5,])
aRGbestGO <- rbind(aRGbestGO, c(NA, "Golgi vesicle transport", NA, NA, NA, 0.05, 0.05, NA, NA, "2mo"))
aRGbestGO$Count <- as.numeric(aRGbestGO$Count)
aRGbestGO$p.adjust <- as.numeric(aRGbestGO$p.adjust)
ggplot(aRGbestGO, aes(x=factor(Age, levels=c("d23","1mo","1.5mo","2mo","3mo","4mo","5mo","6mo")), y=factor(Description, levels=rev(unique(Description))), size=Count, color=p.adjust)) +
  geom_point() +
  scale_color_gradient(low="orange", high="purple")+
  xlab("Age") + ylab("GO term") + ggtitle("aRGs")

aRG_3genes <- c('ARF1', 'KDELR2', 'YIPF5', 'TMED3', 'SEC23IP', 'TMED9', 'GOSR2', 'GOLGA5', 'ARF4', 'COPZ1', 'BET1', 'SEC22A', 'SNX12')

##datasets and batches
#For 1 month:  "#984EA3", "#4DAF4A", "#A8DDB5", "#377EB8" , GM08330  Mito210 c1 bx, Mito210 c1 bx, Mito210 c2
#For 1.5 and 2 months, "#A1D99B", "#A63603", Mito210 c1, PGP1 c1
#4 and 5 months: "#A1D99B", "#A63603", "#F4A582", Mito210 c1, PGP1 c1
#For 3 months: "#984EA3", "#AE017E", "#F768A1", "#4DAF4A", '#377EB8','#FD8D3C', "#CC4C02", GM08330, Mito210 c1 bx, #Mito210 c2 bx, HUES66, HUES66, PGP1 c2 bx, PGP1 c2 bx.
#For 6 months: "#E41A1C", "#984EA3", "#AE017E", "#4DAF4A", '#377EB8','#FD8D3C', "#FEC44F", 11a, GM08330, Mito210 c1, #Mito210 c2, HUES66 , PGP1 c2 cx, PGP1 c2 bx 
##genes to plot 
# c("PLK1", "MCM2", "NEUROD2", "GRIN2B", "PDE1A", SOX5", "INHBA", "SATB2", "APOE", "OLIG1")






