#1/12/21

library(Seurat)
library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

seur = readRDS("~/Documents/Arlotta/Organoids/Ana_Amanda/July/All/1mo_fordeseq_121220/progs_010821.rds")
seur = subset(seur, subset=FinalName %in% c("aRG", "IP"))

#Loop through each cell type in this subset of this seur object
celltypes = c("aRG", "IP")
for (ct in celltypes) {
  
  ctNoSpace = gsub(" ",".", ct)

  sig = read_xlsx(paste0("~/Documents/Arlotta/Organoids/Ana_Amanda/July/All/1mo_fordeseq_121220/progs/", ct, ".vs.Other.xlsx"))

  genes = sig[sig$log2FoldChange>0,]$gene
  genes = genes[1:200]
  genes = sapply(strsplit(genes, ".", fixed=T), "[[", 1)

  seur = AddModuleScore(seur, list(genes), name=ctNoSpace)
  seur$CellTypeOfInterest = ifelse(seur$FinalName==ct, ct, "Other")

  df = seur@meta.data[,c("dataset",paste0(ctNoSpace, "1"),"CellTypeOfInterest", "org")]

  df$org = factor(df$org)
  df$CellTypeOfInterest = factor(df$CellTypeOfInterest, levels = c(ct, "Other"))
  df$dataset = factor(df$dataset)
  dfm = melt(df)
  ReplicateAverages <- dfm %>% group_by(org,CellTypeOfInterest,dataset) %>% summarize(mean(value), .groups = "keep")
  colnames(ReplicateAverages)[[4]] = 'value'
  pdf(paste0("~/Documents/Arlotta/Organoids/Ana_Amanda/July/All/1mo_fordeseq_121220/progs/", ctNoSpace, ".SigModuleScore.pdf"), width=5, height=3)
  print(ggplot(dfm, aes(x=CellTypeOfInterest,y=value)) + 
    geom_line(data=ReplicateAverages, aes(group=org), color="lightgray", position = position_dodge2(0.3)) + 
    geom_violin(fill=NA) + 
    geom_point(data=ReplicateAverages, size=3, cex = 4, aes(group=org,color=dataset), position=position_dodge2(0.3))  + 
    scale_color_manual(name="dataset", values = c("red","green", "blue", "purple", "orange", "yellow","red4")) +
    theme_classic() +
    labs(y="Module Expression Score", x=NULL))
  dev.off()


  ReplicateAverages1 = ReplicateAverages[ReplicateAverages$CellTypeOfInterest!="Other",]
  ReplicateAverages2 = ReplicateAverages[ReplicateAverages$CellTypeOfInterest=="Other",]
  res = aov(value ~ org, data=ReplicateAverages1)
  sink("SigModuleOutput.txt", append = T)
  print(ct)
  print(summary(res))
  print("Other")
  res = aov(value ~ org, data=ReplicateAverages2)
  print(summary(res))
  print("")
  sink()
}
