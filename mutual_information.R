library(aricode)
library(ggplot2)
#AMI: cite Vinh et al "Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance"

setwd("~/Documents/AtlasAna/mutualInformation/")

seurfs = c("../FinalObjects/23days/23d_harmonizedObj_060421.rds",
           "../FinalObjects/1month/harmonizedObj_1mo_111520_final.rds",
           "../FinalObjects/1.5month/1.5mo_harmonizedObj_060421.rds",
           "../FinalObjects/2months/2mo_harmonizedObj_060421.rds",
           "../FinalObjects/3months/3mo_harm_111120.rds",
           "../FinalObjects/4months/4mo_harm_060421.rds",
           "../FinalObjects/5months/5mo_harmonizedObj_060421_3.rds",
           "../FinalObjects/6months/6mo_harm_111120.rds")
ages = c("23d", "1m", "1.5m", "2m","3m","4m", "5m","6m")

mis = data.frame()
for (j in 1:length(seurfs)) {
  seur = readRDS(seurfs[[j]])
  meta = seur@meta.data
  rm(seur)
  
  if ("dataset" %in% colnames(meta)) {
    for (i in levels(factor(meta$dataset))) {
      meta1 = meta[meta$dataset==i,]
      if(!"org" %in% colnames(meta1)) {meta1$org = meta1$orig.ident}
      
      #Mutual information between clusters and organoids
      adjusted_mi = AMI(meta1$FinalName, meta1$org)
  
      miR = data.frame("age" = ages[[j]], "dataset" = i, "mi" = adjusted_mi)
      mis = rbind(mis, miR)
      
      counts = as.matrix(table(meta1$FinalName,meta1$org))
      pdf(paste0(ages[[j]],"_",i, ".pdf"))
      print(ggplot(melt(counts),aes(y=value, x=Var2,fill=Var1)) + 
        geom_bar(stat="identity", position = position_fill()))
      dev.off()
    }
  } else {
    if(!"org" %in% colnames(meta)) {meta$org = meta$orig.ident}
    adjusted_mi = AMI(meta$FinalName, meta$org)
    miR = data.frame("age" = ages[[j]], "dataset" = 1, "mi" = adjusted_mi)
    mis = rbind(mis, miR)
    
    counts = as.matrix(table(meta$FinalName,meta$org))
    pdf(paste0(ages[[j]],"_",i, ".pdf"))
    print(ggplot(melt(counts),aes(y=value, x=Var2,fill=Var1)) + 
      geom_bar(stat="identity", position = position_fill()))
    dev.off()
    
  }
}

mis$tech = "RNA"
mis$cellline = as.character(mis$dataset)
mis$cellline[mis$age=="1m" & mis$dataset==1 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="1m" & (mis$dataset==2 | mis$dataset==3) & mis$tech=="RNA"] = "Mito210 c1"
mis$cellline[mis$age=="1m" & mis$dataset==4 & mis$tech=="RNA"] = "Mito210 c2"
mis$cellline[mis$age=="3m" & mis$dataset==1 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="3m" & mis$dataset==2 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="3m" & mis$dataset==3 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="3m" & mis$dataset==4 & mis$tech=="RNA"] = "Mito210 c1"
mis$cellline[mis$age=="3m" & mis$dataset==5 & mis$tech=="RNA"] = "Mito210 c2"
mis$cellline[mis$age=="3m" & mis$dataset %in% 6:7 & mis$tech=="RNA"] = "PGP1 c2"
mis$cellline[mis$age=="6m" & mis$dataset==1 & mis$tech=="RNA"] = "11a"
mis$cellline[mis$age=="6m" & mis$dataset==2 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="6m" & mis$dataset==3 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="6m" & mis$dataset==4 & mis$tech=="RNA"] = "Mito210 c1"
mis$cellline[mis$age=="6m" & mis$dataset==5 & mis$tech=="RNA"] = "Mito210 c2"
mis$cellline[mis$age=="6m" & mis$dataset %in% 6:7 & mis$tech=="RNA"] = "PGP1 c2"
mis$cellline[mis$cellline=="Mito210c1"] = "Mito210 c1"

saveRDS(mis, "mutualinformations.rds")

#human
seur = readRDS("~/Documents/HumanData/CZI/CZI_060421.rds")
meta = seur@meta.data
czi = AMI(meta$FinalName, meta$orig.ident) ##This is one donor per age?!
czi = 0.02579153

seur = readRDS("~/Documents/HumanData/TrevinoRNA/trevinoObj.rds")
meta = seur@meta.data
trevino = AMI(meta$CellTypes, meta$Tissue.ID)
trevino = 0.06254533

seur = readRDS("~/Documents/HumanData/Geschwind2019/Geschwin_clusteredSeur.rds")
meta = seur@meta.data[!is.na(seur$Cluster),]
poliou = AMI(meta$Cluster, meta$Donor)
poliou = 0.01347987

counts = as.matrix(table(meta$Cluster,meta$Donor))
pdf("poliou.pdf")
print(ggplot(melt(counts),aes(y=value, x=Var2,fill=Var1)) + 
        geom_bar(stat="identity", position = position_fill()))
dev.off()

#plot!
ggplot(mis[mis$tech=="RNA",], aes(age, mi, color = cellline)) +
  geom_jitter(size=3, alpha=0.6, height=0, width=0.2) +
  theme_classic() +
  geom_hline(yintercept=trevino, linetype="dashed", color="#3B3B3B") +
  geom_hline(yintercept=czi, linetype="dashed", color="#696969") +
  geom_hline(yintercept=poliou, linetype="dashed", color="#AEAEAE") +
  scale_color_manual(values=c("#e41a1c", '#984ea3', '#ae017e', "#4daf4a", '#377eb8',  '#fd8d3c',  '#cc4c02')) +
  ylab("AMI")
ggsave("mutal_information.pdf", width=4, height=3)

#Share-seq
mis$platform = "10X"
share = readRDS("../ShareSeq/FinalObjects/seur_share-seq_RNA.rds")
meta = share@meta.data
meta$age = as.character(meta$sample)
meta$age[meta$sample == "day_23"] = "23d"
meta$age[meta$sample == "day_35"] = "1m"
meta$age[meta$sample == "day_59"] = "2m"
meta$age[meta$sample == "day_90"] = "3m"
for (a in unique(meta$age)) {
  meta1 = meta[meta$age==a,]
  adjusted_mi = AMI(meta1$FinalName, meta1$replicate)
  miR = data.frame(age = a, dataset = "PGP1", mi = adjusted_mi, tech = "RNA", cellline="PGP1", platform="ShareSeq")
  mis = rbind(mis, miR)
}
ggplot(mis[mis$tech =="RNA",], aes(age, mi, color = cellline, shape=platform)) +
  geom_jitter(size=3, alpha=0.6, height=0, width=0.2) +
  theme_classic() +
  geom_hline(yintercept=trevino, linetype="dashed", color="#3B3B3B") +
  geom_hline(yintercept=czi, linetype="dashed", color="#696969") +
  geom_hline(yintercept=poliou, linetype="dashed", color="#AEAEAE") +
  scale_color_manual(values=c("#e41a1c", '#984ea3', '#ae017e', "#4daf4a", '#377eb8',  '#fd8d3c',  '#cc4c02')) +
  ylab("AMI")
ggsave("mutal_information-shareseq.pdf", width=4, height=3)


#ATAC
atacs = c("~/Documents/scATAC/SUV_Mito210_1m/wt/Signac/clusteredSign.rds",
          "~/Documents/scATAC/SUV_Mito210_3m/wt/Signac/clusteredSign.rds",
          "~/Documents/scATAC/SUV_Mito210_6m/wt/Signac/clusteredSign.rds")
ages = c("1m","3m", "6m")

for (j in 1:length(atacs)) {
  seur = readRDS(atacs[[j]])
  meta = seur@meta.data
  rm(seur)
  
  adjusted_mi = AMI(meta$CellType, meta$orig.ident)
  
  miR = data.frame("age" = ages[[j]], "dataset" = i, "mi" = adjusted_mi, "cellline"="Mito210", tech="ATAC", platform="10X")
  mis = rbind(mis, miR)
  
  counts = as.matrix(table(meta$CellType,meta$orig.ident))
  pdf(paste0(ages[[j]],"_atac", ".pdf"))
  print(ggplot(melt(counts),aes(y=value, x=Var2,fill=Var1)) + 
          geom_bar(stat="identity", position = position_fill()))
  dev.off()
}

trevino_atac = 0.04999161

ggplot(mis, aes(age, mi, color = tech, shape=platform)) +
  geom_jitter(size=3, alpha=0.6, height=0, width=0.2) +
  theme_classic() +
  geom_hline(yintercept=trevino_atac, linetype="dashed", color="#3B3B3B") +
  scale_color_manual(values=c("#4DAF4A", "gray")) +
  ylab("AMI")
ggsave("mutual_information-atacAndRNA.pdf", width=4, height=3)


saveRDS(mis, "mutualinformations.rds")
