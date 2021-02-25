library(mpmi)
library(ggplot2)

setwd("~/Documents/wt_merge/")

seurfs = c("~/Documents/wt_merge/FinalObjects/23d_111520.rds",
           "~/Documents/wt_merge/FinalObjects/1mo_111520_final.rds",
           "~/Documents/wt_merge/FinalObjects/2mo_111520_final.rds",
           "~/Documents/wt_merge/FinalObjects/3mo_harmonizedObj_102820.rds",
           "~/Documents/wt_merge/FinalObjects/Mito210c1_4mo_022321.rds",
           "~/Documents/wt_merge/FinalObjects/Mito210c1_5mo_022321.rds",
           "~/Documents/wt_merge/FinalObjects/6mo_harmonizedObj_103020.rds")
ages = c("23d", "1m", "2m","3m","4m", "5m","6m")

mis = data.frame()
for (j in 1:length(seurfs)) {
  seur = readRDS(seurfs[[j]])
  meta = seur@meta.data
  rm(seur)
  
  if ("dataset" %in% colnames(meta)) {
    for (i in levels(factor(meta$dataset))) {
      meta1 = meta[meta$dataset==i,]

      #Mutual information between clusters and organoids
      assigns = meta1[,c("FinalName","org")]
      real_mi = dmi(assigns)$mi[1,2]
      #randomized_mi = list()
      #for (i in 1:1000){
      #  randomClusts = sample(assigns$FinalName)
      #  assigns$FinalName = randomClusts
      #  randomized_mi[i] = as.numeric(dmi(assigns)$mi[1,2])
      #}
      #randomized_mi = as.numeric(randomized_mi)
      #z = (real_mi - mean(randomized_mi)) / sd(randomized_mi)
  
      miR = data.frame("age" = ages[[j]], "dataset" = i, "mi" = real_mi)
      mis = rbind(mis, miR)
    }
  } else {
    assigns = meta[,c("FinalName","org")]
    real_mi = dmi(assigns)$mi[1,2]
    miR = data.frame("age" = ages[[j]], "dataset" = 1, "mi" = real_mi)
    mis = rbind(mis, miR)
  }
}

mis$cellline = "Mito210 c1"
mis$tech = "RNA"
mis$cellline[mis$age=="1m" & mis$dataset==1 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="1m" & mis$dataset==4 & mis$tech=="RNA"] = "Mito210 c2"
mis$cellline[mis$age=="3m" & mis$dataset==1 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="3m" & mis$dataset==2 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="3m" & mis$dataset==3 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="3m" & mis$dataset==5 & mis$tech=="RNA"] = "Mito210 c2"
mis$cellline[mis$age=="3m" & mis$dataset %in% 6:7 & mis$tech=="RNA"] = "PGP1"
mis$cellline[mis$age=="6m" & mis$dataset==1 & mis$tech=="RNA"] = "11a"
mis$cellline[mis$age=="6m" & mis$dataset==2 & mis$tech=="RNA"] = "GM"
mis$cellline[mis$age=="6m" & mis$dataset==3 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="6m" & mis$dataset==5 & mis$tech=="RNA"] = "HUES66"
mis$cellline[mis$age=="6m" & mis$dataset %in% 6:7 & mis$tech=="RNA"] = "PGP1"

ggplot(mis[mis$tech=="RNA",], aes(age, mi, color = cellline)) +
  geom_jitter(size=3, alpha=0.6, height=0, width=0.2) +
  theme_classic() +
  geom_hline(yintercept=0.008, linetype="dashed", color="gray72") +
  geom_hline(yintercept=0.064, linetype="dashed", color="gray72") +
  scale_color_manual(values=c("#E41A1C", '#984EA3', "#AE017E", "#4DAF4A", '#377EB8', '#FD8D3C'))
ggsave("mutal_information.pdf", width=4, height=3)

saveRDS(mis, "mutualinformations.rds")


#ATAC
atacs = c("~/Documents/scATAC/SUV_Mito210_1m/wt/Signac/clusteredSign.rds",
          "~/Documents/scATAC/SUV_Mito210_3m/wt/Signac/clusteredSign.rds",
          "~/Documents/scATAC/SUV_Mito210_6m/wt/Signac/clusteredSign.rds")
ages = c("1m","3m", "6m")

for (j in 1:length(atacs)) {
  seur = readRDS(atacs[[j]])
  meta = seur@meta.data
  rm(seur)
  
  assigns = meta[,c("CellType","orig.ident")]
  real_mi = dmi(assigns)$mi[1,2]
  miR = data.frame("age" = ages[[j]], "dataset" = 1, "mi" = real_mi, "cellline"="Mito210", tech="ATAC")
  mis = rbind(mis, miR)
}

ggplot(mis, aes(age, mi, color = tech)) +
  geom_jitter(size=3, alpha=0.6, height=0, width=0.2) +
  theme_classic() +
  geom_hline(yintercept=0.008, linetype="dashed", color="gray72") +
  geom_hline(yintercept=0.064, linetype="dashed", color="gray72") +
  scale_color_manual(values=c("#4DAF4A", "gray"))
ggsave("mutual_information-atacAndRNA.pdf", width=4, height=3)

ggplot(mis[mis$tech=="ATAC",], aes(age, mi, color = cellline)) +
  geom_point(size=3) +
  theme_classic() +
  geom_hline(yintercept=0.008, linetype="dashed", color="gray72") +
  geom_hline(yintercept=0.064, linetype="dashed", color="gray72") +
  scale_color_manual(values=c("#4DAF4A"))
ggsave("mutual_information-atac.pdf", width=4, height=3)


saveRDS(mis, "mutualinformations.rds")
