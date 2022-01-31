library(Seurat)
library(harmony)
library(ggplot2)
library(DropletUtils)
library(factoextra)
library(bioDist)
library(dplyr)
library(reshape2)
library(psych)


#Export for Compass
combined = readRDS("../2Human+3m_met/2human.all3m.harmonized.rds")
Idents(combined) = "dataset"
seur = subset(combined, downsample=1000)
seur = NormalizeData(seur, normalization.method = "RC", scale.factor = 1e6) #CPM - compass needs linear not LogNormalized
write10xCounts(path="2humanData-down1k-CPM/", x=GetAssayData(seur,slot="data"))
write.table(rownames(GetAssayData(seur,slot="data")), "2humanData-down1k-CPM/genes.tsv", quote = F, row.names = F, col.names = F)
saveRDS(seur, "2humanData-down1k-CPM/seur.rds")
#conda activate compass
#compass --data-mtx matrix.mtx genes.tsv barcodes.tsv --num-processes 10 --species homo_sapiens

setwd("2humanData-down1k-CPM/")
seur = readRDS("seur.rds")

#metareactions and scale penalties
penalties = read.table("reactions.tsv", sep="\t", header = T, row.names = 1)
penalties[penalties <= 1e-4] = 0
penalties = penalties[rowSums(penalties)!=0,]
penalties = penalties[rowMaxs(as.matrix(penalties)) - rowMins(as.matrix(penalties)) >= 1e-3,]
dist = spearman.dist(as.matrix(penalties))
hc <- hclust(dist)
clus = cutree(hc, h=0.05)
penalties$meta = clus
meta  = group_by(penalties, meta) %>% summarize_all(mean)
mat = as.matrix(-1 * log1p(meta[,2:ncol(meta)]))
mat = mat - min(mat)

groups = seur$dataset[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$dataset[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
res.pca <- prcomp(t(mat))
saveRDS(res.pca, "res.pca.scaledMetaRxns.rds")
fviz_pca_ind(res.pca,habillage=groups,geom="point", palette=c("#DDAA33","#BB5566","#004488"), title="Compass Score PCA", pointshape=19)
ggsave("compass-pca-dataset-scaledMetaRxns.tiff")
loadings <- res.pca$rotation
top = order(abs(loadings[,1]), decreasing = T)[1:25]
paths = read.table("../rxn_md.csv", sep=",", quote='"', header=T)
toppaths = rownames(penalties[penalties$meta %in% top,])
write.table(paths[paste0(paths$rxn_code_nodirection, "_pos") %in% toppaths | paste0(paths$rxn_code_nodirection, "_neg") %in% toppaths,], 
           "top_PC1_pathways-scaledMetaRxns.txt", sep="\t")

seur$FinalName[is.na(seur$FinalName)] = seur$CellSubTypes[is.na(seur$FinalName)]
seur$FinalName[is.na(seur$FinalName)] = seur$Cluster[is.na(seur$FinalName)]
seur$MatchedCellType = NA
seur$MatchedCellType[seur$FinalName %in% c("GluN1","GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8","CPN", "ExN","ExM","ExM-U")] = "UL Neurons"
seur$MatchedCellType[seur$FinalName %in% c("CFuPN","SP","ExDp", "ExDp1","ExDp2")] = "DL Neurons"
seur$MatchedCellType[seur$FinalName %in% c("IP","nIPC")] = "IPCs"
seur$MatchedCellType[seur$FinalName %in% c("CGE_IN","MGE_IN","Immature IN","InMGE","InCGE")] = "INs"
seur$MatchedCellType[seur$FinalName %in% c("tRG","Early_RG","Late_RG","Cyc_Prog","aRG","oRG","PgS","PgG2M", "vRG")] = "RG"
seur$MatchedCellType[seur$FinalName=="PN"] = "Uns. PNs"
seur$MatchedCellType[is.na(seur$MatchedCellType)] = "Other"
groups = seur$MatchedCellType[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$MatchedCellType[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
fviz_pca_ind(res.pca,habillage=groups,geom="point", title="Compass Score PCA Cell Type", pointshape=19, axes=2:3)
ggsave("compass-pca-celltype-scaledMetaRxns-dim23.png")

groups = seur$nCount_RNA[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$nCount_RNA[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
fviz_pca_ind(res.pca,col.ind=log(groups), geom="point", pointshape=19)
ggsave("compass-pca-nUMI-scaledMetaRxns.png")

#correlate PCs with things
cell_PC_scores = as.data.frame(res.pca$x)[,1:20]
groups = seur$dataset[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$dataset[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
cell_PC_scores$dataset=groups
groups = seur$MatchedCellType[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$MatchedCellType[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
cell_PC_scores$CellType = groups
groups = seur$nCount_RNA[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
groups[is.na(groups)] = seur$nCount_RNA[colnames(mat)[is.na(groups)]] #CCAGCACGCCTG.x and CGTATGAATGCA.y
cell_PC_scores$nCount = groups

rows = c("nCount","dataset","CellType")
columns = paste0("PC",1:20)
rsqs = data.frame()
for (r in rows) {
  for (c in columns) {
    model.lm = lm(substitute(c ~ r, list(c = as.name(c), r = as.name(r))), data = cell_PC_scores)
    rsq = summary(model.lm)$r.squared
    rsqs[r,c] = rsq
  }
}
rsqs$var = rownames(rsqs)
mrsq = melt(rsqs)
colnames(mrsq)[3] = "Rsquared"
ggplot(mrsq, aes(x=variable, y=var, fill=Rsquared)) +
  geom_tile(color="white", size=0.25) +
  labs(x="",y="")+
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0,0), position="top")+
  scale_fill_gradient2(low="lightgray",mid="lightgray", high="mediumorchid4", midpoint=0.02)+
  theme(axis.ticks=element_blank())
ggsave("PCsCorrelationsMetaData.pdf", width=9, height=2)




#trevino- Normalize only metgenes and similar cell type proportions
seur = readRDS("../../2Human+3m_met/2human.all3m.harmonized.rds")
seur = subset(seur, subset=dataset!="Geschwind2019")
seur$FinalName[is.na(seur$FinalName)] = seur$CellSubTypes[is.na(seur$FinalName)]
seur$MatchedCellType = NA
seur$MatchedCellType[seur$FinalName %in% c("GluN1","GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8","CPN", "ExN","ExM","ExM-U")] = "UL Neurons"
seur$MatchedCellType[seur$FinalName %in% c("CFuPN","SP","ExDp", "ExDp1","ExDp2")] = "DL Neurons"
seur$MatchedCellType[seur$FinalName %in% c("IP","nIPC")] = "IPCs"
seur$MatchedCellType[seur$FinalName %in% c("CGE_IN","MGE_IN","Immature IN","InMGE","InCGE")] = "INs"
seur$MatchedCellType[seur$FinalName %in% c("tRG","Early_RG","Late_RG","Cyc_Prog","aRG","oRG","PgS","PgG2M", "vRG")] = "RG"
seur$MatchedCellType[seur$FinalName=="PN"] = "Uns. PNs"
seur$MatchedCellType[is.na(seur$MatchedCellType)] = "Other"
table(seur$MatchedCellType, seur$dataset)
seur$ct_ds = paste0(seur$MatchedCellType, "_", seur$dataset)
Idents(seur) = "ct_ds"
seur = subset(seur, downsample=50)
#compass --list-genes "metGenes.txt" --species homo_sapiens
metgenes = as.character(read.table("metGenes.txt")$V1)
seur = subset(seur, features=metgenes)
seur = subset(seur, subset=nCount_RNA<2500) #removing 3 outliers
seur = NormalizeData(seur, normalization.method = "RC", scale.factor = 1e6) #CPM - compass needs linear not LogNormalized
write10xCounts(path="trevino3m-celltypes-down50-CPM-metgenes/", x=GetAssayData(seur,slot="data"))
write.table(rownames(GetAssayData(seur,slot="data")), "trevino3m-celltypes-down50-CPM-metgenes/genes.tsv", quote = F, row.names = F, col.names = F)
saveRDS(seur, "trevino3m-celltypes-down50-CPM-metgenes/seur.rds")
#conda activate compass
#compass --data-mtx matrix.mtx genes.tsv barcodes.tsv --num-processes 15 --species homo_sapiens

setwd("trevino3m-celltypes-down50-CPM-metgenes/")
seur = readRDS("seur.rds")

#Differential rxns
#scale penalties, no metarxns
penalties2 = read.table("reactions.tsv", sep="\t", header = T, row.names = 1)
mat = as.matrix(-1 * log1p(penalties2))
mat = mat[rowMaxs(mat) - rowMins(mat) >= 1e-3,]
mat = mat - min(mat)

groups = seur$dataset[ifelse(startsWith(colnames(mat), "X"), substr(gsub("\\.","-",colnames(mat)), 2, nchar(colnames(mat))), gsub("\\.","-",colnames(mat)))]
orgCells = mat[,groups=="organoid_3m"]
humanCells = mat[,groups=="TrevinoRNA"]
paths = read.table("../rxn_md.csv", sep=",", quote='"', header=T)

wilcox_results = data.frame()
for (i in 1:nrow(orgCells)) {
  rxn = rownames(orgCells)[[i]]
  wilcox_results[rxn,"rxn"] = rxn
  wilcox_results[rxn,"rxn_sys"] = as.character(paths[paste0(paths$rxn_code_nodirection, "_pos") ==rxn | paste0(paths$rxn_code_nodirection, "_neg") ==rxn,"subsystem"])
  wilcox_results[rxn,"rxn_name_long"] = as.character(paths[paste0(paths$rxn_code_nodirection, "_pos") ==rxn | paste0(paths$rxn_code_nodirection, "_neg") ==rxn,"rxn_name_long"])
  res = wilcox.test(orgCells[i,], humanCells[i,], paired=F, alternative = "two.sided")
  wilcox_results[rxn,"wilcox_p"] = res$p.value
  wilcox_results[rxn,"wilcox_stat"] = res$statistic
  wilcox_results[rxn,"cohen_d"] = cohen.d(mat[i,], group=groups=="organoid_3m")$cohen.d[,2]
}
wilcox_results$adj_p = p.adjust(wilcox_results$wilcox_p, method="BY")

ggplot(wilcox_results, 
       aes(x=cohen_d, y=-1*log10(adj_p)))+
  geom_point()+
  labs(y="-log(Wilcoxon adjusted p)", x="Cohen's d" ) +
  facet_wrap(vars(rxn_sys))
ggsave("volcanoPlot-all-scaled-systems.png", width=15, height=10)

library(ggrepel)
wilcox_results$rxn_sys_short = "Other" 
wilcox_results$rxn_sys_short[wilcox_results$rxn_sys=="Glycolysis/gluconeogenesis"] = "Glycolysis"
wilcox_results$rxn_sys_short[wilcox_results$rxn_sys %in% c("Transport, extracellular", "Transport, lysosomal", "Transport, endoplasmic reticular","Transport, nuclear","Transport, mitochondrial","Transport, golgi apparatus")] = "Cellular Transport"
wilcox_results$rxn_sys_short[wilcox_results$rxn_sys=="Oxidative phosphorylation"] = "Oxidative Phosphorylation"
wilcox_results$rxn_sys_short[wilcox_results$rxn_sys=="Fructose and mannose metabolism"] = "Fructose and mannose metabolism"
wilcox_results$rxn_sys_short[wilcox_results$rxn_sys=="Pentose phosphate pathway"] = "Pentose phosphate pathway"

pdf("volcanoPlot-scaled-select.pdf", width=7, height=5)
ggplot(wilcox_results, aes(x=cohen_d, y=-1*log10(adj_p)))+
  geom_point(aes(col=rxn_sys_short), alpha=0.8)+
  geom_point(data=subset(wilcox_results, rxn_sys=="Oxidative phosphorylation"), aes(x=cohen_d, y=-1*log10(adj_p)), col="#88ccee", alpha=0.8)+
  geom_point(data=subset(wilcox_results, rxn_sys=="Glycolysis/gluconeogenesis"), aes(x=cohen_d, y=-1*log10(adj_p)), col="#332288", alpha=0.8)+
  labs(y="-log(Wilcoxon adjusted p)", x="Cohen's d" ) +
  #geom_label_repel(aes(label=ifelse((-1*log10(adj_p)>15) | (cohen_d<(-0.55)),rxn_name_long,"")),max.overlaps = 50)+
  geom_label_repel(aes(label=ifelse(rxn %in% c("ALAt4_pos","PFK_pos"),rxn_name_long,""))) +
  scale_color_manual(values=c("#ddcc77","#44aa99","#332288","#dddddd","#88ccee","#aa4499")) +
  theme_classic()
dev.off()