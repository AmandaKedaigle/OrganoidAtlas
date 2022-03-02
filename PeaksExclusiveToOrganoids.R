## Get GO terms associated with the peaks called in organoids/tissue exclusively:
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(ggplot2)

suv1Peaks <- readRDS("SUV1_CellTypePeaks.RDS")
suv1Ranges <- readRDS("SUV1_CellTypeRanges.RDS")
for(celltype in names(suv1Peaks)){
  suv1Ranges[[celltype]]$Gene <- suv1Peaks[[celltype]]$gene_name
}
suv3Peaks <- readRDS("SUV3_CellTypePeaks.RDS")
suv3Ranges <- readRDS("SUV3_CellTypeRanges.RDS")
for(celltype in names(suv3Peaks)){
  suv3Ranges[[celltype]]$Gene <- suv3Peaks[[celltype]]$gene_name
}
suv6Peaks <- readRDS("SUV6_CellTypePeaks.RDS")
suv6Ranges <- readRDS("SUV6_CellTypeRanges.RDS")
for(celltype in names(suv6Peaks)){
  suv6Ranges[[celltype]]$Gene <- suv6Peaks[[celltype]]$gene_name
}
TrevPeaks <- readRDS("Trevino_CellTypePeaks.RDS")
TrevRanges <- readRDS("Trevino_CellTypeRanges.RDS")
for(celltype in names(TrevPeaks)){
  TrevRanges[[celltype]]$Gene <- TrevPeaks[[celltype]]$gene_name
}
SubPeaks <- readRDS("Trevino_SubTypePeaks.RDS")
SubRanges <- readRDS("Trevino_SubTypeRanges.RDS")
for(celltype in names(SubPeaks)){
  SubRanges[[celltype]]$Gene <- SubPeaks[[celltype]]$gene_name
}

arg1mo <- suv1Ranges$aRG[! suv1Ranges$aRG %over% TrevRanges$glia,] #Compare the 1-month aRGs to Trevino glia; what's accessible in organoids only?
arg3mo <- suv3Ranges$oRG[! suv3Ranges$oRG %over% TrevRanges$glia,] #Compare the 1-month aRGs to Trevino glia; what's accessible in organoids only?
arg6mo <- suv6Ranges$aRG[! suv6Ranges$aRG %over% TrevRanges$glia,] #Compare the 1-month aRGs to Trevino glia; what's accessible in organoids only?

arg1moGenes <- unique(arg1mo$Gene) #8952 genes
arg3moGenes <- unique(arg3mo$Gene) #7026 genes
arg6moGenes <- unique(arg6mo$Gene) #6736 genes

## What about exclusive to Trevino?
gliTrev <- TrevRanges$glia[! TrevRanges$glia %over% suv1Ranges$aRG,]
gliTrev <- gliTrev[! gliTrev %over% suv3Ranges$oRG,]
gliTrev <- gliTrev[! gliTrev %over% suv6Ranges$aRG,]
gliTrevGenes <- unique(gliTrev$Gene)

arg1moExGenes <- arg1moGenes[! arg1moGenes %in% gliTrevGenes] 
arg3moExGenes <- arg3moGenes[! arg3moGenes %in% gliTrevGenes]
arg6moExGenes <- arg6moGenes[! arg6moGenes %in% gliTrevGenes]
sum(arg1moExGenes %in% arg3moExGenes & arg1moExGenes %in% arg6moExGenes)
aRGorg <- arg1moExGenes[arg1moExGenes %in% arg3moExGenes & arg1moExGenes %in% arg6moExGenes]

arg1moGO <- enrichGO(gene = arg1moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
arg3moGO <- enrichGO(gene = arg3moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
arg6moGO <- enrichGO(gene = arg6moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
aRGorgGO <- enrichGO(gene = aRGorg, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
dotplot(aRGorgGO) 

## How about CPNs?
pn1mo <- suv1Ranges$`Newborn PN`[! suv1Ranges$`Newborn PN` %over% TrevRanges$exNeuron,] #Compare the 1-month PNs to Trevino ; what's accessible in organoids only?
cpn3mo <- suv3Ranges$CPN[! suv3Ranges$CPN %over% TrevRanges$exNeuron,] #Compare the 1-month CPNs to Trevino; what's accessible in organoids only?
cpn6mo <- suv6Ranges$CPN[! suv6Ranges$CPN %over% TrevRanges$exNeuron,] #Compare the 1-month CPNs to Trevino; what's accessible in organoids only?
pn1moGenes <- unique(pn1mo$Gene) 
cpn3moGenes <- unique(cpn3mo$Gene) 
cpn6moGenes <- unique(cpn6mo$Gene) 
## What about exclusive to Trevino?
exTrev <- TrevRanges$exNeuron[! TrevRanges$exNeuron %over% suv1Ranges$`Newborn PN`,]
exTrev <- exTrev[! exTrev %over% suv3Ranges$CPN,]
exTrev <- exTrev[! exTrev %over% suv6Ranges$CPN,]
exTrevGenes <- unique(exTrev$Gene) 

pn1moExGenes <- pn1moGenes[! pn1moGenes %in% exTrevGenes] 
cpn3moExGenes <- cpn3moGenes[! cpn3moGenes %in% exTrevGenes] 
cpn6moExGenes <- cpn6moGenes[! cpn6moGenes %in% exTrevGenes] 
sum(pn1moExGenes %in% cpn3moExGenes & pn1moExGenes %in% cpn6moExGenes) # genes in common between all three organoids, absent from Trevino exNeurons
CPNorg <- pn1moExGenes[pn1moExGenes %in% cpn3moExGenes & pn1moExGenes %in% cpn6moExGenes]

pn1moGO <- enrichGO(gene = pn1moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
cpn3moGO <- enrichGO(gene = cpn3moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
cpn6moGO <- enrichGO(gene = cpn6moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
CPNorgGO <- enrichGO(gene = CPNorg, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
dotplot(CPNorgGO) +
  ggtitle("GO terms enriched in genes accessible in organoid CPNs\nbut not fetal exNeurons") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #

## How about IPs?
ip1mo <- suv1Ranges$IP[! suv1Ranges$IP %over% TrevRanges$IPC,] #Compare the 1-month IPs to Trevino ; what's accessible in organoids only?
ip3mo <- suv3Ranges$IP[! suv3Ranges$IP %over% TrevRanges$IPC,] #Compare the 1-month IPs to Trevino; what's accessible in organoids only?
ip6mo <- suv6Ranges$IP[! suv6Ranges$IP %over% TrevRanges$IPC,] #Compare the 1-month IPs to Trevino; what's accessible in organoids only?
ip1moGenes <- unique(ip1mo$Gene)
ip3moGenes <- unique(ip3mo$Gene) 
ip6moGenes <- unique(ip6mo$Gene)
## What about exclusive to Trevino?
ipTrev <- TrevRanges$IPC[! TrevRanges$IPC %over% suv1Ranges$IP,]
ipTrev <- ipTrev[! ipTrev %over% suv3Ranges$IP,]
ipTrev <- ipTrev[! ipTrev %over% suv6Ranges$IP,] 
ipTrevGenes <- unique(ipTrev$Gene) 

ip1moExGenes <- ip1moGenes[! ip1moGenes %in% ipTrevGenes] 
ip3moExGenes <- ip3moGenes[! ip3moGenes %in% ipTrevGenes] 
ip6moExGenes <- ip6moGenes[! ip6moGenes %in% ipTrevGenes]
sum(ip1moExGenes %in% ip3moExGenes & ip1moExGenes %in% ip6moExGenes) # genes in common between all three organoids, absent from Trevino exNeurons
IPorg <- ip1moExGenes[ip1moExGenes %in% ip3moExGenes & ip1moExGenes %in% ip6moExGenes]

ip1moGO <- enrichGO(gene = ip1moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
ip3moGO <- enrichGO(gene = ip3moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
ip6moGO <- enrichGO(gene = ip6moExGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
IPorgGO <- enrichGO(gene = IPorg, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
dotplot(IPorgGO) +
  ggtitle("GO terms enriched in genes accessible in organoid IP\nbut not fetal IPC") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #

## How about INs?
in6mo <- suv6Ranges$`Immature IN`[! suv6Ranges$`Immature IN` %over% SubRanges$IN2,] #Compare 6-mo INs to Trevino IN2, specifically
doubIN <- names(table(in6mo$Gene)[table(in6mo$Gene)>1]) # genes have more than 1 peak
inGenes <- doubIN[!doubIN %in% SubRanges$IN2$Gene] # genes with multiple peaks in organoids and no peaks in fetal
inGenes2 <- doubIN[!doubIN %in% TrevRanges$intNeuron$Gene]# genes not found in *any* Trevino intNeur

INorgGO <- enrichGO(gene = inGenes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
INorgGO2 <- simplify(INorgGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(INorgGO2) +
  ggtitle("GO terms enriched in genes accessible\nin organoid Immature IN but not fetal IN2") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #

## How 6mo CPN vs GluN6
CPN6mo_GluN6 <- suv6Ranges$CPN[! suv6Ranges$CPN %over% SubRanges$GluN6,] #Compare 6-mo INs to Trevino GluN6
CPN6mo_GluN6_Strong <- names(table(CPN6mo_GluN6$Gene)[table(CPN6mo_GluN6$Gene)>1]) #3936 genes have more than 1 peak
CPN6mo_GluN6_Genes <- CPN6mo_GluN6_Strong[! CPN6mo_GluN6_Strong %in% SubRanges$GluN6$Gene] #525 genes exclusively accessible in 6-mo organoids

CPN6moGO <- enrichGO(gene = CPN6mo_GluN6_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
CPN6moGO2 <- simplify(CPN6moGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(CPN6moGO2) +
  ggtitle("GO terms enriched in genes accessible\nin 6 month CPNs but not fetal GluN6") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
GluN6_CPN6 <- SubRanges$GluN6[! SubRanges$GluN6 %over% suv6Ranges$CPN,] #Compare 6-mo INs to Trevino GluN6
GluN6_CPN6_Strong <- names(table(GluN6_CPN6$Gene)[table(GluN6_CPN6$Gene)>1]) #15221 genes have more than 1 peak
GluN6_CPN6_Genes <- GluN6_CPN6_Strong[! GluN6_CPN6_Strong %in% suv6Ranges$CPN$Gene] #974 genes exclusively accessible in fetal GluN6

GluN6GO <- enrichGO(gene = GluN6_CPN6_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
GluN6GO2 <- simplify(GluN6GO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(GluN6GO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal GluN6 but not 6 month CPNs") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))
## IP to nIPC
IP6mo_IPC <- suv6Ranges$IP[! suv6Ranges$IP %over% SubRanges$nIPC,] #Compare 6-mo IPs to Trevino nIPC
IP6mo_IPC_Strong <- names(table(IP6mo_IPC$Gene)[table(IP6mo_IPC$Gene)>1]) #3936 genes have more than 1 peak
IP6mo_IPC_Genes <- IP6mo_IPC_Strong[! IP6mo_IPC_Strong %in% SubRanges$nIPC$Gene] #512 genes exclusively accessible in 6-mo organoids
IP6moGO <- enrichGO(gene = IP6mo_IPC_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
IP6moGO2 <- simplify(IP6moGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(IP6moGO2) +
  ggtitle("GO terms enriched in genes accessible\nin 6 month IP but not fetal IPC") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
IPC_6moIP <- SubRanges$nIPC[! SubRanges$nIPC %over% suv6Ranges$IP,] #Compare fetal nIPCs to IPs
IPC_6moIP_Strong <- names(table(IPC_6moIP$Gene)[table(IPC_6moIP$Gene)>1]) #18737 genes have more than 1 peak
IPC_6moIP_Genes <- IPC_6moIP_Strong[! IPC_6moIP_Strong %in% suv6Ranges$IP$Gene] #1615 genes exclusively accessible in fetal nIPCs
IPCGO <- enrichGO(gene = IPC_6moIP_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
IPCGO2 <- simplify(IPCGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(IPCGO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal IPC but not 6 month  organoid IP") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))
## Immature IN to IN2
IN6mo_IN2 <- suv6Ranges$`Immature IN`[! suv6Ranges$`Immature IN` %over% SubRanges$IN2,] #Compare 6-mo IPs to Trevino nIPC
IN6mo_IN2_Strong <- names(table(IN6mo_IN2$Gene)[table(IN6mo_IN2$Gene)>1]) #3936 genes have more than 1 peak
IN6mo_IN2_Genes <- IN6mo_IN2_Strong[! IN6mo_IN2_Strong %in% SubRanges$IN2$Gene] #512 genes exclusively accessible in 6-mo organoids
IN6moGO <- enrichGO(gene = IN6mo_IN2_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
IN6moGO2 <- simplify(IN6moGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(IP6moGO2) +
  ggtitle("GO terms enriched in genes accessible\nin 6 month Immature IN but not fetal IN2") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
IN2_IN6mo <- SubRanges$IN2[! SubRanges$IN2 %over% suv6Ranges$`Immature IN`,] #Compare fetal nIPCs to IPs
IN2_IN6mo_Strong <- names(table(IN2_IN6mo$Gene)[table(IN2_IN6mo$Gene)>1]) #18737 genes have more than 1 peak
IN2_IN6mo_Genes <- IN2_IN6mo_Strong[! IN2_IN6mo_Strong %in% suv6Ranges$`Immature IN`$Gene] #1615 genes exclusively accessible in fetal nIPCs
IN2GO <- enrichGO(gene = IN2_IN6mo_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                  minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
IN2GO2 <- simplify(IN2GO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(IN2GO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal IN2 but not 6 month organoid Immature IN") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))
## CFuPN (3-month) to GluN8
CFuPN3mo_GluN8 <- suv3Ranges$CFuPN[! suv3Ranges$CFuPN %over% SubRanges$GluN8,] #Compare 6-mo IPs to Trevino nIPC
CFuPN3mo_GluN8_Strong <- names(table(CFuPN3mo_GluN8$Gene)[table(CFuPN3mo_GluN8$Gene)>1]) #3936 genes have more than 1 peak
CFuPN3mo_GluN8_Genes <- CFuPN3mo_GluN8_Strong[! CFuPN3mo_GluN8_Strong %in% SubRanges$GluN8$Gene] #512 genes exclusively accessible in 6-mo organoids
CFuPN3moGO <- enrichGO(gene = CFuPN3mo_GluN8_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
CFuPN3moGO2 <- simplify(CFuPN3moGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(CFuPN3moGO2) +
  ggtitle("GO terms enriched in genes accessible\nin 3 month CFuPN but not fetal GluN8") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
GluN8_CFuPN3 <- SubRanges$GluN8[! SubRanges$GluN8 %over% suv3Ranges$CFuPN,] #Compare fetal nIPCs to IPs
GluN8_CFuPN3_Strong <- names(table(GluN8_CFuPN3$Gene)[table(GluN8_CFuPN3$Gene)>1]) #18737 genes have more than 1 peak
GluN8_CFuPN3_Genes <- GluN8_CFuPN3_Strong[! GluN8_CFuPN3_Strong %in% suv3Ranges$CFuPN$Gene] #1615 genes exclusively accessible in fetal nIPCs
GluN8GO <- enrichGO(gene = GluN8_CFuPN3_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                  minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
GluN8GO2 <- simplify(GluN8GO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(GluN8GO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal GluN8 but not 3 month organoid CFuPN") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))

aRG1mo_Early <- suv1Ranges$aRG[! suv1Ranges$aRG %over% SubRanges$`Early RG`,]
aRG1mo_Early_Strong <- names(table(aRG1mo_Early$Gene)[table(aRG1mo_Early$Gene)>1]) #3936 genes have more than 1 peak
aRG1mo_Early_Genes <- aRG1mo_Early_Strong[! aRG1mo_Early_Strong %in% SubRanges$`Early RG`$Gene] #512 genes exclusively accessible in 6-mo organoids
aRG1moGO <- enrichGO(gene = aRG1mo_Early_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                       minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
aRG1moGO2 <- simplify(aRG1moGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(aRG1moGO2) +
  ggtitle("GO terms enriched in genes accessible\nin 1 month aRG but not fetal Early RG") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
Early_aRG1 <- SubRanges$`Early RG`[! SubRanges$`Early RG` %over% suv1Ranges$aRG,] #Compare fetal nIPCs to IPs
Early_aRG1_Strong <- names(table(Early_aRG1$Gene)[table(Early_aRG1$Gene)>1]) #18737 genes have more than 1 peak
Early_aRG1_Genes <- Early_aRG1_Strong[! Early_aRG1_Strong %in% suv1Ranges$aRG$Gene] #1615 genes exclusively accessible in fetal nIPCs
EarlyGO <- enrichGO(gene = Early_aRG1_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
EarlyGO2 <- simplify(EarlyGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(EarlyGO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal Early RG but not 1 month organoid aRG") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))

#Unspecified PNs
PN3_exNeur <- suv3Ranges$PN[! suv3Ranges$PN %over% TrevRanges$exNeuron,]
PN6_exNeur <- suv6Ranges$PN[! suv6Ranges$PN %over% TrevRanges$exNeuron,]
PN3_Strong <- names(table(PN3_exNeur$Gene)[table(PN3_exNeur$Gene)>1])
PN6_Strong <- names(table(PN6_exNeur$Gene)[table(PN6_exNeur$Gene)>1])
PN_Strong <- unique(c(PN3_Strong, PN6_Strong))
PN_Genes <- PN_Strong[! PN_Strong %in% TrevRanges$exNeuron$Gene]
PNGO <- enrichGO(gene = PN_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                     minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
PNGO2 <- simplify(PNGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(PNGO2) +
  ggtitle("GO terms enriched in genes accessible\nin organoid PNs but not fetal exNeurons") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) #
## And reverse it
Early_aRG1 <- SubRanges$`Early RG`[! SubRanges$`Early RG` %over% suv1Ranges$aRG,] #Compare fetal nIPCs to IPs
Early_aRG1_Strong <- names(table(Early_aRG1$Gene)[table(Early_aRG1$Gene)>1]) #18737 genes have more than 1 peak
Early_aRG1_Genes <- Early_aRG1_Strong[! Early_aRG1_Strong %in% suv1Ranges$aRG$Gene] #1615 genes exclusively accessible in fetal nIPCs
EarlyGO <- enrichGO(gene = Early_aRG1_Genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP",
                    minGSSize = 10,  pvalueCutoff = 0.05, qvalueCutoff = 0.15)
EarlyGO2 <- simplify(EarlyGO, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(EarlyGO2) +
  ggtitle("GO terms enriched in genes accessible\nin fetal Early RG but not 1 month organoid aRG") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14))


IN6mo <- IN6moGO@result
IN6mo$Comparison <- "Org.\nIm. INs"
IN2 <- IN2GO@result
IN2$Comparison <- "Fetal\nINs"
CPN6mo <- CPN6moGO@result
CPN6mo$Comparison <- "Org.\nCPNs"
GluN6 <- GluN6GO@result
GluN6$Comparison <- "Fetal\nGluN6"
IP6mo <- IP6moGO@result
IP6mo$Comparison <- "Org.\nIPs"
IPC <- IPCGO@result
IPC$Comparison <- "Fetal\nIPCs"
CFuPN <- CFuPN3moGO@result
CFuPN$Comparison <- "Org.\nCFuPNs"
GluN8 <- GluN8GO@result
GluN8$Comparison <- "Fetal\nGluN8"
aRG <- aRG1moGO@result
aRG$Comparison <- "Org.\naRGs"
EarlyRG <- EarlyGO@result
EarlyRG$Comparison <- "Fetal\nEarly RGs"

fullRes <- rbind(IN6mo, IN2, CPN6mo, GluN6, IP6mo, IPC, CFuPN, GluN8, aRG, EarlyRG)
SigTerms <- unique(c(IN6moGO2@result$ID[1:8], IN2GO2@result$ID[1:8], CPN6moGO2@result$ID[1:8], GluN6GO2@result$ID[1:8], IP6moGO2@result$ID[1:8],
                   IPCGO2@result$ID[1:8], CFuPN3moGO2@result$ID[1:8], GluN8GO2@result$ID[1:8], aRG1moGO2@result$ID[1:8], EarlyGO2@result$ID[1:8]))
SigTerms <- SigTerms[!is.na(SigTerms)]
ggplot(fullRes[fullRes$ID %in% SigTerms,], aes(x=Comparison,
                                               y=factor(Description,
                                                        levels=levels(factor(Description))[rev(c(26,20,19,21,17,7,9,8,22,23,1,2,6,5,18,4,3,10,12,13,14,15,16,11,24,25))]),
                                               size=-log10(p.adjust), fill=Count, color=(p.adjust < 0.05))) +
  geom_point(position="identity", pch=21) +
  scale_fill_gradient(low="white", high="darkred") +
  scale_color_manual(values=c("grey","black")) +
  ggtitle("GO enrichment in genes exclusively accessible in\norganoids or fetal tissue") +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_blank()) +
  xlab("Cell type")


genes1 <- readRDS("GenesFromPeaks_up_in_1m_not_Trev.RDS")

myList1 <- unique(genes1$gene_name)

myGO1 <- enrichGO(gene = myList1,
         OrgDb = org.Hs.eg.db, 
         keyType = 'SYMBOL',
         ont = "BP",
         minGSSize = 10,
         pvalueCutoff = 0.05, 
         qvalueCutoff = 0.15)

simGO1 <- simplify(myGO1,cutoff=0.5, by="p.adjust", select_fun=min)

allRes <- simGO@result
allRes$Age <- "1mo"

genes3 <- readRDS("GenesFromPeaks_up_in_3m_not_Trev.RDS")
myList3 <- unique(genes3$gene_name)
myGO3 <- enrichGO(gene = myList3,
                 OrgDb = org.Hs.eg.db, 
                 keyType = 'SYMBOL',
                 ont = "BP",
                 minGSSize = 10,
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.15)
simGO3 <- simplify(myGO3,cutoff=0.5, by="p.adjust", select_fun=min)
res3 <- simGO3@result
res3$Age <- "3mo"
allRes <- rbind(allRes, res3)

genes6 <- readRDS("GenesFromPeaks_up_in_6m_not_Trev.RDS")
myList6 <- unique(genes6$gene_name)
myGO6 <- enrichGO(gene = myList6,
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'SYMBOL',
                  ont = "BP",
                  minGSSize = 10,
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.15)
simGO6 <- simplify(myGO6,cutoff=0.5, by="p.adjust", select_fun=min)
res6 <- simGO6@result
res6$Age <- "6mo"
allRes <- rbind(allRes, res6)

ggplot(allRes, aes(x=Age, y=Description, size=Count, color=qvalue)) +
  geom_point(stat="identity")


TrevGenes <- readRDS("GenesFromPeaks_Tev_No_Orgs.RDS")

TrevGenes <- unique(TrevGenes$gene_name)
genes1Unique <- myList1[!myList1 %in% TrevGenes]
genes3Unique <- myList3[!myList3 %in% TrevGenes]
genes6Unique <- myList6[!myList6 %in% TrevGenes]

unGO1 <- enrichGO(gene = genes1Unique,
                  OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", minGSSize = 20, maxGSSize = 200, pvalueCutoff = 0.05, qvalueCutoff = 0.15)
simUnGO1 <- simplify(unGO1,cutoff=0.5, by="p.adjust", select_fun=min)
unRes <- simUnGO1@result
unRes$Age <- "1mo"

unGO3 <- enrichGO(gene = genes3Unique,
                  OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", minGSSize = 10, pvalueCutoff = 0.05, qvalueCutoff = 0.15)
simUnGO3 <- simplify(unGO3,cutoff=0.5, by="p.adjust", select_fun=min)
unRes3 <- simUnGO3@result
unRes3$Age <- "3mo"
unRes <- rbind(unRes, unRes3)

unGO6 <- enrichGO(gene = genes6Unique,
                  OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP", minGSSize = 10, pvalueCutoff = 0.05, qvalueCutoff = 0.15)
simUnGO6 <- simplify(unGO6,cutoff=0.5, by="p.adjust", select_fun=min)
unRes6 <- simUnGO6@result
unRes6$Age <- "6mo"
unRes <- rbind(unRes, unRes6)
