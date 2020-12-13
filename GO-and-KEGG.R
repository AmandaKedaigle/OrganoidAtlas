##Gene Ontology

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)


#for xlsx, sort by pval threshold:
degs`23` = read_xlsx("file1.xlsx") %>% 
  filter(log2FoldChange > 1 ) %>% 
  top_n(n=-500, wt = padj)
degs1$time = 'a'

degs1 <- read_xlsx("file2.xlsx") %>% 
  filter(log2FoldChange > 1 ) %>% 
  #group_split(Cluster)
  #group_by(Cluster)  %>% 
  top_n(n=-500, wt = padj) 
degs2$time = 'b'

degs3 <- read_xlsx("file3.xlsx") %>% 
  filter(log2FoldChange > 1 ) %>% 
  #group_split(Cluster)
  #group_by(Cluster)  %>% 
  top_n(n=-500, wt = padj) 
degs3$time = 'c'

degs4 <- read_xlsx("file4.xlsx") %>% 
  filter(log2FoldChange > 1 ) %>% 
  #group_split(Cluster)
  #group_by(Cluster)  %>% 
  top_n(n=-500, wt = padj) 
degs4$time = 'd'

degs <- rbind(degs1, degs2, degs3, degs4)  %>% group_split(time)


genes = list()
for (time in 1:length(degs)) {
  genes[[time]] = pull(degs[[time]], "gene") 
}
names(genes) = 0:(length(degs)-1)


#for basic GO

ck = compareCluster(genes, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
#to add a background, make a list of genes, and then add "universe = backgroundlist" above
ck2 = clusterProfiler::simplify(ck,cutoff=0.7, by="p.adjust", select_fun=min)
ck2b <- ck2
levels <- names(genes)[order(names(genes))]
ck2b@compareClusterResult$Cluster <- factor(ck2b@compareClusterResult$Cluster, levels = levels)

finalplot <- dotplot(ck2b, showCategory = 4) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_colour_gradient2(low = "chose", mid = "your", high = "colors")

pdf('filename.pdf', width = 8, height=4.5)
print(finalplot)
dev.off()


#for KEGG, we need entrez ids
genesEntrez = list()
for (l in names(genes)) {
  genesEntrez[[l]] = mapIds(org.Hs.eg.db,
                            keys=genes[[l]],
                            keytype="SYMBOL",
                            column="ENTREZID",
                            multiVals="first")
}

kk = compareCluster(genesEntrez, fun="enrichKEGG", organism= "hsa", pvalueCutoff = 1)
kkb <- kk
levels <- names(genes)[order(names(genes))]
kkb@compareClusterResult$Cluster <- factor(kkb@compareClusterResult$Cluster, levels = levels)

finalplot1 <- dotplot(kkb) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_colour_gradient2(low = "chose", mid = "your", high = "colors")

pdf('filename.pdf', width = 8, height=4.5)
print(finalplot1)
dev.off()

