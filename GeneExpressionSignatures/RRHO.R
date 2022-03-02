#Compare 2 gene lists (DEGs in different conditions) using the RRHO algorithm
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943622/
#https://www.nature.com/articles/s41598-018-27903-2
#Amanda Kedaigle 11/24/2020

library(RRHO2)
library(ggplot2)
library(stats)

# Create "gene" lists:
# First column is the element (possibly gene) identifier, and the second
#is its value on which to sort. For differential gene expression, values are often
#-log10(P-value) * sign(effect)

for (of in list.files(path = "DESeq-organoids-1mo/")) {
  for (hf in list.files(path = "DESeq-Geschwind//")) {
    on = strsplit(of, ".", fixed=T)[[1]][[1]]
    hn = strsplit(hf, ".", fixed=T)[[1]][[1]]
    
    organoid = readxl::read_xlsx(paste0("DESeq-organoids-1mo/", of))
    org.list = data.frame(gene=organoid$gene, value = (-1*log10(organoid$padj)*ifelse(organoid$log2FoldChange>0,1,-1)))
    org.list = org.list[!is.na(org.list$value),]
    org.sig = organoid[organoid$padj<0.015 & organoid$log2FoldChange>1.5,]$gene

    human = readxl::read_xlsx(paste0("DESeq-Geschwind/", hf))
    human.list = data.frame(gene=human$gene, value = (-1*log10(human$padj)*ifelse(human$log2FoldChange>0,1,-1)))
    human.list = human.list[!is.na(human.list$value),]
    human.sig = human[human$padj<0.015 & human$log2FoldChange>1.5,]$gene

    org.list = org.list[org.list$gene %in% human.list$gene,]
    human.list = human.list[human.list$gene %in% org.list$gene,]

    RRHO2_obj <-  RRHO2_initialize(org.list, human.list, labels = c(paste0("Organoid ",on), paste0("Human ", hn)),  boundary=0.03)
    png(paste0("RRHO2-1m/",on, "v", hn, ".png"))
    RRHO2_heatmap(RRHO2_obj)
    dev.off()
    
    }
}


#Explore gene list differences...this is to look at the genes changed in one v the other
uu = RRHO2_obj$genelist_uu #up in both - 1 is up in org, 2 is up in human, 3 is both

write.table(uu[[3]], "RRHO2/PNvExN-upInBoth.txt", quote = F, row.names = F, col.names = F)
write.table(uu[[2]][!uu[[2]] %in% uu[[3]]], "RRHO2/PNvExN-upOnlyInHuman.txt", quote = F, row.names = F, col.names = F)
write.table(uu[[1]][!uu[[1]] %in% uu[[3]]], "RRHO2/PNvExN-upOnlyInOrg.txt", quote = F, row.names = F, col.names = F)

#GO
rrhoaRG <- as.character(read.table("RRHO2/aRGvvRG-upOnlyInOrg.txt")$V1)
rrhoPN <- as.character(read.table("RRHO2/PNvExN-upOnlyInOrg.txt")$V1)

genes = list(rrhoaRG, rrhoPN)
names(genes) = c("aRGsig-OrganoidOnly", "PNsig-OrganoidOnly")

ck = compareCluster(genes, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
ck2 = clusterProfiler::simplify(ck,cutoff=0.7, by="p.adjust", select_fun=min)
ck2b <- ck2
levels <- names(genes)[order(names(genes))]
ck2b@compareClusterResult$Cluster <- factor(ck2b@compareClusterResult$Cluster, levels = levels)

finalplot <- dotplot(ck2b, showCategory = 3) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_colour_gradient2(low = "#c2a5cf", mid = "#984ea3", high = "#ff7f00")

pdf('enrichGO_PNandaRG_cat3.pdf', width = 8.5, height=4.5)
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

finalplot1 <- dotplot(kkb, showCategory=5) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("enrichKEGG Cell Types") + scale_colour_gradient2(low = "#c2a5cf", mid = "#984ea3", high = "#ff7f00")

pdf('enrichKEGG_PNandaRG_cat5.pdf', width = 8.5, height=4.5)
print(finalplot1)
dev.off()

