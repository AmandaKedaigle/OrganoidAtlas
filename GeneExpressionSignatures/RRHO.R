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
write.table(uu[[2]][!uu[[2]] %in% uu[[3]]], "RRHO2/OutofDiagBoth/PNvExN-upOnlyInHuman.txt", quote = F, row.names = F, col.names = F)
write.table(uu[[1]][!uu[[1]] %in% uu[[3]]], "RRHO2/OutofDiagBoth/PNvExN-upOnlyInOrg.txt", quote = F, row.names = F, col.names = F)

#Look at overlap with gene lists
d1 = "~/Documents/Arlotta/Organoids/Ana_Amanda/Paper/Signatures_071321/TolookforRRHO/Metabolic_pathways/"
d2 = "~/Documents/Arlotta/Organoids/Ana_Amanda/Paper/Signatures_071321/RRHO2_3mo_1/OutofDiagBoth/"
fl1 = list.files(path=d1, pattern=".txt")
fl2 = list.files(path=d2, pattern=".txt")
overlaps = data.frame() #save how big the overlap is for each comparison
for (newfile in fl1) {
  newfilepath = file.path(d1, newfile)
  for (upfile in fl2) {
    upfilepath = file.path(d2, upfile)
    newList = as.character(read.table(newfilepath)[,1])
    upList = as.character(read.table(upfilepath)[,1])
    overlap = intersect(newList, upList)
    lapply(overlap, write, paste0(substr(newfile, 1, nchar(newfile)-4),
                                  substr(upfile, 1, nchar(upfile)-4),
                                  "_overlap.txt"), append=TRUE)
    overlaps[newfile,upfile] = length(overlap)
  }
}
#hypergeometric test
ps = data.frame()
for (r in rownames(overlaps)) {
  for (c in colnames(overlaps)) {
    universe = 33694 #This is the number of genes in one of our seurat objects
    rowsize = length(as.character(read.table(file.path(d1, r))[,1]))
    colsize = length(as.character(read.table(file.path(d2, c))[,1]))
    numoverlap = overlaps[r,c]
    p = phyper(numoverlap-1, rowsize, universe-rowsize, colsize, lower.tail = F)
    ps[r,c] = p
  }
}
numTests = sum(1:(nrow(overlaps)-1))
padjs = apply(ps, c(1,2), function (x) min(x*numTests, 1))
padjs = apply(padjs, c(1,2), function (x) max(x, 10^-100))
logps = as.matrix(-1* log10(padjs))
library(corrplot)
tiff('Both_3mo_met_Poliu.tiff', width = 3000, height= 8000, res = 300)
plot.new()
corrplot(logps, is.corr=F, type = "lower", diag = F, tl.col = "black", cl.length=2)
dev.off()
