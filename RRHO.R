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

ps = data.frame()
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
    
    #find the border of lower left quadrant
    row = which.min(!is.na(RRHO2_obj$hypermat[,1]))-1
    col = which.min(!is.na(RRHO2_obj$hypermat[1,]))-1
    quad = RRHO2_obj$hypermat[1:row, 1:col]
    avg = mean(quad)
    normavg = avg/mean(RRHO2_obj$hypermat, na.rm=T)
    pR = data.frame("OrgType" = on, "HumanType" = hn, "avgp" = avg, "normavg" = normavg)
    ps = rbind(ps, pR)
    }
}

#Plot normalized maximum up-up pvalue for each comparison
psC = ps
ps$OrgType = factor(ps$OrgType, levels = c("Subcortical progenitors", "Subcortical neurons", "Cortical hem", "aRG", "oRG", "IP", "Newborn PN", "Newborn DL PN", "Cajal Retzius"))
ps$HumanType = factor(ps$HumanType, levels = c("vRG", "oRG", "IP","ExN", "ExM", "ExM-U", "ExDp", "InCGE", "InMGE", "OPC"))
ggplot(ps, aes(OrgType, HumanType, fill=normavg)) +
  geom_tile(color="white", size=3) +
  scale_fill_gradient2(low="navy", mid="yellow", high="red3", midpoint=1.55) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave("RRHO2-1m/normalizedAvg-uu-pval-allComparisons.pdf", width=5, height=6.1)

saveRDS(ps, "RRHO2-1m/normalizedAvg-uu-pvals.rds")