## Usage: RRHO2.R <input_DARlist1> <input_DARlist2> <output_folder>
## Expects the DAR lists to be .RDS objects, typically produced by DAR_analysis_Green_Clusters.R or similar
library(RRHO2)
library(ggplot2)
library(stats)
args <- commandlinearguments()
orgName <- args[1]
humName <- args[2]
outFolder <- args[3]

dir.create(outFolder)
organoid <- readRDS(orgName)
human <- readRDS(humName)
fixRanges <- function(myRanges){
  outRanges <- c()
  for(thing in myRanges){
    tmp <- unlist(strsplit(thing, "-"))
    tmp[2] <- as.character(as.numeric(tmp[2])+1)
    outRanges <- c(outRanges, paste(tmp, collapse="-"))
  }
  return(outRanges)
}
ProperlyFormat <- function(myStr){
  myStr <- gsub(" ", "_", myStr)
  myStr <- gsub("/", "_", myStr)
  return(myStr)
}
ps = data.frame()
for (orgType in names(organoid)) {
  for (humType in names(human)) {
    orgDat <- organoid[[orgType]]
    org.list = data.frame(peak=fixRanges(orgDat$peak), value = (-1*log10(orgDat$p_val_adj)*ifelse(orgDat$avg_log2FC>0,1,-1)))
    org.list = org.list[!is.na(org.list$value),]
    
    humDat <- human[[humType]]
    human.list = data.frame(peak=humDat$peak, value = (-1*log10(humDat$p_val_adj)*ifelse(humDat$avg_log2FC>0,1,-1)))
    human.list = human.list[!is.na(human.list$value),]
    
    org.list = org.list[org.list$peak %in% human.list$peak,]
    human.list = human.list[human.list$peak %in% org.list$peak,]
    RRHO2_obj <-  RRHO2_initialize(org.list, human.list, labels = c(paste0("Organoid ",ProperlyFormat(orgType)), paste0("Human ", ProperlyFormat(humType))),  boundary=0.03)

    png(paste0(outDir, "/", ProperlyFormat(orgType), "v", ProperlyFormat(humType), ".png"))
    RRHO2_heatmap(RRHO2_obj)
    dev.off()

    row = which.min(!is.na(RRHO2_obj$hypermat[,1]))-1
    col = which.min(!is.na(RRHO2_obj$hypermat[1,]))-1
    quad = RRHO2_obj$hypermat[1:row, 1:col]
    avg = mean(quad)
    normavg = avg/mean(RRHO2_obj$hypermat, na.rm=T)
    pR = data.frame("OrgType" = ProperlyFormat(orgType), "HumanType" = ProperlyFormat(humType), "avgp" = avg, "normavg" = normavg)
    ps = rbind(ps, pR)
  }
}

#Plot maximum up-up pvalue for each comparison
psC = ps
#ps$OrgType = factor(ps$OrgType, levels = c("Subcortical_progenitors", "Subcortical_neurons", "Cortical_hem", "aRG", "IP", "Newborn_PN", "Newborn_DL_PN")) #Specific for 1m
ps$OrgType = factor(ps$OrgType, levels = c("aRG", "oRG", "IP", "PN","CPN", "CFuPN")) #Specific for the 3m sample
#ps$OrgType = factor(ps$OrgType, levels = c("aRG", "oRG", "Glial_precursors", "oRG_Astroglia","Astroglia", "IP", "PN", "CPN","IN_progenitors", "Imature_IN")) #Specific for 6m
ps$HumanType = factor(ps$HumanType, levels = c("oligo", "glia", "IPC", "exNeuron", "intNeuron"))
ggplot(ps, aes(OrgType, HumanType, fill=normavg)) +
  geom_tile(color="white", size=3) +
  scale_fill_gradient2(low="navy", mid="yellow", high="red3", midpoint=1.55) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0(outDir,"/normalizedAvg-uu-pval-allComparisons.pdf"), width=5, height=6.1)
saveRDS(ps, paste0(outDir,"/normalizedAvg-uu-pvals.rds"))
