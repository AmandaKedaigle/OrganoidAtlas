library(ggplot2)
library(writexl)

bg_files = list.files(path = ".", recursive = T, full.names = T, pattern = "*_tBG.rds")
bg_timepts = substr(sapply(strsplit(bg_files, "/"), "[", 2), 1,2)
bg_pathways = substr(sapply(strsplit(bg_files, "/"), "[", 2), 7,100000)
full_files = list.files(path = "../HumanToOrg/", recursive = T, full.names = T, pattern = "tFULL.rds")
full_timepts = substr(sapply(strsplit(full_files, "/"), "[", 4), 1,2)
rm_files = list.files(path = "../HumanToOrg/", recursive = T, full.names = T, pattern = "tRM.rds")
rm_timepts = substr(sapply(strsplit(rm_files, "/"), "[", 4), 1,2)
rm_pathways = substr(sapply(strsplit(rm_files, "/"), "[", 4), 4,100000)

length = length(bg_files)  
t_whole = data.frame(Var1 = rep(NA, length),Var2 = rep(NA, length),value = rep(NA, length), variable = rep(NA, length),timept = rep(NA, length),
                     pathway = rep(NA, length),FullFreq = rep(NA, length),diff = rep(NA, length),p = rep(NA, length),padj = rep(NA, length))
wholecount = 1

for (timept in c("1m", "3m","6m")) {
  
  length = sum(bg_timepts==timept) * 48
  t = data.frame(Var1 = rep(NA, length),Var2 = rep(NA, length),value = rep(NA, length),variable = rep(NA, length),
                 timept = rep(NA, length),pathway = rep(NA, length),FullFreq = rep(NA, length),diff = rep(NA, length),p = rep(NA, length))
  tcount = 1
  
  for (pathway in levels(factor(bg_pathways))) {
  
    bgSub = bg_files[bg_timepts==timept & bg_pathways==pathway]
    if (length(bgSub) < 3000) print(paste0(timept," ",pathway," ", length(bgSub), "< 3,000!"))
    length = length(bgSub) * 48
    bgdist = data.frame(Var1 = rep(NA, length),Var2 = rep(NA, length),variable = rep(NA, length),value = rep(NA, length),
                        FullFreq = rep(NA, length),diff = rep(NA, length))
    t_full = readRDS(full_files[full_timepts==timept])
    jcount = 1
    for (j in 1:length(bgSub)) {
      tBG = readRDS(bgSub[j])
      tBG$FullFreq = t_full$value
      tBG$diff = tBG$value - tBG$FullFreq
      bgdist[jcount:(jcount+47),] =  tBG
      jcount = jcount + 48
    }
  
    rmfile = rm_files[rm_timepts==timept & rm_pathways==pathway]
    if (length(rmfile)!=1) print(paste0(timept," ",pathway," ", length(rmfile), "Not one RM'd file!"))
    t1 = readRDS(rmfile)
    t1$timept = timept
    t1$pathway = pathway
  
    t1$FullFreq = t_full$value
    t1$diff = t1$value - t1$FullFreq
  
    t1$p = NA
    for (k in 1:nrow(t1)) {
      tBGsub = bgdist[bgdist$Var1 == t1$Var1[k] & bgdist$Var2 == t1$Var2[k],]
      t1$p[k] = 1 - (sum(tBGsub$diff < t1$diff[k])/length(bgSub))
    }
  
    t[tcount:(tcount+47),] = t1
    tcount = tcount+48
  }
  t$padj = p.adjust(t$p, method="BY")
  t_whole[wholecount:(wholecount+nrow(t)-1)] = t
  wholecount = wholecount+nrow(t)
  
}
saveRDS(t_whole, "t_whole.rds")

ggplot(t_whole, aes(pathway, -1*log(padj))) +
  geom_point(aes(color= Var1)) +
  facet_grid(timept~Var2) +
  geom_hline(yintercept = -1*log(0.05), linetype="dotted") +
  theme(axis.text.x = element_text(angle=90, hjust=1)) + 
  labs(color = "Human Cell", main = "Organoid Type Assigment")
ggsave("AllPathwaySignificance.pdf", units="in", width=65, height=12, limitsize = F)

ggplot(t_whole[t_whole$Var2 %in% c("aRG", "PN") & t_whole$timept=="3m",], aes(pathway, Var1)) +
  geom_tile(aes(fill= padj)) +
  scale_fill_gradient2(high="gray90", mid="gray80", low="firebrick", midpoint=0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  facet_grid(Var2~.) +
  labs(x="", y="Human Cell Type", main="Cells changing to being classified as aRG")
ggsave("CellsChangingToaRGandPN-3m.pdf", units="in", width=10, height=10) 

pdf(paste0("glycolysis_in_6m.pdf"), width=7.3, height=3.5, useDingbats=F)
ggplot(t_whole[t_whole$pathway=="HALLMARK_GLYCOLYSIS" & t_whole$timept=="6m",], aes(Var1, Var2)) + 
  geom_point(aes(size = padj, fill = diff, color=padj<0.05), pch=21, stroke=0.7) + theme_bw() +
  xlab("Human Cell Types") + 
  ylab("Organoid Cell Types") + 
  scale_color_manual(values=c("gray50","black"))+
  scale_size_continuous(name="Adjusted p value", range = c(6,2)) +
  scale_fill_gradient2(low="cyan4", mid="white", high="brown4", name="Change in Fraction of Human Cells") + 
  theme(axis.text.x=element_text(angle=45, hjust=0.9))
dev.off()

head(t_whole)

colnames(t_whole) = c("Human Fetal Cell Type", "Organoid Label", "var", "Frequency", "Timepoint", "Pathway Genes Removed",
                      "Frequency in Full Model", "Difference", "p", "padj")
t_whole = t_whole[,-3]
write_xlsx(t_whole, "RFAllPathwayPvalues.xlsx")
