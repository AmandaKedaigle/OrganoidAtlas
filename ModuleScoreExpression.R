library(Seurat)
library(dplyr)
library(DescTools)
library(ggplot2)
library(multcompView)

seur = readRDS("~/Documents/wt_merge/FinalObjects/Downsample_modules/3m-noUnkCR.rds")

seur = subset(seur, subset=FinalName!="Unknown")
#seur=subset(seur, subset=FinalName!="Cajal Retzius") #In 3m also remove C-R which has 99 cells

for (column in c("HallmarkHypoxia", "HallmarkGlycolysis", "HallmarkOxphos")) {

dat = data.frame("organoid" = seur$org, CellType = seur$FinalName, val = seur@meta.data[,column])
dat$org_ct = paste(dat$organoid, dat$CellType, sep="_")

dat.avg <- dat %>%
  group_by(org_ct) %>%
  summarize(Mean = mean(val))
dat.avg$ct = gsub("^.*_", "", dat.avg$org_ct)

order = dat.avg %>% group_by(ct) %>% summarize(avg=mean(Mean))
order = order[order(order$avg,decreasing=T),]$ct
seur$FinalName = factor(seur$FinalName, levels=order)
dat.avg$ct = factor(dat.avg$ct, levels=order)

#res = aov(val ~ CellType, data = dat)
res = aov(Mean ~ ct, data=dat.avg)
s=summary(res)
#6m both glyc p<2e-16, only aRG has all of em<0.05 (<2e-16). Oxphos also <2e-16 none w all. Hypox same, aRG all < 4.4e-14
#3m hypox <2e-16, all of aRG, PN, and oRG have all p<0.001. Oxphos p =0.014 none w all. Glyc both p <2e-16, aRG and PN all <0.001
#1m glyc p = 5.11e-15 none w all, hypox = 1.34e-14 with Cortical Hem all of em <1.5e-5, oxphos = 6.54e-9 none w all
t = PostHocTest(res, method = "hsd", ordered=T)
pval = as.data.frame(t$ct)$pval
names(pval) = rownames(t$ct)
letters = multcompLetters(pval)
letters = letters$Letters[levels(seur$FinalName)]
letters = data.frame(names= names(letters), letters = letters)

#t1 = as.data.frame(t[[1]])
#t1$ct1 = sapply(strsplit(rownames(t1),"-"), "[",1)
#t1$ct2 = sapply(strsplit(rownames(t1),"-"), "[",2)

cols = c('#8dd3c7','#bebada', '#fb9a99', '#08519C', '#a6d854', '#fccde5',
         '#c6dbef', '#c7e9c0', '#bf812d', '#dfc27d',
         '#f6e8c3', '#8c510a', '#67000d','#a50f15',
         '#fff7bc','#fee391','#A1D99B', '#D9F0A3','#fcbba1','#80b1d3', '#fdb462',
         '#02818a', '#dd3497','#fa9fb5','#d9d9d9')
levels = c('aRG', 'IP', 'Newborn PN', 'Newborn DL PN', 'Cajal Retzius', 'Cortical hem', 
           'Preplate/Subplate', 'FOXG1- EMX1- neurons', 'Subcortical progenitors', 'Subcortical neurons',
           'Subcortical neuronal precursors', 'Subcortical interneurons','Neural crest','Neural placode', 
           'oRG','oRG II','oRG/Astroglia', 'Astroglia', 'PN', 'CFuPN', 'CPN',
           'Glial precursors', 'IN progenitors','Immature IN','Unknown')
cols = cols[match(levels(seur$FinalName),levels)]

p = VlnPlot(seur, column, group.by = "FinalName", pt.size = 0, cols=cols) + NoLegend() + ggtitle("") + theme(axis.text.x = element_blank(), axis.title.x=element_blank())

ymax = layer_scales(p)$y$range$range[[2]]

p = p + geom_point(data=dat.avg, aes(x=ct, y=Mean), fill="black")
#  annotate(geom="text",label=paste0("p=",signif(s[[1]][["Pr(>F)"]][[1]],digits = 2)), x=order[length(order)], y=ymax, size=5)

#for (ct in levels(factor(seur$FinalName))) {
#  t1s = t1[t1$ct1==ct | t1$ct2==ct,]
#  if (all(t1s$pval<0.001)) { p = p+ annotate(geom="text", label="*", x=ct, y=ymax, size=10) }
#}
p +geom_label(data=letters, aes(x=names, label=letters), y = ymax, fill="NA", label.size=0)
p
ggsave(paste0("1m",column,".pdf"))
}
       