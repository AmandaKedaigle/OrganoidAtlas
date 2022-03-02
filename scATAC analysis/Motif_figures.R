library(ggplot2)
library(XML)
library(dplyr)
library(writexl)

#Gather Motifs together ---
#get html files for top 5 motifs per cell type
file.names = list.files(path = ".", recursive = T, full.names = T, pattern = "motif[1-5].info.html")
motif.files = list.files(path = ".", recursive = T, full.names = T, pattern = "motif[1-5].motif")

celltypes = sapply(strsplit(file.names, "peaks."), "[", 2)
celltypes = sapply(strsplit(celltypes, "/"), "[", 1)
timepts = sapply(strsplit(file.names, "/peaks."), "[", 1)
timepts = substr(timepts, nchar(timepts)-1, nchar(timepts))

tbl = data.frame("TimePoint" = character(0), "CellType" = character(0), "Motif" = character(0), 
                 "pval" = numeric(0), "foldChange" = numeric(0), "Matching TFs" = character(0), stringsAsFactors = F)
tblLong = data.frame("TimePoint" = character(0), "CellType" = character(0), "Motif" = character(0), 
                      "pval" = numeric(0), "foldChange" = numeric(0), "Matching TFs" = character(0), stringsAsFactors = F)

for (i in 1: length(file.names)) {
  f = file.names[[i]]
  h = readHTMLTable(f, header=F)
  h = h[c(1, seq(4, length(h),2))] #first table is motif stats, even numbered after 4 are best known motif matches
  h = lapply(h, transform, V2=as.character(V2))
  fn = readLines(f)[startsWith(readLines(f), "<H4")] #get motif names from titles above tables
  names(h)[2:length(h)] = sapply(strsplit(substr(fn,5,nchar(fn)), "/", fixed=T), "[", 1) #pull out just the name
  p = as.numeric(h[[1]][1,2]) #seq pvalue
  if (p <= 1e-10) {
    fc = as.numeric(substr(h[[1]][5,2],1,nchar(h[[1]][5,2])-1))/as.numeric(substr(h[[1]][7,2],1,nchar(h[[1]][7,2])-1)) #percent in targets over percent in bg
    motif = readLines(motif.files[[i]])[[1]]
    motif = strsplit(substr(motif,2,nchar(motif)), "\t")[[1]][[1]]
    tfs = list()
  
    # get motifs with match score above treshold
    opts = h[2:length(h)]
    tfs = names(opts)[sapply(opts, function(x) { as.numeric(x[2,2]) >= 0.59 })]
    celltype = celltypes[[i]]
    timept = timepts[[i]]

    tbl[nrow(tbl)+1,] = c(timept, celltype, motif, p, fc, paste(tfs, collapse=", "))
    tblLong = rbind(tblLong, data.frame("TimePoint" = timept, "CellType" = celltype, "Motif" = motif, 
                                         "pval" = p, "foldChange" = fc, "Matching TFs" = tfs, stringsAsFactors = F))
  }
}

#add RNA information
seur1 = readRDS("mito210c1_d28_111520_final.rds")
Idents(seur1) = "FinalName"
a1 = AverageExpression(seur1)
colnames(a1$RNA) = c("NewbornPN", "aRG", "CorticalHem", "NewbornDLPN", "Unknown","SubcorticalNeurons",
                     "IP", "CR", "SubcorticalProgenitors") #to match tblLong CellTypes
seur3 = readRDS("Mito210c1_3mo_102820.rds")
Idents(seur3) = "FinalName"
a3 = AverageExpression(seur3)
seur6 = readRDS("Mito210c1_6mo_102820.rds")
Idents(seur6) = "FinalName"
a6 = AverageExpression(seur6)
colnames(a6$RNA) = c("aRG","IP","oRG","oRG II","oRG-Astroglia","Astroglia","PN","CPN",
                     "GlialPrecursors","INProgenitors","ImmatureIN") #to match tblLong CellTypes

tblLong$TF_simpleName = sapply(strsplit(tblLong$Matching.TFs, "(", fixed=T), "[", 1)
PBs = lengths(strsplit(tblLong$Matching.TFs, "_", fixed=T))>1
tblLong$TF_simpleName[PBs] = sapply(strsplit(tblLong$Matching.TFs[PBs], "_", fixed=T), "[", 2)
tblLong$TF_simpleName = toupper(tblLong$TF_simpleName)

tblLong$RNA = NA
tblLong[tblLong$TimePoint=="1m","RNA"] = apply(tblLong[tblLong$TimePoint=="1m",], 1, function(x) {if (x[["TF_simpleName"]] %in% rownames(a1$RNA)) a1$RNA[x[["TF_simpleName"]], x[["CellType"]]] else NA})
tblLong[tblLong$TimePoint=="3m","RNA"] = apply(tblLong[tblLong$TimePoint=="3m",], 1, function(x) {if (x[["TF_simpleName"]] %in% rownames(a3$RNA)) a3$RNA[x[["TF_simpleName"]], x[["CellType"]]] else NA})
tblLong[tblLong$TimePoint=="6m","RNA"] = apply(tblLong[tblLong$TimePoint=="6m",], 1, function(x) {if (x[["TF_simpleName"]] %in% rownames(a6$RNA)) a6$RNA[x[["TF_simpleName"]], x[["CellType"]]] else NA})

#Write to supp file
tblShort = tblLong %>% group_by(TimePoint,CellType,Motif) %>% 
  mutate(Matching.TFs = list(Matching.TFs), ExpressedTFs = list(TF_simpleName[RNA>0.1]), RNA = NULL, TF_simpleName=NULL) %>%
  distinct()

#adapted from Cybernetic on StackOverflow
tibble_with_lists_to_xslx <- function(tibble_object, file_path_name) {
  set_lists_to_chars <- function(x) { 
    if(class(x) == 'list') { y <- paste(unlist(x[1]), sep='', collapse=', ') } else { y <- x  } 
    return(y) }
  new_frame <- data.frame(lapply(tibble_object, set_lists_to_chars), stringsAsFactors = F)
  write_xlsx(new_frame, file_path_name)
}
tibble_with_lists_to_xslx(tblShort, "compareMotifs-homerResults.xlsx")


tblSub = tblLong[tblLong$Matching.TFs %in% c("EMX1","EMX2", "LHX9(Homeobox)", "Lhx2(Homeobox)",
                                             "MEIS2","NeuroD1(bHLH)","NFIX", "NFIA", "NFIC","Sox9(HMG)",
                                             "DLX5(Homeobox)", "DLX1(Homeobox)", "NeuroG2(bHLH)", "PH0095.1_Lhx5"),]
tblSub$TF_simpleName = factor(tblSub$TF_simpleName, levels = c("EMX1", "EMX2","EN1", "LHX9","LHX2","LHX5",
                                                             "NEUROD1","NEUROG2", "NFIC","NFIA","NFIX",
                                                             "MEIS2","SOX9","DLX1","DLX5"))

#plot with RNA levels
ggplot(tblSub, aes(CellType, TF_simpleName, color=log(RNA,10), size=foldChange)) +
  geom_point() +
  scale_color_gradient2(low="white", mid="skyblue", high="royalblue") +
  facet_grid(.~TimePoint, scales="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("homerMotifs-selected-RNA.pdf", height = 7, width=9, units = "in", useDingbats = F)

