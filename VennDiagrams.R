#Venn Diagrams
library(VennDiagram)


# load this function: how many genes in dataset (set of lists) are overlapping between all the lists
likes <- function(dataset) {
  overlap = dataset[[1]]
  if (length(dataset) > 1) {
    for (i in 2:length(dataset)) {
      overlap = intersect(overlap,dataset[[i]])
    }
  }
  length(overlap)
}

#load this function: see how many lists we have and plot a venn diagram with that many (works for 1-4 lists)
plotVenns <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(likes(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                            likes(a[2:3]), likes(a[c(1, 3)]), likes(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                          likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                          likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2,4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), likes(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}

#make two empty lists here you will fill in the loop
allGenesUp = list()
allGenesDown = list()

#Here I am reading in lists of genes that I want to compare -- this code is for reading them from csv files, you would have to change if you have in different formats.

fnames = c("Org 1","Org 2","Org 3") #Assign names to your lists here, which will be displayed on the diagrams
  fs = c("markers_org1.csv","markers_org2.csv", "markers_org3.csv") #set to the file paths for files containing your lists of genes

for (f in 1:length(fnames)) {
  res = read.table(fs[[f]], sep = ',', header = T) #read in the file
  res = res[!is.na(res$p_val_adj),] #remove lines with p value set to NA -- set "padj" to the name of the pvalue column in your file
  geneUp = res$X[res$p_val_adj<0.01 & res$avg_logFC>0] #list of genes with significant pvals and upregulated. Set "padj" and "log2FoldChange" to column names in your file
  geneDown = res$X[res$p_val_adj<0.01 & res$avg_logFC<0] #same for downregulation
  allGenesUp[[fnames[[f]]]] = as.character(geneUp) #populating those empty lists we made earlier
  allGenesDown[[fnames[[f]]]] = as.character(geneDown)
}

#Makes the plots!

tiff('FILENAME.tiff', res = 300,  width = 800, height= 800)
venn.plot = plotVenns(allGenesUp, category = fnames, lty = "blank", fill = c('color1','color2','color3')[1:length(allGenesUp)])
dev.off()

tiff('FILENAME.tiff', res = 300, width = 800, height= 800)
venn.plot = plotVenns(allGenesDown, category = fnames, lty = "blank",fill = c('color1','color2','color3')[1:length(allGenesDown)])
dev.off()

##calculate genes that overlap up
overlap_up = calculate.overlap(allGenesUp)
overlap_up_df = data.frame(gene = as.character(), category =  as.character())
for(i in 1:length(overlap_up)){
  temp = overlap_up[[i]]
  name = rep(names(overlap_up)[[i]], length(temp)) ##name of the list, created as a vector and the length the same as the items within this element
  df = data.frame(gene = temp, category = name)
  overlap_up_df = rbind( overlap_up_df, df)
}
write.csv(overlap_up_df, "overlap_up_df_orgs.csv")

##calculate genes that overlap down 

overlap_down = calculate.overlap(allGenesDown)
overlap_down_df = data.frame(gene = as.character(), category =  as.character())
for(i in 1:length(overlap_down)){
  temp = overlap_down[[i]]
  name = rep(names(overlap_down)[[i]], length(temp)) ##name of the list, created as a vector and the length the same as the items within this element
  df = data.frame(gene = temp, category = name)
  overlap_down_df = rbind( overlap_down_df, df)
}
write.csv(overlap_down_df, "overlap_down_df_orgs.csv")


