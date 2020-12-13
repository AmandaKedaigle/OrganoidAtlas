#UpSet plots for overlaps of gene lists from Venn Diagram code

setwd("~/directory/")

library(UpSetR)

#reading in files
#in overlap files, skipping first column which just lists numbers
dataset_names = c("file1", "file2", "file3", "file4", "file5", "file6")
dataset_files = c("name1", "name2", "name3", "name4", "name5", "name6")

#Just taking the genes that overlap in all datasets
#for a 3-way overlap:
cat = "a5"
#for a 4-way overlap:
#cat = "a6"

genesF = lapply(dataset_files, read.csv, row.names=1)
genesF = lapply(genesF, function(x) { as.character(x[x$category==cat, "gene"]) })
names(genesF) = dataset_names

#Makes the plots!
p0 <- upset(fromList(genesF), nsets = 7, order.by = "degree",  mainbar.y.label = "Upregulated Gene Intersections", sets.x.label = "Genes Per List",  point.size = 2.2, nintersects = 7,  main.bar.color = "#737373", matrix.color = "#737373", sets.bar.color = c("chose", "your", "colors"))  

pdf("filename.pdf",  width = 6, height = 5)
print(p0)
dev.off()

#There are lots of formatting options, see https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
