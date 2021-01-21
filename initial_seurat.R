#Amanda Kedaigle 7/12/18
#generally copied from Sean's dir10X function in load_Seurat.R, with some input to make extra graphs from the Seurat tutorial

dir = '../round1_data/'
outdir = '../round1_results/'

library(Seurat)
library(stringr)
library(Matrix)

#loads all 10X lanes from a given directort
dir10X<-function(dir="",outdir="",dat=NULL,lst=c(),makeSeurat=T,minGenes=500,num=c())
{
	if(length(lst)==0) {	
		print(paste("ls ",dir,"/{7,8,9}/*/outs/filt* | grep : | sed 's/://g'",sep=""))
		lst=system(paste("ls ",dir,"/{7,8,9}/*/outs/filt* | grep : | sed 's/://g'",sep=""),intern=T)
	}	
	if(length(num)>0){
		lst=lst[num]
	}
	print(lst)
	if(is.null(dat)){
		print("Read in!")
		dat=Read10X(lst)
		print(head(dat))
		print("Fix colnames!")
		cols=colnames(dat)
		cols_new=c()
		for(col in cols){
			start=str_sub(col,1,1)
			cur=col
			if(start %in% c("A","T","G","C")){cur=paste("1_",cur,sep="")}
			cols_new<-c(cols_new,cur)
		}

		colnames(dat)=cols_new
		#print("Return!")
		saveRDS(dat, file = paste0(outdir,"dat.rds"))
	}
	print(paste("Dims: ",toString(dim(dat))))
	if(!makeSeurat){return(dat)}

	print("Make object!")
	#Seuratv3
	seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)


	#QC from tutorial
	# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = seur), value = TRUE)
	percent.mito <- Matrix::colSums(GetAssayData(seur, slot='counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
	seur[['percent.mito']] <- percent.mito
	ribo.genes <- grep("^RP[S,L]",rownames(seur), value = TRUE)
	percent.ribo <- Matrix::colSums(GetAssayData(seur, slot='counts')[ribo.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
	seur[['percent.ribo']] <- percent.ribo
	pdf(paste0(outdir,'QC_Rplots.pdf'))
	par(mfrow = c(1, 3))
	print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.mito"))
	print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.ribo"))
	print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
	dev.off()

	print("Normalize data!")
	seur = NormalizeData(seur, normalization.method = "LogNormalize", scale.factor=1000000)

	print("Get variable genes!")
	seur<-FindVariableFeatures(seur,selection.method='mean.var.plot')

	print("Regress out!")
	seur = CellCycleScoring(seur, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
	seur$CC.Difference <- seur$S.Score - seur$G2M.Score
	seur<-ScaleData(seur,features=VariableFeatures(seur),vars.to.regress=c("nCount_RNA","CC.Difference"))

	print("Run PCA!")
	seur<-RunPCA(seur,features = VariableFeatures(seur),verbose=F)

	print("Save object!")
	saveRDS(seur, file = paste0(outdir,"initial_seur.rds"))

	pdf(paste0(outdir,'pca_Rplots.pdf'))
	print(VizDimLoadings(object = seur, dims = 1:2))
	print(DimPlot(object = seur))
	DimHeatmap(object = seur, dims = 1:9, cells = 500, balanced = TRUE)
	print(ElbowPlot(object = seur, 30))
	dev.off()
}

dir10X(dir=dir, outdir=outdir)
