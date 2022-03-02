# code by Sean Simmons

library(matrixStats)
library(RcppML)
library(Matrix)
library(Seurat)

#Modified cNMF

#A modified version of the cNMF algorithm, implemented in R

# param
# seur Seurat object
# numSamp Number of cells to select (if <1 no downsampling occurs)
# k The number of factors
# numRep The number of times to run NMF internally before making the consensus
# rho Internal parameter for cNMF
# t The distance used for filtering out outliers
# getSeruat If true returns Seurat object with cNMF added, else returns results as list

# return
# out list of results (only returned if getSeurat=F)
# seur Seurat object with cNMF added (only returned if getSeurat=T)
Repro_NMF=function(seur,numSamp=0,k=15,numRep=100,rho=.3,t=.03,getSeurat=F)
{
print("Get Data")
dat=seur@assays$RNA@data[seur@assays$RNA@var.features,]
if(numSamp>0)
{
dat=dat[,sample(1:dim(dat)[2],numSamp)]
getSeurat=F
}
mn=rowMeans(as.matrix(dat)>0)
dat=dat[mn>1/500,]

print(dim(dat))
print("Scale Data")
dat_orig=seur@assays$RNA@data[,colnames(dat)]
#dat=t(dat)
#SD=as.numeric(lapply(colnames(dat),function(x){sd(dat[,x])}))
#dat@x <- dat@x / rep.int(SD, diff(dat@p))
#dat=t(dat)



print(class(dat))
print("Run NMF multiple times")
out=lapply(1:numRep,function(x){
print(x);

print("Run NMF")
nmf_cur=RcppML::nmf(dat,k)
w=nmf_cur$w
print(dim(w))

return(w)
})

print("Combine!")
W=do.call(cbind,out)
W=data.frame(W)
for(i in 1:dim(W)[2]){W[i]=W[,i]/sqrt(sum(W[,i]^2))}
W=as.matrix(W)
print("Get values")
Dist=2-2*(t(W) %*% W)
L=rho*numRep
aveDist=apply(Dist,1,function(x){
x=x[order(x,decreasing=F)]
x=x[1:L]
return(mean(x))
})

W=W[,aveDist<t]

print("Cluster")
clusts=kmeans(t(W),k)
clusts=clusts$cluster
clusters=unique(clusts)


print("Combine Clusters")
vals=lapply(clusters,function(x){if(sum(clusts==x)==1){return(W[,clusts==x])};rowMedians(W[,clusts==x])})
W_new=do.call(cbind,vals)
W_new=data.frame(W_new)
for(i in 1:dim(W_new)[2]){W_new[i]=W_new[,i]/sum(abs(W_new[,i]))}
W_new=as.matrix(W_new)
print("Get Cell loadings")
H=project(dat,w=W_new)

print("Get Final Gene Loadings")
W_fin=project(dat_orig,h=H)

ret=list()
ret[["H"]]=t(H)
ret[["W"]]=t(W_fin)
H=t(H)
if(getSeurat)
{
print("Add to Seurat")
colnames(H)=sub("^","cNMF_",1:k)
rownames(H)=rownames(seur@meta.data)
seur[["cnmf"]] <- CreateDimReducObject(embeddings = H, key = "cNMF_", assay = DefaultAssay(seur))
return(seur)
}

return(ret)
}
