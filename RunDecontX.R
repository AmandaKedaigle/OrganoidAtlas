library(celda)
options(stringsAsFactors=F)


RunDecontX <- function(seur, out.prefix) {
counts <- as.matrix(seur[['RNA']]@counts)

# Cell population annotation
Idents(seur) <- seur$FinalName
pop.label <- seq(1:length(levels(seur)))
names(pop.label) <- levels(seur)
cat("Cell population annotation:\n")
print(pop.label)
seur <- RenameIdents(seur, pop.label)
z <- as.numeric(as.character(Idents(seur)))

# Batch info
Idents(seur) <- as.factor(seur$org)
orgid <- seq(1:length(levels(seur)))
names(orgid) <- levels(seur)
cat("Batch information:\n")
print(orgid)
seur <- RenameIdents(seur, orgid)
batch <- as.numeric(as.character(Idents(seur)))

# construct decontX model using the above cell population and batch information
decontxModel <- decontX(counts = counts, z=z, batch=batch, maxIter=100, logfile=paste0("decontX_log_", out.prefix, ".txt"))

saveRDS(decontxModel, paste0("decontx_res_", out.prefix, ".rds"))

pdf(paste0("model_convergence", out.prefix, ".pdf"))
plot(decontxModel$resList$logLikelihood)
dev.off()
}


