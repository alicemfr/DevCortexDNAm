

## 1. WGCNA - soft thresholding ##

library(WGCNA)
library(data.table)
library(cluster)
library(data.table)


#1. Load data used in GPMethylation =============================================================================================
load(paste0(dataPath, "FetalBrain_Betas_noHiLowConst.RData"))
betas <- betas.dasen


#2. Betas pre-processing ========================================================================================================
# Remove probes annotated to sex chromosomes
probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
epicMan <- epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
betas <- betas[-which(probes$CHR %in% c('X','Y')),]

# Subset betas to nonlinear probes
gpstats <- fread(paste0(resultsPath, "gpstats.tsv"), header=T, sep='\t', data.table=F) #543796
nonlinear_metrics <- readRDS(paste0(postprocessPath, "nonlinearMetrics.rds")) #73035
gpstats.nonlinear <- gpstats[which(gpstats$Probe %in% rownames(nonlinear_metrics)),]
nonlinear <- gpstats.nonlinear$Probe
betas <- betas[which(rownames(betas) %in% nonlinear),] #73035
sigma <- apply(betas, 1, sd)
betas.var <- betas[order(sigma, decreasing = TRUE),]
tbetas <- t(betas.var)


#3. Determine soft threshold ====================================================================================================
powers <- 1:30
block <- blockSize(ncol(tbetas), rectangularBlocks = TRUE, maxMemoryAllocation = NULL, overheadFactor = 3) 
sft <- pickSoftThreshold(tbetas, powerVector = powers, verbose = 5, networkType = "unsigned", blockSize = block)
saveRDS(sft, file=paste0(WGCNAPath,"softThresh.rds"))