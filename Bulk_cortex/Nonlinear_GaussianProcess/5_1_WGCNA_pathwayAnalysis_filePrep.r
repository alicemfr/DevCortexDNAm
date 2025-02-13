

## 5. WGCNA - Prepare nonlinear probes for pathway analysis ##

library(data.table)

uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }


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


#3. Nonlinear probes ============================================================================================================
nonlinProbes <- read.csv(paste0(postprocessPath,"nonlinearProbes.csv"), row.names=1)
nonlinear <- nonlinear[-which(nonlinear$Module=='Untested'),]
nonlinear$Module <- as.factor(nonlinear$Module)


#4. Make dummy results table ====================================================================================================
# pathway analysis relies on p-value to determine significance
# here we are using a dummy p-value column to indicate whether a probe is nonlinear or not

# res$P = all nonlinear probes vs everything else
res <- data.frame(probe=rownames(betas),P=rep(1,nrow(betas)))
rownames(res) <- rownames(betas)
res$P[rownames(res) %in% rownames(nonlinear)] <- 0 # nonlinear probes will be considered significant

# res${Module}.P e.g. res$Blue.P = each module vs everything else (including other nonlinear modules)
for(mod in levels(nonlinear$Module)){ # nonlinear probes of each module
	res[,paste0(mod,'.P')] <- rep(1, nrow(res))
	res[rownames(res) %in% rownames(nonlinear[nonlinear$Module==mod,]), paste0(mod,'.P')] <- 0
}

epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c("CHR","UCSC_RefGene_Name")]
res <- cbind(res, as.data.frame(epicMan))
res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))


#5. Save for input into pathway analysis script =================================================================================
write.csv(res, file="nonlinearProbes_forPathwayAnalysis.csv")