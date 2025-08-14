

## 5. WGCNA - run pathway analysis on WGCNA nonlinear modules ##

library(WGCNA)
library(data.table)
library(missMethyl)


#1. Load data used in GPMethylation ===============================================================================================
load(paste0(dataPath, "FetalBrain_Betas_LOOZ_EX3_noHiLowconst_newmethod_Feb23.RData"))
betas <- betas.dasen
pheno <- read.csv(paste0(dataPath, "FetalBrain_Pheno_23pcw_EX3.csv"), row.names=1, stringsAsFactors=F)
identical(rownames(pheno),colnames(betas))


#2. Betas pre-processing ==========================================================================================================
# Remove probes annotated to sex chromosomes
probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
epicMan<-epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
betas <- betas[-which(probes$CHR %in% c('X','Y')),]


#3. Load WGCNA results ============================================================================================================
load(paste0(WGCNAPath,"Blockwise_sft12_min1000_max10000_signed_EX3_23_nonlinear_lengthscale_LLR.Rdata"))
module_df <- data.frame(probe = names(netwk$colors), colors = labels2colors(netwk$colors))


#4. Make dummy results table ======================================================================================================
modules <- c('turquoise','blue','brown','yellow','green','red','All')

for(module in modules){
	print(module)
	
	res <- data.frame(probe=rownames(betas),P=rep(1,nrow(betas)))
	rownames(res) <- rownames(betas)

	if(module=='All'){
	  print(paste('Applying pathway analysis on all',nrow(module_df),'nonlinear probes'))
	  res$P[which(rownames(res) %in% module_df$probe)] <- 0
	}else{
	  submod <- module_df[which(module_df$colors %in% module),]
	  print(paste('Applying pathway analysis on',nrow(submod),'probes in',module,'module'))
	  res$P[which(rownames(res) %in% submod$probe)] <- 0
	}

	res.sig <- res[which(res$P==0),]

	gst <- gometh(sig.cpg=res.sig$probe, all.cpg=res$probe, collection='GO', array.type="EPIC")
	gst <- gst[order(gst$P.DE),]
	
	bonf <- 0.05/nrow(gst)
	gst.sig <- gst[gst$P.DE<bonf,]
	write.csv(gst.sig, paste0(outputPath,"missMethylPathwayAnalysis_nonlinear_fetalBrain_",module,".csv"))
}