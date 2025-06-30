

## 3. WGCNA -  visualise modules ##

library(WGCNA)
library(data.table)
library(tidyverse)
library(magrittr)


#1. Load data used in GPMethylation =============================================================================================
load(paste0(dataPath, "FetalBrain_Betas_noHiLowConst.RData"))
betas <- betas.dasen
pheno <- read.csv(paste0(dataPath, "FetalBrain_Pheno_23pcw_EX3.csv"), row.names=1, stringsAsFactors=F)
identical(rownames(pheno),colnames(betas))


#2. Sample pre-processing =======================================================================================================
# calculate predicted age
source(refPath, "FetalClockFunction.R"))
pheno$PredAge <- FetalClock(betas)/7


#3. Betas pre-processing ========================================================================================================
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
betas <- betas[which(rownames(betas) %in% nonlinear),] #73035 (72310 minus X&Y)

	
#4. Load WGCNA blockwiseModules results =========================================================================================
load(paste0(WGCNAPath,"Blockwise_sft12_min1000_max10000_signed.Rdata")) #72310
module_df <- data.frame(probe = names(netwk$colors), colors = labels2colors(netwk$colors)) # match labels to module colours

# annotate nonlinearProbes.csv with module name
nonlinProbes <- read.csv(paste0(postprocessPath,"nonlinearProbes.csv"), row.names=1) #73035
module_df_match <- module_df[match(rownames(nonlinProbes),module_df$probe),]
identical(rownames(nonlinProbes)[-which(is.na(module_df_match$probe))], as.character(module_df_match$probe)[-which(is.na(module_df_match$probe))])
# for the 725 chrX probes that were not included in WGCNA, set module to "Untested"
module_df_match$module <- as.character(module_df_match$colors)
module_df_match$module[which(is.na(module_df_match$module))] <- 'Untested'
# append module to nonlinProbes
nonlinProbes <- cbind(nonlinProbes, Module=module_df_match$module)
write.csv(nonlinProbes, file=paste0(postprocessPath,"nonlinearProbes.csv"))


#5. Plot modules ================================================================================================================
library(ggplot2)
library(cowplot)
library(gridExtra)

source(paste0(scriptsPath,"4_3_WGCNA_PC1vsAge.r"))

# Using module colours on plots. Make colour exceptions for some colours that won't show up well on white background.
colourExceptions <- c('lightcyan','lightyellow','white','ivory','honeydew1','floralwhite')
colourAlternatives <- c('cyan','oldlace','lightgrey','lightgrey','honeydew2','lightgrey')


# Plot PC1 vs Age for each module
tb <- table(module_df$colors)               # Table of how many probes assigned to each module 
tb <- tb[order(tb, decreasing=TRUE)]        # Order table highest to lowest
modules <- names(tb)                        # Extract module names
modules <- modules[-which(modules=='grey')] # Exclude the grey module, as this is a non-specific module

myplots <- list()                           # List of ggplots for plotAllModuleProbes. Must be called 'myplots'.
for(m in modules){
	i <- which(modules==m)
	module <- m
	if(module %in% colourExceptions){       # Change module's plot colour if it is an exception
		indx <- which(colourExceptions %in% module)
		colour <- colourAlternatives[indx]
	}else{
		colour <- module
	}
	plotAllModuleProbes(betas=betas, pheno, module=module, module_df, colour=colour, i=i, axisTextSize=12, axisTitleSize=12, plotTitleSize=15, xAxisCol='Age')
}

pdf(paste0(plotPath,"fetalNonlinear_WGCNAmodules_PC1vsAge.pdf"), width=10, height=6)
plot_grid(plotlist=myplots,nrow=2, ncol=3)  # Plot the 9 modules in 3x3 matrix	
dev.off()