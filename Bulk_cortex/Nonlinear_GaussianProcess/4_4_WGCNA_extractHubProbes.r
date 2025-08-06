

## 4. WGCNA -  extract hub probes ##

library(WGCNA)
library(data.table)
library(tidyverse)
library(magrittr)
library(cowplot)


#1. Load data used in GPMethylation =============================================================================================
load(paste0(dataPath, "FetalBrain_Betas_noHiLowConst.RData"))
betas <- betas.dasen
pheno <- read.csv(paste0(dataPath, "FetalBrain_Pheno_23pcw_EX3.csv"), row.names=1, stringsAsFactors=F)
identical(rownames(pheno),colnames(betas))


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
betas <- betas[which(rownames(betas) %in% nonlinear),] #73035 (72310 minus X&Y)

	
#3. Load WGCNA blockwiseModules results =========================================================================================
load(paste0(WGCNAPath,"Blockwise_sft12_min1000_max10000_signed.Rdata")) #72310
module_df <- data.frame(probe = names(netwk$colors), colors = labels2colors(netwk$colors)) # match labels to module colours


#4. Correlate a probe's DNAm against the eigengene (PC1) for its module =========================================================
tb <- table(module_df$colors)               # Table of how many probes assigned to each module 
tb <- tb[order(tb, decreasing=TRUE)]        # Order table highest to lowest
modules <- names(tb)                        # Extract module names
modules <- modules[-which(modules=='grey')] # Exclude the grey module, as this is a non-specific module


#5. Identify probes with greatest correlation to the module eigengene ===========================================================
hub.probes <- data.frame() # df of top 20 hub probes for all modules

for(module in modules){
	print(module)
	
	submod <- module_df[which(module_df$colors %in% module),]                                   # Subset the module_df to the module of interest
	betas_interest <- betas[which(rownames(betas) %in% submod$probe),]                          # Match the probes in submod to the DNAm matrix, betas
	betas_interest <- betas_interest[match(submod$probe, rownames(betas_interest)),]
	if(!identical(as.character(rownames(betas_interest)),as.character(submod$probe))){ stop() } # Check probes match between new betas (betas_interest) and probes in module (submod$probe)
	betas_interest <- data.frame((betas_interest)*100)                                          # Multiply betas by 100 to convert to %DNAm
	
	pca.res <- prcomp(t(na.omit(betas_interest)))                                               # Run PCA on DNAm matrix of probes in module
	pc1 <- pca.res$x[,c(1,2)]
	rownames(pc1) <- gsub("X","",rownames(pc1))                                                 # Remove the 'X' it appends to start of probe names

	# Correlate PC1 and DNAm values for probes in module. If mean correlation < 0, invert PC1 values so that PC1 is a proxy for DNAm.
	pc1.betas.cor <- apply(betas_interest, 1, function(x){cor(x, pc1[,1], use = 'pairwise.complete.obs')})
	pc1.betas.cor <- abs(pc1.betas.cor)
	
	hub <- pc1.betas.cor[order(pc1.betas.cor, decreasing=T)]                                    # order probes from highest to lowest correlation
	hub <- hub[1:20]                                                                            # select top 20 hub probes
	hub <- data.frame(Probe=names(hub), Cor=as.numeric(hub), Module=rep(module, length(hub)))
	hub.probes <- rbind(hub.probes, hub)
}

saveRDS(hub.probes, file=paste0(WGCNAPath, "hubProbes.rds"))


#6. Plot hub probes for each module =============================================================================================
hub <- readRDS(paste0(WGCNAPath, "hubProbes.rds"))
epicMan <- epicManifest[match(hub$Probe, epicManifest$IlmnID),c("CHR","UCSC_RefGene_Name")]
hub <- cbind(hub, as.data.frame(epicMan))
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
hub$Gene <- unlist(lapply(hub$UCSC_RefGene_Name, uniqueAnno))
hub$UCSC_RefGene_Name <- NULL

# Using module colours on plots. Make colour exceptions for some colours that won't show up well on white background.
colourExceptions <- c('lightcyan','lightyellow','white','ivory','honeydew1','floralwhite')
colourAlternatives <- c('cyan','oldlace','lightgrey','lightgrey','honeydew2','lightgrey')


tb <- table(hub$Module)                     # Table of how many probes assigned to each module 
modules <- names(tb)                        # Extract module names
module_df <- hub
						
for(module in modules){
	print(module)

	if(module %in% colourExceptions){       # Change module's plot colour if it is an exception
		indx <- which(colourExceptions %in% module)
		colour <- colourAlternatives[indx]
	}else{
		colour <- module
	}
	
	submod <- module_df[which(module_df$Module %in% module),]            # Subset the module_df to the module of interest
	print(dim(submod))
	submod <- submod[order(submod$Cor, decreasing=T),]                   # Order probes largest to smallest correlation
	
	write.csv(submod, file=paste0(WGCNAPath, module, "_hubProbes.csv"))

	betas_interest <- betas[which(rownames(betas) %in% submod$Probe),]   # Match the probes in submod to the DNAm matrix, betas
	betas_interest <- betas_interest[match(submod$Probe, rownames(betas_interest)),]
	if(!identical(as.character(rownames(betas_interest)),as.character(submod$Probe))){ stop("rownames do not match") } # Check probes match between new betas (betas_interest) and probes in module (submod$Probe)
	betas_interest <- data.frame((betas_interest)*100)
	colnames(betas_interest) <- gsub("X","",colnames(betas_interest))
	
	if(nrow(betas_interest)>9){
		n.probes <- 9
	}else{
		n.probes <- nrow(betas_interest)
	}
	
	myplots <- list()
	
	for(i in 1:n.probes){
	
		# make title: probe ID, CHR
		ttl <- paste0(submod$Probe[i], "\n")
		if(submod$Gene[i]!=''){	ttl <- paste0(ttl, submod$Gene[i], " - ")   }
		ttl <- paste0(ttl, 'CHR',submod$CHR[i])
		
		df <- data.frame(Betas=as.numeric(betas_interest[i,]), Age=pheno[match(colnames(betas_interest),rownames(pheno)),'Age'])
			
		p <- ggplot(df, aes(x=Age, y=Betas))+														
				geom_point(size=2)+
				geom_smooth(method='loess',se=FALSE, linewidth=0.5, colour=colour) +
				theme_minimal()+
				ggtitle(ttl)+
				labs(y='DNA methylation (%)', x='Age (pcw)')+
				theme(axis.text=element_text(size=17), axis.title=element_text(size=18))+
				theme(plot.title = element_text(size=19))
		myplots[[i]] <-p
	}
	
	pdf(paste0(plotPath,module,"_hubProbes.pdf"), width=12, height=12)
	show(plot_grid(plotlist=myplots, nrow=3, ncol=3))
	dev.off()
}