

## Plot first principal component (PC1) vs Age for each WGCNA module ##

plotAllModuleProbes <- function(betas, pheno, module, module_df, returnSubmod=FALSE, colour, toPlot=TRUE, loop=TRUE, i=NULL, axisTextSize=19, axisTitleSize=22, plotTitleSize=23, xAxisCol='Age'){

	#1. Extract DNAm data for all probes in module ===================================================================================================
	submod <- module_df[which(module_df$colors %in% module),]                                   # Subset the module_df to the module of interest
	if(returnSubmod){ return(submod) }
		
	betas_interest <- betas[which(rownames(betas) %in% submod$probe),]                          # Match the probes in submod to the DNAm matrix, betas
	betas_interest <- betas_interest[match(submod$probe, rownames(betas_interest)),]
	if(!identical(as.character(rownames(betas_interest)),as.character(submod$probe))){ stop() } # Check probes match between new betas (betas_interest) and probes in module (submod$probe)
	betas_interest <- data.frame((betas_interest)*100)
	n.probes <- nrow(betas_interest)                                                            # Make note of number of probes in module for plot below
	
	
	#2. Plot DNAm data for probes ====================================================================================================================
	# plot first principal component of the DNAm (module eigengene) vs Age
	pca.res <- prcomp(t(na.omit(betas_interest)))                                               # Run PCA on DNAm matrix of probes in module
	pc1 <- pca.res$x[,c(1,2)]
	rownames(pc1) <- gsub("X","",rownames(pc1))                                                 # Remove the 'X' it appends to start of probe names

	# Make df of first column of pca results (PC1) and Age
	df <- data.frame(value=pc1[,1], Age=pheno[match(rownames(pc1),rownames(pheno)),xAxisCol])
	
	# Correlate PC1 and DNAm values for probes in module. If mean correlation < 0, invert PC1 values so that PC1 is a proxy for DNAm.
	pc1.betas.cor <- apply(betas_interest, 1, function(x){cor(x, pc1[,1], use = 'pairwise.complete.obs')})
	if(mean(as.numeric(pc1.betas.cor))<0){
		df$value <- (-df$value)
	}

	if(toPlot){
		p <- ggplot(df, aes(x=Age, y=value))+
		geom_point()+
		geom_smooth(method='loess',se=TRUE, size=0.5, colour=colour) +
		theme_minimal()+
		ggtitle(paste0(module,'   (n=', n.probes,')'))+
		labs(y='Module eigengene', x='Age (pcw)')+
		theme(axis.text=element_text(size=axisTextSize), axis.title=element_text(size=axisTitleSize))+
		theme(plot.title = element_text(size=plotTitleSize))
	
		if(loop){
			myplots[[i]] <<- p
		}else{
			show(p)
		}
	}
}