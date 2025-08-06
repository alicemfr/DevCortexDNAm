

## Functions to calculate PCA and plot ##

#1. screePlot() = calculates variance explained by principle components of pca results. Option to plot barplot.
#2. pca.Gene() = runs PCA on betas. [...] = additional arguments such as scale=TRUE to be supplied to prcomp.
#3. plotPCA() = scatter plots of 2 PCs. Uses screePlot() to display % variance explained on axes. Optional use of PCAproject() to display projected samples on pca plot.
#4. PCAproject() = using the loadings of one set of pca results, project the loadings for other samples onto this space.
#5. pcaCorPlot() = heatmap showing significance of correlations between PCs and variables of interest.
#6. pcaCorPlot.withCov() = heatmap showing significance of linear regression coefficients between PCs and variables of interest. Ability to control for variables.


gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n]}
n=2; cols = gg_color_hue(n)
	
screePlot <- function(pca.res, nPCs, pcA, pcB, toPlot=TRUE, return.pcVars=FALSE){
  names(pca.res$sdev)=paste('PC',1:nPCs,sep='')
  pcAvar <- signif(((pca.res$sdev[pcA]^2)/sum(pca.res$sdev^2))*100, digits=3)
  pcBvar <- signif(((pca.res$sdev[pcB]^2)/sum(pca.res$sdev^2))*100, digits=3)
  if(return.pcVars){return(list(pcAvar, pcBvar))}
  if(toPlot){
	barplot((pca.res$sdev[1:nPCs])^2,beside=TRUE,las=2,names.arg=names(pca.res$sdev)[1:nPCs], col=cols[2])
	legend('topright',paste('PC',pcA,'=',pcAvar,'%\nPC',pcB,'=',pcBvar,'%'),cex=1.5, bty='n', title='Variance explained')
  }
}


pca.Gene <- function(QCmetrics.sub, betas, varProbes=FALSE, nVarprobes=NULL, probeList=NULL, ...){
  # epicManifest must be loaded ahead of time
  # function does not check for NAs; na.omit() must be run ahead of time on the betas
  # rownames(pheno) must match colnames(betas)
  
  if(rownames(pheno)!=colnames(betas)){
    stop("rownames(pheno) must match colnames(betas)")
  }
  
  betas.pca <- betas
  betas.sub <- betas.pca[,which(colnames(betas.pca) %in% rownames(QCmetrics.sub))]
  betas.sub <- betas.sub[,match(rownames(QCmetrics.sub), colnames(betas.sub))]
  identical(rownames(QCmetrics.sub), colnames(betas.sub))
  
  if(!is.null(probeList)){
    betas.sub <- betas.sub[which(rownames(betas.sub) %in% probeList),]
  }
  
  if(varProbes){
    sigma<-apply(betas.sub, 1, sd)
    betas.sub <- betas.sub[order(sigma, decreasing = TRUE)[1:nVarprobes],]
  }
  
  probes <- data.frame(Probe=rownames(betas.sub))
  rownames(probes) <- rownames(betas.sub)

  tbetas.sub <- t(betas.sub)
  pca.sub <- prcomp(tbetas.sub, ...)
  return(PCA.sub=pca.sub)

}


plotPCA <- function(pca.res, pheno, pcA=NULL, pcB=NULL, colourBy=NULL, pointBy=NULL, pchs=NULL, legend=TRUE, legendPos=NULL, cex.legend=1, title=NULL, xlim=NULL, ylim=NULL, projection=FALSE, projectData, projectPheno, colourBy.project, returnPCA=FALSE, colours=NULL, legendOrder=NULL, cex.point.condition.sizes=NULL, cex.point=1, cex.lab=1, cex.main=1, cex.axis=1, xlab.perc=NULL, ylab.perc=NULL, colourBy.legendSigFig=NULL, pointBy.legendSigFig=NULL, nPCs=10){
  
  pcVars <- screePlot(pca.res, nPCs, pcA, pcB, toPlot=FALSE, return.pcVars=TRUE)
  
  if(is.null(colourBy)){colourBy <- 'Cell_Type'}
  if(is.null(legendPos)){legendPos <- "topleft"}
  if(is.null(pcA)){pcA <- 1}
  if(is.null(pcB)){pcB <- 2}
  if(is.null(xlim)){xlim=c(min(pca.res$x[,pcA]),max(pca.res$x[,pcA]))}
  if(is.null(ylim)){ylim=c(min(pca.res$x[,pcB]),max(pca.res$x[,pcB]))}
  if(is.null(pchs)){pchs=c(16,3,2,10,5,8,4,0,18,20)}
  
  if(projection){
    print('Checking project data')
    projectData <- projectData[,c(pcA,pcB)] # extract the same PCs
    projectPheno <- projectPheno[which(rownames(projectPheno) %in% rownames(projectData)),] # ensure correct sample info for projectData
    
    pca.plot <- pca.res$x[,c(pcA,pcB)] # orginal PCA data
    print('Adding projection data to original data')
    pca.plot <- rbind(pca.plot, projectData) # add the projection data
    
    # combine the original pheno with projection pheno. If pheno doesn't already have column called colourBy.project, add a dummy column by that name.
    if(!any(colnames(pheno)==colourBy.project)){
      pheno[,paste0(colourBy.project)] <- rep('NA',nrow(pheno))
    }

    print('Creating the plotfactor column')
    pheno$colourBy <- as.character(pheno[,paste0(colourBy)])
    projectPheno$colourBy <- as.character(projectPheno[,paste0(colourBy.project)])
    pca.pheno <- data.frame(colourBy=c(pheno$colourBy, projectPheno$colourBy))
    
	if(grepl('Age',colourBy)){
		pca.pheno[,'colourBy'] <- as.numeric(as.character(pca.pheno[,'colourBy']))
		plotfactor<-factor(pca.pheno[,'colourBy'], levels = unique(pca.pheno[,'colourBy'])[order(unique(pca.pheno[,'colourBy']), decreasing=F)])
	}else{
		plotfactor<-factor(pca.pheno$colourBy, levels = c(as.character(unique(pca.pheno$colourBy))))
	}
    
	if(is.null(pointBy)){ points <- 16 }else{
		pheno$pointBy <- as.character(pheno[,paste0(pointBy)])
		projectPheno$pointBy <- as.character(projectPheno[,paste0(pointBy)])
		pca.pheno$pointBy <- c(pheno$pointBy, projectPheno$pointBy)
		pca.pheno$pointBy <- as.factor(pca.pheno$pointBy)
		points <- pchs[1:nlevels(pca.pheno$pointBy)][pca.pheno$pointBy]
	}
	pheno <- pca.pheno
  }else{
  
	if(is.null(pointBy)){ points <- 16 }else{
		pheno$pointBy <- factor(pheno[,paste0(pointBy)])
		points <- pchs[1:nlevels(pheno$pointBy)][pheno$pointBy]
	}
    pca.plot <- pca.res$x[,c(pcA,pcB)]
    
	if(grepl('Age',colourBy)){
		plotfactor<-factor(pheno[,paste0(colourBy)], levels = unique(pheno[,paste0(colourBy)])[order(unique(pheno[,paste0(colourBy)]), decreasing=F)])
	}else{
		plotfactor<-factor(pheno[,paste0(colourBy)], levels = c(as.character(unique(pheno[,paste0(colourBy)]))))
	}
	
  }
  
  if(!is.null(colours)){ colour <- colours }else{ colour <- rainbow(nlevels(plotfactor)) }
  
  
  # Plot PCA #
  if(is.null(xlab.perc)){
		plot(pca.plot[,1], pca.plot[,2], col=colour[factor(plotfactor)], pch=points, main=title, xlab=paste0('PC',pcA,' - ',pcVars[[1]],'%'), ylab=paste0('PC',pcB,' - ',pcVars[[2]],'%'), xlim=xlim, ylim=ylim, cex=cex.point.condition.sizes[pheno$pointBy], cex.lab=cex.lab, cex.main=cex.main, cex.axis=cex.axis)

  }else{ # if wanting to display own percentage values (or equivalent) on axes
	plot(pca.plot[,1], pca.plot[,2], col=colour[factor(plotfactor)], pch=points, main=title, xlab=paste0('PC',pcA,' - ',xlab.perc), ylab=paste0('PC',pcB,' - ',ylab.perc), xlim=xlim, ylim=ylim, cex=cex.point, cex.lab=cex.lab, cex.main=cex.main, cex.axis=cex.axis)
  }
  
  # Add legend #
  if(legend){
		
	if(is.null(legendOrder)){
		#legendSigFig = if displaying long numeric values on legend, option to round to significant number of figures.
		ifelse(!is.null(colourBy.legendSigFig), colourBy.levels <- signif(as.numeric(levels(factor(plotfactor))),digits=colourBy.legendSigFig), colourBy.levels <- levels(factor(plotfactor)))
		ifelse(!is.null(pointBy.legendSigFig), pointBy.levels <- signif(as.numeric(levels(factor(pheno$pointBy))),digits=pointBy.legendSigFig), pointBy.levels <- levels(factor(pheno$pointBy)))
		legendNames <- c(colourBy.levels, pointBy.levels)
		colourOrder <- c(colour, rep('black',nlevels(factor(pheno$pointBy))))
		if(!is.null(pointBy)){ 	pchOrder <- c(rep(16,nlevels(plotfactor)), pchs[1:nlevels(pheno$pointBy)])  }
	}else{
		# re-order legend labels, pch and colours to match order requested in legendOrder
		legendNames <- c(levels(factor(plotfactor)), levels(factor(pheno$pointBy)))	
		colourOrder <- c(colour, rep('black',nlevels(factor(pheno$pointBy))))[match(legendOrder, legendNames)]
		if(!is.null(pointBy)){ 	pchOrder <- c(rep(16,nlevels(plotfactor)), pchs[1:nlevels(pheno$pointBy)])[match(legendOrder, legendNames)]  }
		legendNames <- legendOrder
	}
	
	if(!is.null(pointBy)){  
		legend(paste(legendPos), legend=legendNames, col=colourOrder, pch=pchOrder, cex=cex.legend)  
	}else{
		legend(paste(legendPos), legend=legendNames, col=colourOrder, pch=points, cex=cex.legend) 
	}
  }
  
  if(returnPCA){
    return(pca.plot)
  }
}


PCAproject <- function(pca.res, testdata){ 
  # pca.res = PCA(dataset1)
  # testdata = betas of dataset2. probes in testdata need to be identical to the probes used in pca.res 
  # Note: the predict() line below is the same as doing scale(testdata, pca.res$center, pca.res$scale) %*% pca.res$rotation
  testdata <- t(testdata) # transpose betas of dataset2 as this is what was done for PCA of dataset1 in pca.Gene()
  if(identical(colnames(testdata), names(pca.res$center))){
	project.new = predict(pca.res, newdata=testdata)  
	return(project.new)
  }else{
    print("Probes in testdata need to be identical to the probes used in pca.res")
  }
}


library(corrplot)
pcaCorPlot <- function(samplesheet, columns, pca.res, nPCs=10, colLabels=NULL, returnP=FALSE){  # credit Dr Josh Harvey
	# columns = character string of column names of samplesheet containing the variables of interest
	# pca.res = pca results from prcomp() run on same samples as in samplesheet.
	# nPCs = max number of principle components to plot out
	# colLabels = optional: specify character string for renaming headers of heatmap if different to `columns`
	# returnP = logical whether to return p-values of correlation statistics
	PCAres <- pca.res
	PCs <- PCAres$x[,1:nPCs]
	testCor <- cor.mtest(data.frame(samplesheet[,columns],PCs),conf.level=0.95)
	testCorp <- testCor$p 
	testCorp <- testCorp[1:length(columns),(length(columns)+1):ncol(testCorp)]
	corPCA <- cor(PCs,samplesheet[,columns])
	if(!is.null(colLabels)){
		colnames(corPCA) <- colLabels
	}

	corrplot(corPCA,p.mat = t(testCorp), method = 'color', insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5, pch.col = 'grey20', tl.col='black', tl.cex=1.5)

	if(returnP){
		testCorp <- t(testCorp)
		colnames(testCorp) <- colnames(corPCA)
		rownames(testCorp) <- paste0('PC',1:nPCs)
		return(testCorp)
	}
}


pcaCorPlot.withCov <- function(samplesheet, columns, pca.res, nPCs=10, colLabels=NULL, covariate=NULL, returnP=FALSE){
	# columns = character string of column names of samplesheet containing the variables of interest
	# pca.res = pca results from prcomp() run on same samples as in samplesheet.
	# nPCs = max number of principle components to plot out
	# colLabels = optional: specify character string for renaming headers of heatmap if different to `columns`
	# covariate = covariate to control for. If there is a clash and the same variable is present in `columns` and `covariate`,
	    # it will become NA in `columns` and will continue to be controlled for.
	# returnP = logical whether to return p-values of correlation statistics
	PCAres <- pca.res
	PCs <- PCAres$x[,1:nPCs]
	pvals <- matrix(NA, nrow=nPCs, ncol=length(columns))
	rownames(pvals) <- paste0('PC',1:nPCs)
	colnames(pvals) <- columns
	
	if(!is.null(covariate)){
		if(covariate %in% columns){
			print("Chosen covariate is in your columns. Returning NAs for this covariate.")
		}
	}
	
	for(i in 1:nPCs){
		for(j in 1:length(columns)){
			
			if(is.null(covariate)){
				mod <- lm(PCs[,i]~samplesheet[,columns[j]])
				p <- coef(summary(mod))[2,'Pr(>|t|)']
			}else{
				if(covariate==columns[j]){
					p <- NA
				}else{
					Y <- PCs[,i]
					X <- samplesheet[,columns[j]]
					C <- samplesheet[,covariate]
					mod <- lm(reformulate(paste("X", "C", sep='+'), "Y"))
					p <- coef(summary(mod))['X','Pr(>|t|)']
				}
			}
			pvals[i,j] <- p
		}
	}
	corPCA <- cor(PCs,samplesheet[,columns])
	if(!is.null(colLabels)){
	colnames(corPCA) <- colLabels
	}
	corrplot(corPCA,p.mat = pvals, method = 'color', insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5, pch.col = 'grey20', tl.col='black', tl.cex=1.5)

	if(returnP){
		return(pvals)
	}	
}