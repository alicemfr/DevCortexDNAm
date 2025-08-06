
# Contains the functions to plot DNA methylation against a continuous variable e.g. Age.
# Option to colour points (`colourBy`) and change plot character (`pchBy`) by a discrete variable e.g. Sex.
# Option to add loess fitted lines over data. Multiple lines can be fitted based on discrete variable specified in `loessBy` argument, e.g. Cell-type

# Example:

	# age.group='Fetal'
	# colourBy <- 'CellType'
	# pchBy <- NULL
	# pch <- 16
	# colours <- c(plasma(4)[2],viridis(4)[3])
	# loessBy <- 'CellType'
	# p.age <- 'fetal'
	# res.plot <- res.ord
	# res.column <- paste0(age.group,'.P')

	# if(nrow(res.plot)<50){ N <- nrow(res.plot) }else{ N <- 50 }
	
	# pdf(".pdf", height=8, width=12)
	# for(i in 1:N){
		# plotAges(res=res.plot, pheno=SampleSheet, betas=betas, plot.column='age.rescale', res.column=res.column, p.age=p.age, i=i, gene=res.plot[i,'Gene'], pchBy=pchBy, pch=pch, chr=res.plot[i,'CHR'], colourBy=colourBy, colours=colours, xaxt='n', plot.loess=TRUE, loessSpan=NULL, loessBy=loessBy, pchByLegend=TRUE, loessByLegend=FALSE, colourByLegend=TRUE)
	# }
	# dev.off()
	
	
plotLoess <- function(pheno, betas, probe, colourBy, loessSpan, plot.column, colour, nLoessPoints, Lty, lwd){
			
	for(i in 1:nlevels(factor(pheno[,paste0(colourBy)]))){
	
		if(is.null(loessSpan)){ loess.span <- NULL }else{ loess.span <- loessSpan }
	
		group <- levels(factor(pheno[,paste0(colourBy)]))[i]
		samples.col <- which(pheno[,paste0(colourBy)]==group)
		cols <- colour[i]
		Xvals <- pheno[samples.col,paste0(plot.column)]
		Yvals <- as.numeric(betas[probe,samples.col])*100
		mean.betas <- c()
		
		# calculate mean y-val at each unique x-val
			for(p in 1:length(unique(Xvals))){
				X <- unique(Xvals)[p]
				indx <- which(Xvals==X)
				if(length(indx)==1){
					mean.betas <- c(mean.betas, Yvals[indx])
				}else{
					mean.betas <- c(mean.betas, mean(Yvals[indx]))
				}
			}
			
		mean.betas <- as.numeric(mean.betas)
		perc90 <- as.numeric(sort(mean.betas)[0.90*length(mean.betas)])	#90th percentile
		perc10 <- as.numeric(sort(mean.betas)[0.10*length(mean.betas)])	#10th percentile
		mean.betas.range <- perc90 - perc10

		
		# choosing loess span
		if(is.null(loess.span)){
			if(mean.betas.range>=60){ 
				loess.span=0.4 
			}
			if((40 <= mean.betas.range) & (mean.betas.range<60)){ 
				loess.span=0.8
			}
			if(mean.betas.range<40){ 
				loess.span=1
			}
		}else{
			loess.span=loessSpan
		}
		
		plotting.mean.Yvals <- mean.betas[order(unique(Xvals))]
		plotting.mean.Xvals <- unique(Xvals)[order(unique(Xvals))] #not actually mean x-values - just called this to be synchronous with plotting.mean.Yvals
		
		df.loess <- loess(plotting.mean.Yvals ~ plotting.mean.Xvals, span=loess.span)
		x <- data.frame(Age=seq(from=(min(plotting.mean.Xvals)+1),to=(max(plotting.mean.Xvals)-3),length.out=nLoessPoints))	# exluding the extreme x-vals as these can dramatically influence the line
		y <- predict(df.loess, newdata=x$Age)
		lines(as.numeric(x$Age), as.numeric(y), col=cols, lty=Lty, lwd=lwd)

	}
}
			
			
plotAges <- function(res, pheno, betas, plot.column, res.column, i, colourBy=NULL, colourByLegend=TRUE, p.age=NULL, colours=NULL, pchBy=NULL, pchByLegend=TRUE, gene=NULL, chr=NULL, subsetTo=NULL, subsetToCol=NULL, xlab='', cexPoint=1.2, cexMain=1.8, cexLab=1.7, cexAxis=1.5, pch=16, ablineOn=TRUE, ablinePos=c(40), axisPos=c(40,60,75), txtLabels=c('Pre-natal','Post-natal'),txtXPos=c(20,78),txtYPos=c(105,105), extraDetails=TRUE, extraDetailsLabels=c('6',c('40 / 0 ','8','30'),'95'), ylim=c(0,105), extendXlim=TRUE, cexLegend=1.4, legendPos='topright', plot.loess=FALSE, loessSpan=NULL, nLoessPoints=50, lwd=2, loessBy=NULL, loessByLty=c(1,2), loessByLegend=TRUE, ...){

	par(mar=c(5.1, 5, 4.1, 2.1))
	
	probe <- rownames(res)[i]                                                   # probe to plot
	p <- signif(res[i,paste0(res.column)], digits=3)                            # extract p-value from results

	# create plot title
	if(!is.null(p.age)){ 	age <- paste0('.',p.age)   }                        # age group this p-value refers to
	ttl <- probe                                                                # title always has probe ID
	if(!is.null(gene)){ if(gene!=''){	ttl <- paste0(ttl, ' - ',gene)   } }    # append gene name if present
	if(!is.null(chr)){ 		ttl <- paste0(ttl, ' - CHR ',chr)   }               # append chr number
	ttl <- paste0(ttl,	"\np", age, " = ", p)                                   # append p-value

	# subset data if applicable
	if(!is.null(subsetTo)){
		pheno <- pheno[which(pheno[,paste0(subsetToCol)] %in% subsetTo),]
		pheno <- droplevels(pheno)
		betas <- betas[,which(colnames(betas) %in% rownames(pheno))]
	}
	
	# establish x-axis limits
	Xlim <- c(min(pheno[,paste0(plot.column)]), (max(pheno[,paste0(plot.column)]))) # min to max value of plot.column
	if(extendXlim){	Xlim[2] <- Xlim[2]+20 }                                         # extend the x-axis so legend doesn't block figure
	
	# data point effects
	if(is.null(colourBy)){								
		pheno$dummy <- rep(1, nrow(pheno)) # if not colouring points, create dummy variable column so code below still works
		colourBy <- 'dummy'
	}
	if(!is.null(colours)){ 	colour <- colours 	}else{	colour <- rainbow(nlevels(factor(pheno[,paste0(colourBy)])))	} # rainbow() unless manual colours specified
	if(!is.null(pchBy)){ Pch=pch[factor(pheno[,paste0(pchBy)])]}else{ Pch=pch }
	
	#1. Plot --------------------------------------------------------------------------------------------------------------------
	plot(pheno[,paste0(plot.column)], as.numeric(betas[probe,])*100, ylim=ylim,ylab='DNA methylation (%)',xlab=xlab,main=ttl,
	col=colour[factor(pheno[,paste0(colourBy)])], pch=Pch, cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexMain, xlim=Xlim, cex=cexPoint, ...)

	#2. Add loess lines ---------------------------------------------------------------------------------------------------------
	if(plot.loess==TRUE){
		
		# check if loessBy exists
		# check if loessBy and colourBy are the same
		if(!is.null(loessBy)){
			loess.NULL <- FALSE
			loess.col.same <- loessBy==colourBy
		}else{
			loess.NULL <- TRUE
			if(is.null(colourBy)){
				loess.col.same <- TRUE
			}else{
				loess.col.same <- FALSE
			}
		}
		
		# if loessBy exists and is not the same as colourBy...
		if(!loess.NULL & !loess.col.same){
			for(l in 1:nlevels(factor(pheno[,paste0(loessBy)]))){                   # for each level of the loessBy factor...
				group.loess <- levels(factor(pheno[,paste0(loessBy)]))[l]

				samples.col.loess <- which(pheno[,paste0(loessBy)]==group.loess)
				pheno.loess <- pheno[samples.col.loess,]                            # ...extract the samples...
				betas.loess <- betas[,samples.col.loess]                            # ...and corresponding DNAm.
			
				plotLoess(pheno=pheno.loess, betas=betas.loess, probe=probe, colourBy=colourBy, loessSpan=loessSpan, plot.column=plot.column, colour=colour, nLoessPoints=nLoessPoints, Lty=loessByLty[l], lwd=lwd)
				
			}
		}else{ # either loessBy is NULL OR loessBy and colourBy are identical
            # run DNAm for all samples together and plot one line
			
			plotLoess(pheno=pheno, betas=betas, probe=probe, colourBy=colourBy, loessSpan=loessSpan, plot.column=plot.column, colour=colour, nLoessPoints=nLoessPoints, Lty=1, lwd=lwd)
		
		}
	}

	#3. Add legend --------------------------------------------------------------------------------------------------------------
	if(is.null(pchBy)){ pchByLegend=FALSE }
	if(is.null(loessBy)){ loessByLegend=FALSE }
	
	if(colourByLegend){
		if(pchByLegend==FALSE & loessByLegend==FALSE){
			legend(legendPos, legend=levels(factor(pheno[,paste0(colourBy)])), col=colour, pch=16, cex=cexLegend)
		}
		if(pchByLegend==TRUE & loessByLegend==FALSE){
			legend(legendPos, legend=c(levels(factor(pheno[,paste0(colourBy)])),levels(factor(pheno[,paste0(pchBy)]))), 
			col=c(colour, rep('black',nlevels(factor(pheno[,paste0(pchBy)])))), 
			pch=c(rep(15, nlevels(factor(pheno[,paste0(colourBy)]))),unique(pch)), cex=cexLegend)
		}
		if(pchByLegend==FALSE & loessByLegend==TRUE){
			legend(legendPos, legend=c(levels(factor(pheno[,paste0(colourBy)])),levels(factor(pheno[,paste0(loessBy)]))), 
			col=c(colour, rep('black',nlevels(factor(pheno[,paste0(loessBy)])))), 
			pch=c(rep(15, nlevels(factor(pheno[,paste0(colourBy)]))),rep(NA, nlevels(factor(pheno[,paste0(loessBy)])))), 
			lty=c(rep(NA, nlevels(factor(pheno[,paste0(colourBy)]))),loessByLty), cex=cexLegend)
		}
		if(pchByLegend==TRUE & loessByLegend==TRUE){
			if(pchBy==loessBy){
				legend(legendPos, legend=c(levels(factor(pheno[,paste0(colourBy)])),levels(factor(pheno[,paste0(pchBy)]))), 
				col=c(colour, rep('black',nlevels(factor(pheno[,paste0(pchBy)])))), 
				pch=c(rep(15, nlevels(factor(pheno[,paste0(colourBy)]))),unique(pch)), 
				lty=c(rep(NA, nlevels(factor(pheno[,paste0(colourBy)]))),loessByLty), cex=cexLegend, merge = FALSE)
			}else{
				legend(legendPos, legend=c(levels(factor(pheno[,paste0(colourBy)])),levels(factor(pheno[,paste0(pchBy)])),levels(factor(pheno[,paste0(loessBy)]))), 
				col=c(colour, rep('black',nlevels(factor(pheno[,paste0(pchBy)]))), rep('black',nlevels(factor(pheno[,paste0(loessBy)])))), 
				pch=c(rep(15, nlevels(factor(pheno[,paste0(colourBy)]))),unique(pch), rep(NA, nlevels(factor(pheno[,paste0(loessBy)])))), 
				lty=c(rep(NA, nlevels(factor(pheno[,paste0(colourBy)]))),rep(NA, nlevels(factor(pheno[,paste0(pchBy)]))),loessByLty), cex=cexLegend)
			}
		}
	}else{
		if(pchByLegend==TRUE & loessByLegend==TRUE){
			if(pchBy==loessBy){
				legend(legendPos, legend=levels(factor(pheno[,paste0(pchBy)])), col='black', pch=unique(pch), lty=loessByLty, cex=cexLegend, merge = FALSE)
			}else{
				legend(legendPos, legend=c(levels(factor(pheno[,paste0(pchBy)])),levels(factor(pheno[,paste0(loessBy)]))), 
				col='black', pch=c(unique(pch), rep(NA,nlevels(factor(pheno[,paste0(loessBy)])))), 
				lty=c(rep(NA,nlevels(factor(pheno[,paste0(pchBy)]))),loessByLty), cex=cexLegend)
			}
		}else{
		if(pchByLegend){
			legend(legendPos, legend=levels(factor(pheno[,paste0(pchBy)])), col='black', pch=unique(pch), cex=cexLegend)
		}
		if(loessByLegend){
			legend(legendPos, legend=levels(factor(pheno[,paste0(loessBy)])), col='black', lty=loessByLty, cex=cexLegend)
		}
		}
	}
	
	#4. Finishing touches -------------------------------------------------------------------------------------------------------
	if(extraDetails!=FALSE){
		axis(side=1, at=c(0,axisPos,100), labels=extraDetailsLabels, cex.axis=cexAxis)
			if(ablineOn==TRUE){
				abline(v=ablinePos, lty=3, lwd=1, col='black')
			}
		text(x=txtXPos, y=txtYPos, labels=txtLabels, cex=cexAxis)
	}
	
	if(xlab==''){
		mtext('Age (pcw)', side=1, line=2.5, at=20, cex=cexAxis)
		mtext('Age (yrs)', side=1, line=2.5, at=70, cex=cexAxis)
	}

}