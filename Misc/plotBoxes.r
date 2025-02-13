

## Boxplot script for lifecourse ##

plotBoxes <- function(res, pheno, betas, plot.column, res.column, p.age=NULL, i, colourBy, colours=NULL, pchBy=NULL, gene=NULL, chr=NULL, subsetTo=NULL, subsetToCol=NULL, cexAxis=1.7, cexMain=1.8, cexLab=1.5, pch=16, ablinePos=40,	ylim=c(0,105), extraDetails=TRUE, inter=FALSE, interCol=NULL, xlab='', pointColBy=NULL, boxNames=NULL, labsPos=NULL){
	
	#----------------------------------------------------
	# temporarily suppress warning messages about NAs
	defaultW <- getOption("warn") 
	options(warn = -1) 
	#----------------------------------------------------

	par(mar=c(5.1, 5, 4.1, 2.1))
	par(cex.axis=cexLab)
	
	'%ni%' <- Negate('%in%')
	gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n]}
	
	probe=rownames(res)[i]
	p=signif(res[i,paste0(res.column)], digits=3)
	
	# create plot title
	if(!is.null(p.age)){ 	age <- paste0('.',p.age)   }						# age group this p-value refers to
	ttl <- probe																# title always has probe ID
	if(!is.null(gene)){ if(gene!=''){	ttl <- paste0(ttl, ' - ',gene)   } }	# append gene name if present
	if(!is.null(chr)){ 		ttl <- paste0(ttl, ' - CHR ',chr)   }				# append chr number
	ttl <- paste0(ttl,	"\np", age, " = ", p)									# append p-value

	# subset data if applicable
	if(!is.null(subsetTo)){
	pheno <- pheno[which(pheno[,paste0(subsetToCol)] %in% subsetTo),]
	pheno <- droplevels(pheno)
	betas <- betas[,which(colnames(betas) %in% rownames(pheno))]
	}

	if(!is.null(colours)){ 	colour <- colours 	}else{	colour <- rainbow(nlevels(pheno[,paste0(plot.column)]))	}
	
	if(is.null(boxNames)){  boxNames <- levels(factor(pheno[,paste0(plot.column)]))  }else{  boxNames <- boxNames  }
	
	if(inter){
		tb <- table(factor(pheno[,paste0(plot.column)]),pheno[,paste0(interCol)])
		colour.seq <- seq(from=1,to=(ncol(tb)*nrow(tb)),2)
		colour2 <- rep(colour[1], (ncol(tb)*nrow(tb)))
		colour2[colour.seq] <- colour[2]
		boxplot(as.numeric(betas[probe,])*100 ~ factor(pheno[,paste0(plot.column)])*pheno[,paste0(interCol)], col=colour2, ylab='DNA methylation (%)', main=ttl, cex.lab=cexLab, cex.main=cexMain, xlab=xlab, outline=FALSE, ylim=ylim, names=boxNames)
		if(!is.null(pointColBy)){
			legend("topright", legend=levels(factor(pheno[,paste0(pointColBy)])), col=gg_color_hue(nlevels(factor(pheno[,paste0(pointColBy)]))), pch=16)
		}
		
		if(!is.null(pchBy)){
			for(i in 1:nlevels(factor(pheno[,paste0(pchBy)]))){
			group <- levels(factor(pheno[,paste0(pchBy)]))[i]
			samples.pch <- which(pheno[,paste0(pchBy)]==group)
		
						
				if(!is.null(pointColBy)){
					for(j in 1:nlevels(factor(pheno[,paste0(pointColBy)]))){
						group <- levels(factor(pheno[,paste0(pointColBy)]))[j]
						samples.col <- which(pheno[,paste0(pointColBy)]==group)
						cols <- gg_color_hue(nlevels(factor(pheno[,paste0(pointColBy)])))[j]
						samples.pch.col <- samples.pch[which(samples.pch %in% samples.col)]
						if(length(samples.pch.col)!=0){

						pheno.plot <- pheno[samples.pch.col,which(colnames(pheno) %in% c(plot.column, interCol))]
						betas.plot <- betas[probe,samples.pch.col]
						if(!all(levels(pheno[,paste0(plot.column)]) %in% unique(pheno.plot[,paste0(plot.column)]))){
							plot.column.missing <- levels(pheno[,paste0(plot.column)])[which(levels(pheno[,paste0(plot.column)]) %ni% unique(pheno.plot[,paste0(plot.column)]))]
							pheno.plot <- rbind(pheno.plot, c(plot.column.missing,levels(pheno[,paste0(plot.column)])[1]))
							betas.plot <- cbind(betas.plot, NA)
						}
						if(!all(levels(pheno[,paste0(interCol)]) %in% unique(pheno.plot[,paste0(interCol)]))){
							interCol.column.missing <- levels(pheno[,paste0(interCol)])[which(levels(pheno[,paste0(interCol)]) %ni% unique(pheno.plot[,paste0(interCol)]))]
							pheno.plot <- rbind(pheno.plot, c(interCol.column.missing,levels(pheno[,paste0(interCol)])[1]))
							betas.plot <- cbind(betas.plot, NA)
						}
						
						
						
						stripchart(as.numeric(betas.plot)*100 ~ factor(pheno.plot[,paste0(plot.column)])*pheno.plot[,paste0(interCol)], method = "jitter",pch=pch[i],add = TRUE, vertical=TRUE, col=cols)
						}
					}
				}else{
					# Take into account when the samples don't have all of the levels required for the boxplot (and assigning levels doesn't work!)
					pheno.plot <- pheno[samples.pch,which(colnames(pheno) %in% c(plot.column, interCol))]
					betas.plot <- betas[probe,samples.pch]
					if(!all(levels(pheno[,paste0(plot.column)]) %in% unique(pheno.plot[,paste0(plot.column)]))){
						plot.column.missing <- levels(pheno[,paste0(plot.column)])[which(levels(pheno[,paste0(plot.column)]) %ni% unique(pheno.plot[,paste0(plot.column)]))]
						pheno.plot <- rbind(pheno.plot, c(plot.column.missing,levels(pheno[,paste0(plot.column)])[1]))
						betas.plot <- cbind(betas.plot, NA)
					}
					if(!all(levels(pheno[,paste0(interCol)]) %in% unique(pheno.plot[,paste0(interCol)]))){
						interCol.column.missing <- levels(pheno[,paste0(interCol)])[which(levels(pheno[,paste0(interCol)]) %ni% unique(pheno.plot[,paste0(interCol)]))]
						pheno.plot <- rbind(pheno.plot, c(interCol.column.missing,levels(pheno[,paste0(interCol)])[1]))
						betas.plot <- cbind(betas.plot, NA)
					}
					
			
					stripchart(as.numeric(betas.plot)*100 ~ factor(pheno.plot[,paste0(plot.column)])*pheno.plot[,paste0(interCol)], method = "jitter",pch=pch[i],add = TRUE, vertical=TRUE) 
				}
			}

		}else{
			if(!is.null(pointColBy)){
				for(i in 1:nlevels(factor(pheno[,paste0(pointColBy)]))){
						group <- levels(factor(pheno[,paste0(pointColBy)]))[i]
						samples.col <- which(pheno[,paste0(pointColBy)]==group)
						cols <- gg_color_hue(nlevels(factor(pheno[,paste0(pointColBy)])))[i]
						stripchart(as.numeric(betas[probe,samples.col])*100 ~ factor(pheno[samples.col,paste0(plot.column)])*pheno[samples.col,paste0(interCol)], method = "jitter",pch=pch[i],add = TRUE, vertical=TRUE, col=cols)
				}
			}else{
			stripchart(as.numeric(betas[probe,])*100 ~ factor(pheno[,paste0(plot.column)])*pheno[,paste0(interCol)], method = "jitter",pch=pch,add = TRUE, vertical=TRUE) 
			}
		}
		
	}else{
		boxplot(as.numeric(betas[probe,])*100 ~ factor(pheno[,paste0(plot.column)]), col=colour, ylab='DNA methylation (%)', main=ttl, cex.lab=cexLab, cex.main=cexMain, xlab=xlab, outline=FALSE, ylim=c(min(as.numeric(betas[probe,])*100),max(as.numeric(betas[probe,])*100)), names=boxNames)
		stripchart(as.numeric(betas[probe,])*100 ~ factor(pheno[,paste0(plot.column)]), method = "jitter",pch=pch,add = TRUE, vertical=TRUE)
	}
	
	if(xlab=='' & inter==TRUE){
		labs <- levels(pheno[,paste0(interCol)])
		for(lb in 1:length(labs)){
			mtext(labs[lb], side=1, line=2.5, at=labsPos[lb], cex=cexAxis)
		}
	}
	
	#----------------------------------------------------
	# turn back on warning messages
	options(warn = defaultW)
	#----------------------------------------------------
}