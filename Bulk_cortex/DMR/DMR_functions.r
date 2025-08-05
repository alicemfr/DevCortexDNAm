	
	
## Functions to identify and plot DMRs from EWAS results ##

# run.dmrff() - applies dmrff() to EWAS results
# plotDMR() - calls miniman() to plot DMR results from run.dmrff()
# plotDMR_wrap() - a wrapper for plotDMR()


library(dmrff)
library(data.table)

source("DMR_miniman.r") # contains miniman function
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE)
	

run.dmrff <- function(betas, res, ES.col, SE.col, P.col, sig.thresh=9e-8, dmr.filename, stats.filename){

	# rename columns
	colnames(res)[which(colnames(res)==ES.col)] <- 'Estimate'
	colnames(res)[which(colnames(res)==SE.col)] <- 'Error'
	colnames(res)[which(colnames(res)==P.col)] <- 'P'

	# subset to significant probes
	res.sig <- res[which(res[,'P'] < sig.thresh), ]
	res.ord <- res.sig[order(res.sig[,'P']), ]
	
	# reorder betas to match top results
	methylation <- betas[match(rownames(res.ord), rownames(betas)),] 
		
	# add annotation data
	epicMan <- epicManifest[match(rownames(res.ord), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO")]
	res.ord <- cbind(res.ord, as.data.frame(epicMan))

	stats <- res.ord[,c('Estimate','Error','P','CHR','MAPINFO')]
	
	# remove probes not annotated to a gene
	identical(rownames(methylation), rownames(stats))
	if(any(is.na(stats$CHR))){
		methylation <- methylation[-which(is.na(stats$CHR)),]
		stats <- stats[-which(is.na(stats$CHR)),]
	}
	
	# sort probe position (MAPINFO) within chromosome (CHR)
	stats2 <- stats[order(stats$CHR),]
	stats.OrderedPosWithinChr <- NULL
	
	for(chr in unique(stats2$CHR)){
		stats.chr <- stats2[stats2$CHR==chr,]
		stats.chr.orderPos <- stats.chr[order(stats.chr$MAPINFO),]
		stats.OrderedPosWithinChr <- rbind(stats.OrderedPosWithinChr, stats.chr.orderPos)
	}
	stats <- stats.OrderedPosWithinChr
	saveRDS(stats, paste0(stats.filename, '.rds'))

	# run dmrff()
	print("Starting dmrff")
	dmrs <- dmrff(estimate=stats$Estimate,
				  se=stats$Error,
				  p.value=stats$P,
				  methylation=methylation,
				  chr=stats$CHR,
				  pos=stats$MAPINFO,
				  maxgap=500,
				  verbose=T)
	print("Saving results")
	saveRDS(dmrs, paste0(dmr.filename, '.rds'))
}


plotDMR <- function(res, chr, pos.range, result, P.col, ES.col, ttl, pos.shaded.regions, pad, ...){
	miniman(data = res, chr = chr, range = pos.range, result=result, Pcol=P.col, ESdat=ES.col, pad = pad, ESlevel = 0, multiply = 1, negcol = "red", cexgene = 0.6, cpgcol ="forestgreen",  genelines = 5,pch=1,col = "black", cex = 1,nullcol = "black", poscol = "blue", P.thresh=9e-8, nonsigcol="black", shadeDMR=TRUE, main=ttl, pos.shaded.regions=pos.shaded.regions, ...)
}


plotDMR_wrap <- function(dmrs, chr, pad=100000, result, P.col, ES.col, pdfName=NULL, buffer=1000, manual.buffer=NULL, ttl=NULL, pdfOn=TRUE, plotCHR=TRUE, plotLoop=FALSE, rw=NULL, pdf.width=7, pdf.height=7, ...){
	
	# pad = base pairs to plot either side of DMR on x-axis
	if(length(pad)==1) { pad <- c(pad, pad) } # if only 1 pad provided, duplicate so LHS padding is pad[1] and RHS padding is pad[2]
	
	if(plotCHR){
		dmr <- dmrs[dmrs$chr==chr,]
		dmr <- dmr[which(dmr$start>(dmr$start[1]-pad[1]) & dmr$end<(dmr$end[1]+pad[2])),]	# if plotting top DMR for a CHR, use index 1
	}		
	
	if(plotLoop){
		dmr <- dmrs[which(dmrs$start>(dmrs$start[rw]-pad[1]) & dmrs$end<(dmrs$end[rw]+pad[2])),] # if plotting many results in a loop, use the correct index of the dmrs object
	}	

	pos.shaded.regions <- as.numeric(matrix(unlist(t(dmr[,c('start','end')])),nrow=1,ncol=(nrow(dmr)*2))) # turn the start and end values into a string in correct order for plotting

	if(!is.null(buffer)){
		manual.buffer <- rep(buffer, length(pos.shaded.regions/2)) # multiply all shaded regions by same amount, n.shaded regions = length(pos.shaded.regions)/2
	}
	
	if(pdfOn){pdf(paste0(pdfName, ".pdf"), width=pdf.width, height=pdf.height)} # when plotting in a loop, set pdfOn=FALSE to create pdf of all plots tog in one file rather than indiv

	plotDMR(res=res, chr=chr, pos.range=dmr[1,c('start','end')], result=result, P.col=P.col, ES.col=ES.col, pad=pad, pos.shaded.regions=pos.shaded.regions, manual.buffer=manual.buffer, ttl=ttl, ...)
	if(pdfOn){dev.off()}
}