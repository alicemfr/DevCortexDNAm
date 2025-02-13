
## Created by Dr Rebecca Smith ##
  # Edited by Alice Franklin #

miniman <- function (data = NULL, chr = NULL, cpg = NULL, range = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, pch = 1, main = NULL, col = "black", cex = 1, result = NULL, pad = 30000, multiply = NULL, nullcol = "black", negcol = "black", poscol = "black", ESlevel = 0, ESdat = NULL, cexgene = 0.5, geneplot = TRUE, genelines = NULL, chrcol = NULL, mapcol = NULL, cpgcol = "forestgreen", cex.axis = 1, cex.lab = 1,P.thresh=9e-8, nonsigcol="black",colour.below.Pthresh=TRUE,shadeDMR=FALSE,pos.shaded.regions=NULL,shaded.buffer=TRUE,manual.buffer=NULL,buffer.perc=10, verbose=FALSE, ESplot = TRUE, Pcol, sigProbesOnly=FALSE, linePlot=FALSE) {

# data should just be a dataframe
# chrcol and mapcol will automatically use 'CHR' and 'MAPINFO' but can be 
# replaced. chrcol column should just be numeric
# Can either provide a chr and range or cpg (must be rowname and have 
# mapinfo and chr in the file)
# result is the column whose p-values you want to plot
# pad is amount of padding to add ot the coordinates given
# cexgene is text size for labelling genes
# geneplot will be added by default but can be stopped if = "FALSE"
# genelines will autimatically be calculated by number of transcripts but 
# can be adjusted. If more transtcripts than genelines, several will be plotted 
# on the same line
# ESdat is the effect size column to colour points by
# ESlevel is target effect size to colour by
# nullcol is colour to plot if data doesn't meet ESlevel
# negcol is colour to plot if data meets -ESlevel 
# poscol is colour to plot if data meets +ESlevel 
# cpgcol is colour for cpg island track
# pos.shaded.regions should be supplied as a string of positions for which to shade the plot.
# first number = start of shade, second = end of shade. If multiple shaded regions, add these as extra values, keeping with the first-second, start-end rule.

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Homo.sapiens)
require(AnnotationHub)
ESplot <- TRUE
if(missing(data))        {
stop("please provide data")
                }	

if(is.null(xlab)) {
	xlab = "Genomic Position"
	}

if(is.null(ylab)) {
	ylab = "-log10(p)"
	}	
	
if(is.null(main)) {
	main = ""
	}

if(is.null(chrcol)){
    chrcol="CHR"
    }
    
if(is.null(mapcol)){
    mapcol="MAPINFO"
    }    
      
if(missing(chr)) {
                        
if(missing(cpg)){

stop("please provide either 'chr' and 'range' or 'cpg'")
                
                }	
		          }

if(! is.null(cpg)){

	data[cpg,chrcol]->chr
	data[cpg,mapcol]->MI
	range = MI:MI
	
                  }	
if(verbose){print(pad)}
if(is.null(xlim)) {
      xlim <- c((min(range) - pad[1]),(max(range) + pad[2]))
                  }

if(missing(result)){

stop("please provide result you would like to plot")
                
                   }	

par(xpd=TRUE)

data[which(data[,chrcol] == chr & data[,mapcol] >= min(range) - 
pad[1] & data[,mapcol] <= max(range) + pad[2]),]->y

print(colour.below.Pthresh)

yval <- y[,result] # define y-axis values
pvals <- y[,Pcol] # define p-values
sigprobes <- which(pvals<P.thresh) # identify significant probes
nonsigprobes <- which(pvals>P.thresh) # identify non-significant probes
if(!is.null(ESdat)){  effSz <- y[,ESdat]  } # define effect size if applicable
if(result==Pcol){ # if y-axis data equal to p-value data:
	YequalsP <- TRUE # flag
	yval <- -log10(yval) # take -log10(yvals)
}else{
	YequalsP <- FALSE
}


if(is.null(ylim)) {
	ylim <- c(min(yval),max(yval))
	}
	
if(is.null(chr)) {
	chr = y[,chrcol]
	}
	
if(is.null(cexgene)) {
	cexgene = 0.5
	}
	
if(is.null(geneplot)) {
	geneplot = TRUE
	}

ESlevel <- ESlevel/multiply

if(geneplot == "TRUE"){
layout(matrix(c(1,2), 2, 1, byrow = TRUE))
par(mar=c(3,3,1,2)+ 0.1, mgp=c(2,1,0)) 


if(!is.null(ESdat)){

	if(colour.below.Pthresh==FALSE){   
	# plot non-significant points  
	plot(yval[nonsigprobes]~y[nonsigprobes,mapcol], xlim = xlim, ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, col=nonsigcol, cex.axis = cex.axis, cex.lab = cex.lab, main=main)

	# plot significant points that have small ES
	points(yval[sigprobes][which(y[sigprobes,ESdat] > -ESlevel & y[sigprobes,ESdat] < ESlevel)]~y[sigprobes,mapcol][which(y[sigprobes,ESdat] > -ESlevel & y[sigprobes,ESdat] < ESlevel)], col=nullcol, cex = cex)

	# significant and +ve ES
	points(yval[sigprobes][which(y[sigprobes,ESdat] >= ESlevel)]~y[sigprobes,mapcol][which(y[sigprobes,ESdat] >= ESlevel)], col=poscol, cex = 1.2, pch=16)

	# significant and -ve ES
	points(yval[sigprobes][which(y[sigprobes,ESdat] <= -ESlevel)]~y[sigprobes,mapcol][which(y[sigprobes,ESdat] <= -ESlevel)], col=negcol, cex = 1.2, pch=16)

		if(linePlot){
			d <- data.frame(X=y[,mapcol], Y=yval)
			d.ord <- d[order(d$X),] #order from smallest position to largest
			lines(d.ord$Y ~ d.ord$X)
		}

	}
	else{
		
		if(sigProbesOnly){
			if(verbose){print("Only plotting significant probes")}
			# +ve ES
			plot(yval[sigprobes][which(effSz[sigprobes] >= ESlevel)]~y[sigprobes,mapcol][which(effSz[sigprobes] >= ESlevel)], col=poscol, cex = cex, xlim = xlim, ylim = ylim, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, cex.axis = cex.axis, cex.lab = cex.lab, main=main)
			# -ve ES
			points(yval[sigprobes][which(effSz[sigprobes] <= -ESlevel)]~y[sigprobes,mapcol][which(effSz[sigprobes] <= -ESlevel)], col=negcol, cex = cex)
			
			if(linePlot){
				d <- data.frame(X=y[sigprobes,mapcol], Y=yval[sigprobes])
				d.ord <- d[order(d$X),] #order from smallest position to largest
				lines(d.ord$Y ~ d.ord$X)
			}
			
		}else{
			if(verbose){print("Plotting all probes")}
			plot(yval[which(effSz > -ESlevel & effSz < ESlevel)]~y[which(effSz > -ESlevel & effSz < ESlevel),mapcol], xlim = xlim, ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, col=nullcol, cex.axis = cex.axis, cex.lab = cex.lab, main=main)
			
			# +ve ES
			points(yval[which(effSz >= ESlevel)]~y[which(effSz >= ESlevel),mapcol], col=poscol, cex = cex)
			# -ve ES
			points(yval[which(effSz <= -ESlevel)]~y[which(effSz <= -ESlevel),mapcol], col=negcol, cex = cex)
			
			if(linePlot){
				d <- data.frame(X=y[,mapcol], Y=yval)
				d.ord <- d[order(d$X),] #order from smallest position to largest
				lines(d.ord$Y ~ d.ord$X)
			}
		}
	}
	
	if(YequalsP){
		abline(h=-log10(P.thresh),lty=3,xpd = FALSE)
	}else{
		abline(h=0,lty=3,xpd = FALSE)
	}

	if(verbose){print(paste0("pos.shaded.regions=",pos.shaded.regions))}

	if(shadeDMR==TRUE){
		ShadeCol <- rgb(255, 255, 0, max = 255, alpha = 40, names = "paleyellow")
		if(is.null(pos.shaded.regions)){
			rect(min(range), -1e6, max(range), 1e6, col=ShadeCol, border=NA, xpd=F)
		}else{
			index.regions.start = seq(from=1,to=length(pos.shaded.regions), by=2)
			index.regions.end <- seq(2,length(pos.shaded.regions), 2) 
			if(shaded.buffer){
				shaded.range <- pos.shaded.regions[index.regions.end] - pos.shaded.regions[index.regions.start]
				if(!is.null(manual.buffer)){
					suppressWarnings({
						pos.shaded.regions[index.regions.start] <- pos.shaded.regions[index.regions.start] - manual.buffer
						pos.shaded.regions[index.regions.end] <- pos.shaded.regions[index.regions.end] + manual.buffer
					})
					if(verbose){print(paste0("pos.shaded.regions=",pos.shaded.regions))}
				}else{
					if(verbose){print(shaded.range)}
					buffer <- (buffer.perc/100)*shaded.range
					if(verbose){print(buffer)}
					suppressWarnings({
						pos.shaded.regions[index.regions.start] <- pos.shaded.regions[index.regions.start] - buffer
						pos.shaded.regions[index.regions.end] <- pos.shaded.regions[index.regions.end] + buffer
					})
				}
			}
			if(verbose){print(paste('plotting',length(pos.shaded.regions)/2,'shaded points'))}
			rect(pos.shaded.regions[index.regions.start], -100, pos.shaded.regions[index.regions.end], 100, col=ShadeCol, border=NA, xpd=F)
		}
	}

}

print('Starting tracks')
chr2=paste("chr", chr, sep="")

hub <- AnnotationHub()
query(hub, c("cpg","hg19"))
cpgs <- as.data.frame(hub[["AH5086"]])
cpgs[which(cpgs$seqnames == chr2 & cpgs$start >=min(range)-pad[1] & 
cpgs$start <= max(range)+pad[2]),]->cpg
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gr <- GRanges(seqnames = chr2, ranges = IRanges(start = min(range)-pad[1], 
end = max(range)+pad[2]))
print('made GRanges')
subsetByOverlaps(transcripts(txdb), gr)->trans

res <- as.data.frame(transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL")))
as.matrix(res[res$TXNAME %in% as.data.frame(trans)[,"tx_name"],"SYMBOL"])->m
colnames(m)<-"SYMBOL"
cbind(as.data.frame(trans),m)->trans
trans[is.na(trans$SYMBOL) == "FALSE",]->trans
if(verbose){print(trans)}
if(verbose){print('Trans 1')}

subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)->exons
trans[rownames(trans) %in% names(exons),]->trans
if(verbose){print('Trans 2')}
if(is.null(genelines)){
	genelines = nrow(trans)
	                  }

rep(1:max(genelines), times=ceiling(nrow(trans)/max(genelines)))->gl
gl[1:nrow(trans)]->gl
trans$gl <- gl

exons[rownames(trans)]->exons
if(verbose){print('Trans 3')}
par(mar=c(1,3,1,2)+ 0.1) 
plot(0,0,type="n", xlim=xlim, ylim=c(0,genelines +1), axes=FALSE, xlab="", 
ylab="", xpd=FALSE)

if ((nrow(cpg) >= 1)){
	for(g in 1:nrow(cpg)){
		polygon(c(cpg[g, "start"], cpg[g, "end"], cpg[g, "end"], cpg[g, "start"]), 
		c((genelines+1)-0.2, (genelines+1)-0.2, (genelines+1)+0.2,(genelines+1)+0.2),
		col=cpgcol, xpd=FALSE, border=cpgcol)
	}
}

for(i in 1:nrow(trans)){

	lines(c(as.data.frame(trans)[i,"start"], as.data.frame(trans)[i,"end"]), 
	c(trans$gl[i],trans$gl[i]), xpd=FALSE)
	text(as.data.frame(trans)[i,"start"]+((as.data.frame(trans)[i,"end"] - 
	as.data.frame(trans)[i,"start"])/2) ,trans$gl[i] + 0.4 , trans$SYMBOL[i], 
	cex=cexgene, xpd=FALSE)

	for(x in 1:nrow(as.data.frame(exons[[i]]))){

		polygon(c(as.data.frame(exons[[i]])[x,"start"], 
		as.data.frame(exons[[i]])[x,"end"], as.data.frame(exons[[i]])[x,"end"], 
		as.data.frame(exons[[i]])[x,"start"]), 
		c(trans$gl[i]-0.2, trans$gl[i]-0.2, trans$gl[i]+0.2,trans$gl[i]+0.2),
		col="black", xpd=FALSE)
	}
}

}

if(geneplot == "FALSE"){
	plot.new()
	if(!is.null(ESdat)){
		par(mar=c(3,3,1,2)+ 0.1, mgp=c(2,1,0)) 

		plot(-log10(y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),result])~
		y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),mapcol], xlim = xlim, 
		ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, 
		col=nullcol, cex.axis = cex.axis, cex.lab = cex.lab)
		points(-log10(y[which(y[,ESdat] >= ESlevel),result])~
		y[which(y[,ESdat] >= ESlevel),mapcol], col=poscol, cex = cex)
		points(-log10(y[which(y[,ESdat] <= -ESlevel),result])~
		y[which(y[,ESdat] <= -ESlevel),mapcol], col=negcol, cex = cex)

	}else{
		plot(-log10(y[,result])~y[,mapcol], xlim = xlim, ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, col=nullcol, cex.axis = cex.axis, cex.lab = cex.lab)
	}
}

}

