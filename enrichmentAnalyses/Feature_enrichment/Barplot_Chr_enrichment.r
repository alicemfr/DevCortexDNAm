

## Enrichment of chromosomes within EWAS results ##

library(data.table)
library(dplyr)


resultsFile <- paste0(AnalysisPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")
res <- readRDS(resultsFile)
colP <- 'P.Age'
colBeta <- 'Beta.Age'


#1. Load results & annotate =====================================================================================================
if(!any(grepl('CHR|chr|Chr', colnames(res)))){
	epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE)
	epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),'CHR']
	res <- cbind(res, as.data.frame(epicMan))
}else{
	# if chromosome col not already called 'CHR', change
	colnames(res)[grep('chr|Chr', colnames(res))] <- 'CHR'
}

if(any(res$CHR %in% 'Y')){ res <- res[-which(res$CHR=='Y'),] }   # remove Y chr
res$CHR[which(res$CHR=='X')] <- '23'                             # change chr X to chr 23
if(any(is.na(res$CHR))){ res <- res[-which(is.na(res$CHR)),] }   # remove NA chrs

res.sig <- res[which(res[,colP]<9e-8),]	                         # significant results
res.nonsig <- res[-which(rownames(res) %in% rownames(res.sig)),] # non-significant results


#2. Calculate statistics for chr distribution ===================================================================================

ChrCounts <- function(mat, HypoHyper=FALSE){
	chr_feature <- mat$CHR									     # chr in results
	chr.feature <- 1:23										     # names of possible chrs
	chr_feature_counts <- matrix(data=NA, ncol=length(chr.feature), nrow=1)
	colnames(chr_feature_counts) <- chr.feature
	for(i in 1:length(chr.feature)){
		ln <- length(which(chr_feature==chr.feature[i]))	     # number for each possible chr in results
		chr_feature_counts[1,i] <- ln
	}
	
	# if further splitting into hypo- or hyper-methylated probes
	if(HypoHyper==TRUE){
		chr_feature_hypohyper_counts <- matrix(data=NA, ncol=length(chr.feature), nrow=2)
		colnames(chr_feature_hypohyper_counts) <- chr.feature
		rownames(chr_feature_hypohyper_counts) <- c('Hypo','Hyper')
		for(i in 1:length(chr.feature)){
			mat.beta <- mat[which(chr_feature==chr.feature[i]),colBeta]    # extract each chr and corresponding effect size
			hypo <- length(which(mat.beta<0))
			hyper <- length(which(mat.beta>0))
			chr_feature_hypohyper_counts[1:2,i] <- c(hypo,hyper)
		}
	}else{
		chr_feature_hypohyper_counts=NULL
	}
	
	return(list(Chr.Features=chr_feature_counts, Chr.HypoHyper.Features=chr_feature_hypohyper_counts))
}

sig.counts <- ChrCounts(res.sig, HypoHyper=T)  # DMPs
nonsig.counts <- ChrCounts(res.nonsig)         # non-DMPs
total.counts <- ChrCounts(res)                 # all array sites


#3a. Barplot of chr distribution for all probes =================================================================================
background.sig.chr <- nrow(res.sig)-sig.counts$Chr.Features                 # number of significant probes that are not in that chr
background.nonsig.chr <- nrow(res.nonsig)-nonsig.counts$Chr.Features        # number of non-significant probes that are not in that chr

Total.probes <- nrow(res)
perc.chr.counts <- (sig.counts$Chr.Features/total.counts$Chr.Features)*100
perc.sig <- (nrow(res.sig)/Total.probes)*100                                # background significance rate across all chrs
counts <- c(perc.sig, as.numeric(perc.chr.counts))


#3b. Calculate enrichment statistics ============================================================================================
enr_test <- function(features, resCol, sig.counts, nonsig.counts, background.sig, background.nonsig){
	enr <- matrix(NA,ncol=2, nrow=length(features))
	rownames(enr) <- features
	colnames(enr) <- c('P','Beta')
		for(i in 1:length(features)){
			f <- features[i] # chr of interest
			sig.f <- sig.counts[[resCol]][which(colnames(sig.counts[[resCol]])==f)]	           # significant counts of this chr
			nonsig.f <- nonsig.counts[[resCol]][which(colnames(nonsig.counts[[resCol]])==f)]   # non-significant counts of this chr
			sig.o <- background.sig[which(colnames(background.sig)==f)]	                       # significant counts of other chrs
			nonsig.o <- background.nonsig[which(colnames(background.nonsig)==f)]               # non-significant counts of other chrs
			
			m <- matrix( c( c(sig.f,sig.o),c(nonsig.f,nonsig.o) ), nrow=2, byrow=T)
			chi <- chisq.test(m)
			obs <- chi$observed
			expd <- chi$expected
			ES <- obs/expd	                        # effect size is observed/expected
			ES[which(ES<1)] <- -(1/ES[which(ES<1)])	# for depletions, take the reciprocal and convert to negative. Easier to interpret.
			p <- chi$p.value
			enr[i,'P'] <- p
			enr[i,'Beta'] <- ES[1,1]
		}
	return(enr)
}

gene_enr <- enr_test(features=1:23, resCol='Chr.Features', sig.counts=sig.counts, nonsig.counts=nonsig.counts, background.sig=background.sig.chr, background.nonsig=background.nonsig.chr)


#3c. Plot =======================================================================================================================

# find which chrs are significant
sigs <- gene_enr[,'P']
sigs.adj <- sigs*length(sigs)                           # correct for number of tests performed
if(any(sigs.adj>1)){sigs.adj[which(sigs.adj>1)] <- 1}   # set any corrected p-values that are >1 to 1
sigf <- sigs.adj<5e-2                                   # T/F whether test is significant (corrected p < 0.05)

# colour bars by magnitude of direction of effect
gene_enr <- as.data.frame(gene_enr)
gene_enr$Chr <- rownames(gene_enr)
gene_enr$Sig <- sigf
gene_enr$result <- rep('Not significant',nrow(gene_enr))					
gene_enr$result[which(gene_enr$Sig==T & gene_enr$Beta>1)] <- 'Over-enriched'		
gene_enr$result[which(gene_enr$Sig==T & gene_enr$Beta<1)] <- 'Under-enriched'
gene_enr <- arrange(gene_enr, result)
colours <- c(
	rep('lightgrey',length(which(gene_enr$result=='Not significant'))), 
	rep('firebrick',length(which(gene_enr$result=='Over-enriched'))),
	rep('steelblue3',length(which(gene_enr$result=='Under-enriched')))
	)
gene_enr$colours <- colours
gene_enr <- gene_enr[match(1:23,gene_enr$Chr),]

# plot
pdf(paste0(plotPath, plotFilename, '.pdf'), width=12,height=8)
par(mar=c(7, 5, 4.1, 2.1))
barplot(counts, main="", ylab='Significant probes (%)', names.arg=c('All', 1:22, "X"), xlab='Chromosome', ylim=c(0,(max(counts)+0.1)), cex.lab=1.8, cex.names=1.25, cex.axis=1.5, col=c('white',gene_enr$colours))
abline(h=counts[1], lty=3)
dev.off()

# Save results
gene_enr$Perc <- as.numeric(perc.chr.counts)
gene_enr$colours <- NULL
gene_enr$Chr <- c(1:22,'X')
gene_enr$P.adj <- as.numeric(sigs.adj)
gene_enr$Total.sites <- as.numeric(total.counts$Chr.Features)
gene_enr$Sig.sites <- as.numeric(sig.counts$Chr.Features)
gene_enr <- gene_enr[,c(which(colnames(gene_enr)=='Chr'), which(colnames(gene_enr)=='Sig.sites'), which(colnames(gene_enr)=='Total.sites'), which(colnames(gene_enr)=='Perc'), which(colnames(gene_enr) %in% c('P','Beta')), which(colnames(gene_enr)=='P.adj'))]
if(any(gene_enr$Perc==0)){
	gene_enr[which(gene_enr$Perc==0),c('P','Beta','P.adj')] <- NA
}
write.csv(gene_enr, file=paste0(resTablePath, resTableFilename, ".csv"))


#4a. Barplot of chr distribution for significant probes  - split by hypo-/hyper-methylated ======================================
chr.feature <- 1:23

HypoHyper.chr <- matrix(data=NA,ncol=length(chr.feature), nrow=2)
colnames(HypoHyper.chr) <- chr.feature
rownames(HypoHyper.chr) <- c('Hypo','Hyper')
HypoHyper.chr['Hypo',] <- (sig.counts$Chr.HypoHyper.Features['Hypo',]/sig.counts$Chr.Features)*100
HypoHyper.chr['Hyper',] <- (sig.counts$Chr.HypoHyper.Features['Hyper',]/sig.counts$Chr.Features)*100

perc.sig.hypo = (length(which(res.sig[,colBeta]<0))/nrow(res.sig))*100
perc.sig.hyper = (length(which(res.sig[,colBeta]>0))/nrow(res.sig))*100

counts <- cbind(Allprobes=c(perc.sig.hypo,perc.sig.hyper),HypoHyper.chr)
# change order of rows so figure shows hyper then hypo
counts.swap <- counts
counts.swap[1,] <- counts.swap[2,]
counts.swap[2,] <- counts[1,]
rownames(counts.swap) <- c('Hyper','Hypo')

counts.swap <- counts.swap[,!is.nan(counts.swap[1,])] # remove any cols with NaN values


#4b. Calculate enrichment statistics ============================================================================================

# binomial test for All Probes hypo vs hyper enrichment
res.dir <- sign(res.sig[,colBeta])
lgrp <- names(table(res.dir))[which(table(res.dir)==max(table(res.dir)))] # most common direction
allprb <- binom.test(length(which(res.dir==lgrp)), length(res.dir), alternative='greater')$p.value

prop.hypo <- length(which(res.dir==lgrp))/length(res.dir) # proportion of hypomethylated probes across all chromosomes

enr_test <- function(features, res, prob){
	enr <- c()
		for(i in 1:length(features)){
			f <- features[i] # chr of interest
			tot.chr <- res$Chr.Features[,which(colnames(res$Chr.Features)==f)] # extract total n probes in chr
			tot.chr <- as.numeric(tot.chr)
			if(tot.chr==0){ # if no DMPs for this feature, skip
				enr <- c(enr, NA)
				next
			}
			hh.chr <- res$Chr.HypoHyper.Features[,which(colnames(res$Chr.HypoHyper.Features)==f)] # hypo/hypermethylated distribution for chr
			hypo.chr <- as.numeric(hh.chr['Hypo']) # extract n hypo probes in chr
			
			binom <- binom.test(hypo.chr, tot.chr, p=prob, alternative='greater') # prob = hypothesised probability of success, e.g. 0.5
			enr <- c(enr,binom$p.value)
		}
	return(enr)
}

chr_enr.50 <- enr_test(features=1:23, res=sig.counts, prob=0.5)    # probability set to 50:50
chr_enr <- enr_test(features=1:23, res=sig.counts, prob=prop.hypo) # probability set to proportion hypomethylated across all chromosomes combined


#4c. Plot =======================================================================================================================
enr <- chr_enr
prob <- prop.hypo
test <- "bg" # probability = background rate

enr <- chr_enr.50
prob <- 0.5
test <- "0.5" # probability = 50:50

# find which chrs are significant
sigs <- c(allprb, enr)                          # combine all p-values
sigs.adj <- sigs*length(sigs)                   # correct for number of tests performed
if(any(sigs.adj[!is.na(sigs.adj)]>1)){sigs.adj[!is.na(sigs.adj)][which(sigs.adj[!is.na(sigs.adj)]>1)] <- 1}
sigf <- sigs.adj<5e-2                           # T/F whether test is significant
sigf <- sigf[which(!is.na(sigf))]               # NA when no DMPs annot to that chr

# change chr 23 to chr X
colnames(counts.swap)[colnames(counts.swap)==23] <- 'X'

# colours based on significance
col.sig <- c('red','turquoise3')
col.nonsig <- c('rosybrown1','lightblue1')
col_list <- list()
for(i in 1:length(sigf)){
	if(sigf[i]==TRUE){
		col_list[[i]] <- col.sig
	}else{
		col_list[[i]] <- col.nonsig
	}
}
cols <- unlist(col_list)
plotnames <- c('All', paste(1:22)[1:22 %in% colnames(counts.swap)], paste('X')['X' %in% colnames(counts.swap)]) # for boxplot, keep appropriate chrs

# Plot
pdf(paste0(plotPath, plotFilename_HypoHyper, "_prob=", test, ".pdf"), width=12,height=8)
par(mar=c(7, 5, 4.1, 2.1))
bxplt <- barplot(counts.swap, main="",ylab='Percentage of Significant probes (%)', names.arg=plotnames, xlab='Chromosome', ylim=c(0,109), beside=TRUE, cex.lab=1.8, cex.names=1.25, cex.axis=1.5, col=cols)
abline(h=prob*100, lty=3)
dev.off()

# Save results
enr <- as.data.frame(enr)
colnames(enr) <- 'P'
enr$Chr <- c(1:22,'X')
enr$P.adj <- as.numeric(sigs.adj[-1]) # exclude the "All probes" value in first position
enr <- cbind(enr, Perc=t(HypoHyper.chr))
enr <- enr[,c(which(colnames(enr)=='Chr'), grep('Perc',colnames(enr)), which(colnames(enr)=='P'), which(colnames(enr)=='P.adj'))]
write.csv(enr, file=paste0(resTablePath, resTableFilename_HypoHyper, "_prob=", test, ".csv"))