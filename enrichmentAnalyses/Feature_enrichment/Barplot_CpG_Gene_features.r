

## Enrichment of genic features within EWAS results ##

library(data.table)
library(dplyr)

# prenatal bulk cortex
resultsFile <- paste0(AnalysisPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")
res <- readRDS(resultsFile)
pThresh <- 9e-8
colP <- 'P.Age'
colBeta <- 'Beta.Age'

# postnatal bulk cortex
resultsFile <- paste0(AnalysisPath, "postnatalAgeReg_bulkFetalDMPs.rds")
res <- readRDS(resultsFile)
pThresh <- 1.2e-6
colP <- 'Age.P'
colBeta <- 'Age.ES'


#1. Load results & annotate =====================================================================================================
if(!any(grepl('UCSC_RefGene_Group', colnames(res))) & !any(grepl('Relation_to_UCSC_CpG_Island', colnames(res)))){ # check if required cols are present
	epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE)
	epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c("UCSC_RefGene_Name", "CHR", "UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")]
	res <- cbind(res, as.data.frame(epicMan))
}

table(res$Relation_to_UCSC_CpG_Island)
res$Relation_to_UCSC_CpG_Island[which(res$Relation_to_UCSC_CpG_Island %in% c('N_Shore','S_Shore'))] <- 'Shore'	# change all shores to Shore
res$Relation_to_UCSC_CpG_Island[which(res$Relation_to_UCSC_CpG_Island %in% c('N_Shelf','S_Shelf'))] <- 'Shelf'  # change all shelves to Shelf
table(res$Relation_to_UCSC_CpG_Island)

if(any(res$CHR %in% 'Y')){ res <- res[-which(res$CHR=='Y'),] }   # remove Y chr
res.sig <- res[which(res[,colP]<pThresh),]	                     # significant results
res.nonsig <- res[-which(rownames(res) %in% rownames(res.sig)),] # non-significant results


#2. Calculate statistics for feature distribution ===============================================================================

FeatureCounts <- function(mat, HypoHyper=FALSE){
	cpg_feature <- mat$Relation_to_UCSC_CpG_Island                                  # cpg features in results
	gene_feature <- mat$UCSC_RefGene_Group                                          # gene features in results

	# gene feature matrix
	gene.feature = c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR")           # names of possible gene features
	gene_feature_counts <- matrix(data=NA, ncol=length(gene.feature), nrow=1)
	colnames(gene_feature_counts) <- gene.feature
	for(i in 1:length(gene.feature)){
		ln <- length(grep(gene.feature[i], gene_feature))                           # number of entries for each possible gene.feature in results
		gene_feature_counts[1,i] <- ln
	}

	# cpg feature matrix
	cpg.feature = c("Island","Shore","Shelf")                                       # names of possible cpg features
	cpg_feature_counts <- matrix(data=NA, ncol=length(cpg.feature), nrow=1)
	colnames(cpg_feature_counts) <- cpg.feature
	for(i in 1:length(cpg.feature)){
		ln <- length(grep(cpg.feature[i], cpg_feature))                             # number of entries for each possible cpg.feature in results
		cpg_feature_counts[1,i] <- ln
	}
	
	# if further splitting into hypo- or hyper-methylated probes
	if(HypoHyper==TRUE){
		gene_feature_hypohyper_counts <- matrix(data=NA, ncol=length(gene.feature), nrow=2)
		colnames(gene_feature_hypohyper_counts) <- gene.feature
		rownames(gene_feature_hypohyper_counts) <- c('Hypo','Hyper')
		for(i in 1:length(gene.feature)){
			mat.beta <- mat[grep(gene.feature[i], gene_feature),colBeta]            # extract each gene.feature and corresponding effect size
			hypo <- length(which(mat.beta<0))
			hyper <- length(which(mat.beta>0))
			gene_feature_hypohyper_counts[1:2,i] <- c(hypo,hyper)
		}
		
		cpg_feature_hypohyper_counts <- matrix(data=NA, ncol=length(cpg.feature), nrow=2)
		colnames(cpg_feature_hypohyper_counts) <- cpg.feature
		rownames(cpg_feature_hypohyper_counts) <- c('Hypo','Hyper')
		for(i in 1:length(cpg.feature)){
			mat.beta <- mat[grep(cpg.feature[i], cpg_feature),colBeta]              # extract each cpg.feature and corresponding effect size
			hypo <- length(which(mat.beta<0))
			hyper <- length(which(mat.beta>0))
			cpg_feature_hypohyper_counts[1:2,i] <- c(hypo,hyper)
		}
		mat.beta.other <- mat[which(cpg_feature==''),colBeta]
	}else{gene_feature_hypohyper_counts=NULL; cpg_feature_hypohyper_counts=NULL}
	
	return(list(Gene.Features=gene_feature_counts, CpG.Features=cpg_feature_counts, Gene.HypoHyper.Features=gene_feature_hypohyper_counts, CpG.HypoHyper.Features=cpg_feature_hypohyper_counts))
}

sig.counts <- FeatureCounts(res.sig, HypoHyper=T)  # DMPs
nonsig.counts <- FeatureCounts(res.nonsig)         # non-DMPs
total.counts <- FeatureCounts(res)                 # all array sites


#3a. Barplot of CpG/Gene feature distribution for all probes ====================================================================
background.sig.gene <- nrow(res.sig)-sig.counts$Gene.Features				# Genic: number of significant probes that are not part of that feature
background.nonsig.gene <- nrow(res.nonsig)-nonsig.counts$Gene.Features		# Genic: number of non-significant probes that are not part of that feature
background.sig.cpg <- nrow(res.sig)-sig.counts$CpG.Features					# CpG: 	 number of significant probes that are not part of that feature
background.nonsig.cpg <- nrow(res.nonsig)-nonsig.counts$CpG.Features		# CpG: 	 number of non-significant probes that are not part of that feature

Total.probes <- nrow(res)
perc.gene.counts <- (sig.counts$Gene.Features/total.counts$Gene.Features)*100
perc.cpg.counts <- (sig.counts$CpG.Features/total.counts$CpG.Features)*100
perc.sig = (nrow(res.sig)/Total.probes)*100                                 # background significance rate across all sites
counts <- c(perc.sig, as.numeric(perc.cpg.counts), as.numeric(perc.gene.counts))


#3b. Calculate enrichment statistics ============================================================================================
enr_test <- function(features, resCol, sig.counts, nonsig.counts, background.sig, background.nonsig){
	enr <- matrix(NA,ncol=2, nrow=length(features))
	rownames(enr) <- features
	colnames(enr) <- c('P','Beta')
		for(i in 1:length(features)){
			f <- features[i]                                                                   # feature of interest
			sig.f <- sig.counts[[resCol]][which(colnames(sig.counts[[resCol]])==f)]            # significant counts of this feature
			nonsig.f <- nonsig.counts[[resCol]][which(colnames(nonsig.counts[[resCol]])==f)]   # non-significant counts of this feature
			sig.o <- background.sig[which(colnames(background.sig)==f)]                        # significant counts of other features
			nonsig.o <- background.nonsig[which(colnames(background.nonsig)==f)]               # non-significant counts of other features
			
			m <- matrix( c( c(sig.f,sig.o),c(nonsig.f,nonsig.o) ), nrow=2, byrow=T)
			chi <- chisq.test(m)
			obs <- chi$observed
			expd <- chi$expected
			ES <- obs/expd                           # effect size is observed/expected
			ES[which(ES<1)] <- -(1/ES[which(ES<1)])  # for depletions, take the reciprocal and convert to negative. Easier to interpret.
			p <- chi$p.value
			enr[i,'P'] <- p
			enr[i,'Beta'] <- ES[1,1]
		}
	return(enr)
}

gene_enr.gene <- enr_test(features=c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR"), resCol='Gene.Features', sig.counts=sig.counts, nonsig.counts=nonsig.counts, background.sig=background.sig.gene, background.nonsig=background.nonsig.gene)

gene_enr.cpg <- enr_test(features=c("Island","Shore","Shelf"), resCol='CpG.Features', sig.counts=sig.counts, nonsig.counts=nonsig.counts, background.sig=background.sig.cpg, background.nonsig=background.nonsig.cpg)

gene_enr <- rbind(gene_enr.cpg, gene_enr.gene)


#3c. Plot =======================================================================================================================

# find which features are significant
sigs <- gene_enr[,'P']
sigs.adj <- sigs*length(sigs)                          # correct for number of tests performed
if(any(sigs.adj>1)){sigs.adj[which(sigs.adj>1)] <- 1}  # set any corrected p-values that are >1 to 1
sigf <- sigs.adj<5e-2                                  # T/F whether test is significant (corrected p < 0.05)

# colour bars by magnitude of direction of effect
gene_enr <- as.data.frame(gene_enr)
gene_enr$Feature <- rownames(gene_enr)
gene_enr$Sig <- sigf
gene_enr$result <- rep('Not significant',nrow(gene_enr))					
gene_enr$result[which(gene_enr$Sig==T & gene_enr$Beta>1)] <- 'Over-enriched'		
gene_enr$result[which(gene_enr$Sig==T & gene_enr$Beta<1)] <- 'Under-enriched'	
fts <- gene_enr$Feature
gene_enr <- arrange(gene_enr, result)
colours <- c(
	rep('lightgrey',length(which(gene_enr$result=='Not significant'))), 
	rep('firebrick',length(which(gene_enr$result=='Over-enriched'))),
	rep('steelblue3',length(which(gene_enr$result=='Under-enriched')))
	)
gene_enr$colours <- colours
gene_enr <- gene_enr[match(fts,gene_enr$Feature),]

# plot
pdf(paste0(plotPath, plotFilename, '.pdf'))
par(mar=c(7, 5, 4.1, 2.1))
barplot(counts, main="",ylab='Significant probes (%)', names.arg=c('All probes', 'Island','Shore','Shelf','TSS1500','TSS200',"5'UTR",'1stExon','Gene body',"3'UTR"), ylim=c(0,(max(counts)+0.05)), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=c('white',gene_enr$colours))
abline(h=counts[1], lty=3)
dev.off()

# Save results
gene_enr <- as.data.frame(gene_enr)
gene_enr$Feature <- rownames(gene_enr)
gene_enr$P.adj <- as.numeric(sigs.adj)
gene_enr$Perc <- c(as.numeric(perc.cpg.counts), as.numeric(perc.gene.counts))
gene_enr$Total.sites <- c(total.counts$CpG.Features, total.counts$Gene.Features)
gene_enr$Sig.sites <- c(sig.counts$CpG.Features, sig.counts$Gene.Features)
gene_enr <- gene_enr[,c(which(colnames(gene_enr)=='Feature'), which(colnames(gene_enr)=='Sig.sites'), which(colnames(gene_enr)=='Total.sites'), which(colnames(gene_enr)=='Perc'), which(colnames(gene_enr) %in% c('P','Beta')),which(colnames(gene_enr)=='P.adj'))]
write.csv(gene_enr, file=paste0(resTablePath, resTableFilename, ".csv"))


#4a. Barplot of CpG/Gene feature distribution for significant probes  - split by hypo-/hyper-methylated =========================
gene.feature <- c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR")
cpg.feature <- c("Island","Shore","Shelf")

HypoHyper.gene <- matrix(data=NA,ncol=length(gene.feature), nrow=2)
colnames(HypoHyper.gene) <- gene.feature
rownames(HypoHyper.gene) <- c('Hypo','Hyper')
HypoHyper.gene['Hypo',] <- (sig.counts$Gene.HypoHyper.Features['Hypo',]/sig.counts$Gene.Features)*100
HypoHyper.gene['Hyper',] <- (sig.counts$Gene.HypoHyper.Features['Hyper',]/sig.counts$Gene.Features)*100

HypoHyper.cpg <- matrix(data=NA,ncol=length(cpg.feature), nrow=2)
colnames(HypoHyper.cpg) <- cpg.feature
rownames(HypoHyper.cpg) <- c('Hypo','Hyper')
HypoHyper.cpg['Hypo',] <- (sig.counts$CpG.HypoHyper.Features['Hypo',]/sig.counts$CpG.Features)*100
HypoHyper.cpg['Hyper',] <- (sig.counts$CpG.HypoHyper.Features['Hyper',]/sig.counts$CpG.Features)*100

perc.sig.hypo = (length(which(res.sig[,colBeta]<0))/nrow(res.sig))*100
perc.sig.hyper = (length(which(res.sig[,colBeta]>0))/nrow(res.sig))*100

counts <- cbind(Allprobes=c(perc.sig.hypo,perc.sig.hyper),HypoHyper.cpg, HypoHyper.gene)
#swap hypo and hyper rows
counts.swap <- counts
counts.swap[1,] <- counts.swap[2,]
counts.swap[2,] <- counts[1,]
rownames(counts.swap) <- c('Hyper','Hypo')


#4b. Calculate enrichment statistics ============================================================================================

# binomial test for All Probes hypo vs hyper enrichment
res.dir <- sign(res.sig[,colBeta])
lgrp <- names(table(res.dir))[which(table(res.dir)==max(table(res.dir)))] # most common direction
allprb <- binom.test(length(which(res.dir==lgrp)), length(res.dir), alternative='greater')$p.value

prop.hypo <- length(which(res.dir==lgrp))/length(res.dir) # proportion of hypomethylated probes across all chromosomes

enr_test <- function(features, res, prob, type){
	enr <- c()
		for(i in 1:length(features)){
			f <- features[i]
			tot.f <- res[[paste(type,'Features',sep='.')]][,which(colnames(res[[paste(type,'Features',sep='.')]])==f)] # extract total n probes in feature
			tot.f <- as.numeric(tot.f)
			if(tot.f==0){ # if no DMPs for this feature, skip
				enr <- c(enr, NA)
				next
			}
			hh.f <- res[[paste(type,'HypoHyper.Features',sep='.')]][,which(colnames(res[[paste(type,'HypoHyper.Features',sep='.')]])==f)] # hypo/hypermethylated distribution for feature
			dir.f <- hh.f[which(hh.f==max(hh.f))] # extract n probes for most common direction
			if(names(dir.f)=='Hyper'){ Prob <- 1-prob }else{ Prob <- prob }
			binom <- binom.test(as.numeric(dir.f), tot.f, p=Prob, alternative='greater') # prob = hypothesised probability of success, e.g. 0.5
			enr <- c(enr,binom$p.value)
		}
	return(enr)
}

# probability assumes 50:50 hypo:hyper
cpg_enr.50 <- enr_test(features=cpg.feature, res=sig.counts, prob=0.5, type='CpG')
gene_enr.50 <- enr_test(features=gene.feature, res=sig.counts, prob=0.5, type='Gene')

# probability set to proportion hypomethylated across all sites combined
cpg_enr <- enr_test(features=cpg.feature, res=sig.counts, prob=prop.hypo, type='CpG')
gene_enr <- enr_test(features=gene.feature, res=sig.counts, prob=prop.hypo, type='Gene')


#4c. Plot =======================================================================================================================
enr <- c(cpg_enr, gene_enr)
prob <- prop.hypo
test <- "bg"  # probability = background rate

enr <- c(cpg_enr.50, gene_enr.50)
prob <- 0.5
test <- "0.5" # probability = 50:50

# find which features are significant
sigs <- c(allprb, enr)                          # combine all p-values
sigs.adj <- sigs*length(sigs)                   # correct for number of tests performed
if(any(sigs.adj>1)){sigs.adj[which(sigs.adj>1)] <- 1}
sigf <- sigs.adj<5e-2                           # T/F whether test is significant
sigf <- sigf[which(!is.na(sigf))]               # NA when no DMPs annot to that feature

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

feat <- c('Island','Shore','Shelf', 'TSS1500','TSS200',"5'UTR",'1stExon','Body',"3'UTR")
feat.lab <- c('Island','Shore','Shelf', 'TSS1500','TSS200',"5'UTR",'1stExon','Gene body',"3'UTR")
plotnames <- c('All probes', paste(feat.lab)[feat %in% colnames(counts.swap)])

# plot
pdf(paste0(plotPath, plotFilename_HypoHyper, "_prob=", test, ".pdf"))
par(mar=c(7, 5, 4.1, 2.1))
barplot(counts.swap, main="",ylab='Percentage of Significant probes (%)', names.arg=plotnames, ylim=c(0,108), las=2, beside=TRUE, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=cols)
abline(h=prob*100, lty=3)
legend(x=20,y=108, legend=c('Hypermethylated','Hypomethylated'), col=c('red','turquoise3'), pch=15, cex=1.2, bty='n')
dev.off()

# Save results
enr <- as.data.frame(enr)
names(enr) <- 'P'
enr$Feature <- c(cpg.feature,gene.feature)
enr$P.adj <- as.numeric(sigs.adj[-1]) # exclude the "All probes" value
enr <- cbind(enr, Perc=rbind(t(HypoHyper.cpg), t(HypoHyper.gene)))
enr <- enr[,c(which(colnames(enr)=='Feature'), grep('Perc',colnames(enr)), which(colnames(enr)=='P'), grep('Beta',colnames(enr)), which(colnames(enr)=='P.adj'))]
write.csv(enr, file=paste0(resTablePath, resTableFilename_HypoHyper, "_prob=", test, ".csv"))

