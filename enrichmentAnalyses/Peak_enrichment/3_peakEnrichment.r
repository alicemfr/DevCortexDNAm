

## Perform logistic regression to test enrichment of cell-type peaks in DMPs, whilst controlling for contributions from other cell-types ##


library(dplyr)

cells <- readRDS(paste0(refPath,"tissueNames_peakEnrichment.rds"))


#1. Load EWAS results annotated to cell-type-specific peaks =====================================================================

# Run either...

# bulk fetal
res <- readRDS("ageReg_fetalBrain_annot_allTissuePeaks.rds") 
outFile <- "fetalBulk_peakEnrichment_logRegStats.csv"
celltype <- ''

# FANS fetal
res <- readRDS("FANS_AgeCellSpecific_annot_allTissuePeaks.rds")
celltype <- 'Non.neuronal'
res$DMP <- res$DMP.Neuronal==FALSE & res$DMP.Non.neuronal==TRUE # change as appropriate to isolate neuronal or non-neuronal DMPs
outFile <- paste0("FANS_",celltype,"specific","_peakEnrichment_logRegStats.csv")


#2. Calculate stats =============================================================================================================

# overlap of DMPs and cell-type peaks
peak_dmps <- c()
for(i in 1:length(cells)){
	peak_dmps[i] <- length(which(res[,cells[i]]==TRUE & res$DMP==TRUE))
}

# test association between a site being a DMP and being in a peak, whilst controlling for all other cell type peaks
f <- as.formula(paste("DMP", "~", paste(cells, collapse='+')))
model <- glm(f, family="binomial", data=res)
sum_mod <- summary(model)$coefficients
sum_mod <- sum_mod[-which(rownames(sum_mod)=='(Intercept)'),]


# prepare df for plotting
stats <- data.frame(Celltype=cells, DMPs.in.Peaks=peak_dmps, Estimate=sum_mod[,'Estimate'], P=sum_mod[,'Pr(>|z|)'])
stats$P.adj <- stats$P*nrow(stats)
stats$P.adj[which(stats$P.adj>1)] <- 1 # p>1 are not meaningful, cap at 1
stats$logP.adj <- -log10(stats$P.adj)
stats$OR <- exp(stats$Estimate) # calc odds ratio from estimate
stats <- stats[,c('Celltype','DMPs.in.Peaks','Estimate','OR','P','P.adj','logP.adj')]

# cell-types with zero DMPs in peaks violate assumptions of logistic regression - change stats to NA
for(i in 1:nrow(stats)){
	if(stats[i,'DMPs.in.Peaks']==0){
		stats[i,c(which(colnames(stats)=='Estimate'):ncol(stats))] <- NA
	}
}

# save output
write.csv(stats, file=outFile, row.names=F)


#3. Plot ========================================================================================================================

# colour bars by magnitude of direction of effect
stats$Sig <- factor(stats$P.adj<5e-2)
stats$result <- rep('Not significant',nrow(stats))					
stats$result[which(stats$Sig==T & stats$Estimate>0)] <- 'Over-enriched'		
stats$result[which(stats$Sig==T & stats$Estimate<0)] <- 'Under-enriched'	
fts <- stats$Celltype
stats <- arrange(stats, result)
colours <- c(
	rep('lightgrey',length(which(stats$result=='Not significant'))), 
	rep('firebrick',length(which(stats$result=='Over-enriched'))),
	rep('steelblue3',length(which(stats$result=='Under-enriched')))
	)
stats$colours <- colours
stats <- stats[match(fts,stats$Celltype),]


# barplot: effect size vs cell type

## bulk ##
pdf(paste0("Barplot_DMP",celltype,"specific","_Peak_n54_logReg_new.pdf"), width=20,height=7)
par(mar=c(10, 5, 3, 2.1))
barplot(stats$Estimate, main="",ylab='Effect size', names.arg=stats$Celltype, ylim=c(-0.5,1.5), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=stats$colours)
dev.off()

## Neuronal & non-neuronal ##
pdf(paste0("Barplot_DMP",celltype,"specific","_Peak_n54_logReg_new.pdf"), width=20,height=7)
par(mar=c(10, 5, 3, 2.1))
barplot(stats$Estimate, main="",ylab='Effect size', names.arg=stats$Celltype, ylim=c(-2,3),las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=stats$colours)
dev.off()


# volcano plot: effect size vs p-value
stats$result <- factor(stats$result, levels=c("Not significant","Over-enriched","Under-enriched"))
stats$colours <- factor(stats$colours, levels=c("lightgrey","firebrick","steelblue3")) # to match order of result levels

## bulk ##
pdf("Volcano_DMP_Peak_n54_logReg_new.pdf", width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj, na.rm=T),(max(stats$logP.adj, na.rm=T)+10)), xlim=c(-1,2), ylab='-log10(p.adj)', xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes','Excitatory'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes','Excitatory'),'logP.adj'], labels=c('Astrocytes','Excitatory\nneurons'), pos=4, cex=0.8, col='black') 
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()

## Neuronal ##
pdf(paste0("Volcano_DMP",celltype,"specific","_Peak_n54_logReg_new.pdf"), width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj, na.rm=T),(max(stats$logP.adj, na.rm=T)+5)), ylab='-log10(p.adj)', xlim=c(-2,4), xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes','Excitatory','Stromal'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes','Excitatory','Stromal'),'logP.adj'], labels=c('Astrocytes','Excitatory\nneurons','Stromal cells'), pos=c(4,4,3), cex=0.8, col='black')
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()

## non-neuronal ##
pdf(paste0("Volcano_DMP",celltype,"specific","_Peak_n54_logReg_new.pdf"), width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj, na.rm=T),(max(stats$logP.adj, na.rm=T)+5)), xlim=c(-2,4), ylab='-log10(p.adj)', xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes'),'logP.adj'], labels=c('Astrocytes'), pos=4, cex=0.8, col='black')
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()