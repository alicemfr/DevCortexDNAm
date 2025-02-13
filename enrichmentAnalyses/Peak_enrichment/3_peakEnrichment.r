

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


#2. Logistic regression =========================================================================================================
# test association between a site being a DMP and being in a peak, whilst controlling for all other cell-type peaks
f <- as.formula(paste("DMP", "~", paste(cells, collapse='+')))
model <- glm(f, data=res)
sum_mod <- summary(model)$coefficients
sum_mod <- sum_mod[-which(rownames(sum_mod)=='(Intercept)'),]


#3. Plot ========================================================================================================================
# prepare df for plotting
stats <- data.frame(Celltype=cells, Estimate=sum_mod[,'Estimate'], P=sum_mod[,'Pr(>|t|)'])
stats$P.adj <- stats$P*nrow(stats)
stats$P.adj[which(stats$P.adj>1)] <- 1 # p>1 are not meaningful, cap at 1
stats$logP.adj <- -log10(stats$P.adj)
stats$OR <- exp(stats$Estimate)        # calc odds ratio from estimate
stats <- stats[,c('Celltype','Estimate','OR','P','P.adj','logP.adj')]
write.csv(stats, file=outFile, row.names=F)


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


# Barplot: effect size vs cell-type

# bulk
pdf(paste0("Barplot_DMP",celltype,"specific","_Peak_n54_logReg.pdf"), width=20,height=7)
par(mar=c(10, 5, 3, 2.1))
barplot(stats$Estimate, main="",ylab='Effect size', names.arg=stats$Celltype, ylim=c(min(stats$Estimate)-0.05,(max(stats$Estimate)+0.05)), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=stats$colours)
dev.off()

# neuronal
pdf(paste0("Barplot_DMP",celltype,"specific","_Peak_n54_logReg.pdf"), width=20,height=7)
par(mar=c(10, 5, 3, 2.1))
barplot(stats$Estimate, main="",ylab='Effect size', names.arg=stats$Celltype, ylim=c(-0.005,0.05), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=stats$colours)
dev.off()

# non-neuronal
pdf(paste0("Barplot_DMP",celltype,"specific","_Peak_n54_logReg.pdf"), width=20,height=7)
par(mar=c(10, 5, 3, 2.1))
barplot(stats$Estimate, main="",ylab='Effect size', names.arg=stats$Celltype, ylim=c(-0.005,0.05), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2, col=stats$colours)
dev.off()


# Volcano plot: effect size vs p-value
stats$result <- factor(stats$result, levels=c("Not significant","Over-enriched","Under-enriched"))
stats$colours <- factor(stats$colours, levels=c("lightgrey","firebrick","steelblue3")) # to match order of result levels

# bulk
pdf("Volcano_DMP_Peak_n54_logReg.pdf", width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj),(max(stats$logP.adj)+10)), xlim=c(-0.1,0.23), ylab='-log10(p.adj)', xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes','Excitatory'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes','Excitatory'),'logP.adj'], labels=c('Astrocytes','Excitatory\nneurons'), pos=4, cex=0.8, col='black') 
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()

# neuronal
pdf(paste0("Volcano_DMP",celltype,"specific","_Peak_n54_logReg.pdf"), width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj),(max(stats$logP.adj)+10)), xlim=c(-0.02,0.07), ylab='-log10(p.adj)', xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes','Excitatory','Stromal'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes','Excitatory','Stromal'),'logP.adj'], labels=c('Astrocytes','Excitatory\nneurons','Stromal cells'), pos=c(4,4,3), cex=0.8, col='black')
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()


# non-neuronal
pdf(paste0("Volcano_DMP",celltype,"specific","_Peak_n54_logReg.pdf"), width=5, height=5)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(stats$Estimate, stats$logP.adj, pch=c(1,16)[stats$Sig], col=levels(stats$colours)[stats$result], ylim=c(min(stats$logP.adj),(max(stats$logP.adj)+5)), xlim=c(-0.01, 0.032), ylab='-log10(p.adj)', xlab='Effect size', cex.lab=1.5, cex.axis=1.2, cex=c(0.5,1.2)[stats$Sig])
text(x=stats[stats$Celltype %in% c('Astrocytes','Schwann'),'Estimate'], y=stats[stats$Celltype %in% c('Astrocytes','Schwann'),'logP.adj'], labels=c('Astrocytes','Schwann cells'), pos=4, cex=0.8, col='black') 
abline(h=-log10(5e-2), lty=3)
abline(v=0, lty=3)
dev.off()