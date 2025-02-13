

## Comparison of Age linear regression statistics in an independent Illumina EPIC 450K fetal cortex cohort ##

# Fetal 450K cohort: Spiers et al. (2015). DOI:10.1101/gr.180273.114


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


#1. Load results =============================================================================================
res450 # Spiers Age EWAS results table
resEPIC <- data.frame(readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))) # Age EWAS results from current study


#2. Combine results ==========================================================================================
res450 <- res450[which(rownames(res450) %in% rownames(resEPIC)),]
resEPIC <- resEPIC[which(rownames(resEPIC) %in% rownames(res450)),]
identical(rownames(resEPIC), rownames(res450))

# 450K Age results are in days post conception. EPIC results are in weeks post conception.
# Multiply 450K Age results effect sizes by 7 to convert from days to weeks.

res <- data.frame(Beta.Age.450=(res450$Age.estimate*7), Beta.Age.EPIC=resEPIC$Beta.Age, P.Age.450=res450$Age.Pvalue, P.Age.EPIC=resEPIC$P.Age)
rownames(res) <- rownames(resEPIC)
res.age.sig.EPIC <- res[which(res$P.Age.EPIC<9e-8),] #20502 - dDMPs in current study that we have stats for in Spiers et al. results table


#3. Effect size scatter plot =================================================================================
ES.plot <- function(res.df, col1, col2, xlim, ylim, xlab, ylab, cor.pos.x, cor.pos.y, fctr=100){
	X <- (res.df[,col1]*100)
	Y <- (res.df[,col2]*100)
	cor.xy <- signif(cor.test(X,Y)$estimate, digits=3)
	cor.p <- signif(cor.test(X,Y)$p.value, digits=3)
	print(cor.xy)
	print(cor.p)
	# very small p-values appear as 0; calc smallest possible value to quote
	print(paste0("Cor.p < 1e-320"[cor.p<1e-320]))
	print(paste0("Cor.p < 2.2e-16"[cor.p<2.2e-16]))

	ggplot(res.df, aes(x=X, y=Y)) +
		geom_bin2d(bins = nrow(res.df)/fctr) +
		scale_fill_continuous(type = "viridis") +
		ggtitle('')+
		xlab(xlab)+
		ylab(ylab)+
		ylim(ylim*100)+xlim(xlim*100)+
		geom_abline(intercept=0, linetype='dashed') +
		geom_vline(xintercept=0, size=1) +
		geom_hline(yintercept=0, size=1) +
		theme_minimal()+
		theme(axis.text=element_text(size=16), axis.title=element_text(size=22))+
		theme(plot.title = element_text(size=23))
}


pdf(paste0(plotPath,"EffectSizeComparison_450vsEPIC_sigEPIC.pdf"), width=7, height=7)
ES.plot(res.df=res.age.sig.EPIC, col1='Beta.Age.EPIC', col2='Beta.Age.450', ylab='Spiers et al.\nEffect size (%)', xlab='Current study\nEffect size (%)', xlim=c(-0.06,0.06), ylim=c(-0.06, 0.06), cor.pos.x=-0.054, cor.pos.y=0.051)
dev.off()




