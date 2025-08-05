

## Comparison of Age linear regression statistics in an independent Illumina EPIC 27K fetal cortex cohort ##

# Fetal 27K cohort: Numata et al. (2012). DOI:10.1016/j.ajhg.2011.12.020


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


#1. Load results =============================================================================================
res27 # Numata Age EWAS results table
resEPIC <- data.frame(readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))) # Age EWAS results from current study


#2. Combine results ==========================================================================================
res27 <- res27[which(rownames(res27) %in% rownames(resEPIC)),]
resEPIC <- resEPIC[which(rownames(resEPIC) %in% rownames(res27)),]
resEPIC <- resEPIC[match(rownames(res27),rownames(resEPIC)),]
identical(rownames(resEPIC), rownames(res27))

# 27K Age results are in years. EPIC results are in weeks post conception.
# Convert 27K Age results effect sizes from years to weeks.
# 52 weeks in year. Effect size per year is 52 times greater than per week.

res27$Beta.Age <- res27$beta.regression.coefficient/52
res <- data.frame(Beta.Age.27=res27$Beta.Age, Beta.Age.EPIC=resEPIC$Beta.Age, P.Age.27=res27$p.value, P.Age.EPIC=resEPIC$P.Age)
rownames(res) <- rownames(resEPIC)
res.age.sig.EPIC <- res[which(res$P.Age.EPIC<9e-8),] #66 - dDMPs in current study that we have stats for in Numata et al. results table


#3. Effect size scatter plot =================================================================================
ES.plot <- function(res.df, col1, col2, xlim, ylim, xlab, ylab){
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
		geom_point() +
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


pdf(paste0(plotPath,"EffectSizeComparison_27vsEPIC_sigEPIC.pdf"), width=7, height=7)
ES.plot(res.df=res.age.sig.EPIC, col1='Beta.Age.EPIC', col2='Beta.Age.27', ylab='Numata et al.\nEffect size (%)', xlab='Current study\nEffect size (%)', xlim=c(-0.08,0.08), ylim=c(-0.08, 0.08))
dev.off()

pdf(paste0(plotPath,"EffectSizeComparison_27vsEPIC_sig27K.pdf"), width=7, height=7)
ES.plot(res.df=res, col1='Beta.Age.EPIC', col2='Beta.Age.27', ylab='Numata et al.\nEffect size (%)', xlab='Current study\nEffect size (%)', xlim=c(-0.08,0.08), ylim=c(-0.08, 0.08))
dev.off()