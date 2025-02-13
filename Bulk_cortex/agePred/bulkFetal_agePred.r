

## Apply epigenetic clock to early fetal samples aged between 6 - 23 pcw ##

# Fetal Clock function: https://github.com/LSteg/EpigeneticFetalClock


library(ggplot2)

#1. Load data ===================================================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))


#2. Run Fetal clock =============================================================================================================
source(paste0(scriptsPath, "FetalClockFunction.R"))
pheno$PredAge <- FetalClock(betas)/7 # predict age and convert units from days to weeks post-conception (pcw)
write.csv(pheno[,c('Individual_ID','Sample_ID','Age','PCW','PredAge','Sex')], file=paste0(outputPath, "bulkFetal_predAge.csv"))


#3. Correlate reported vs predicted ages ========================================================================================
df <- data.frame(ReportedAge=pheno$PCW, PredictedAge=pheno$PredAge)

RepvsPred <- function(df, predCol, xlab="Reported age (pcw)", ylab="Predicted age (pcw)"){
	d <- data.frame(ReportedAge=df$ReportedAge, PredictedAge=df[,predCol])
	cor.age <- signif(cor.test(d$ReportedAge, d$PredictedAge)$estimate, digits=3)
    print(cor.age)
	maxVal <- max( max(d$ReportedAge), max(d$PredictedAge) ) # max plot value
	minVal <- min( min(d$ReportedAge), min(d$PredictedAge) ) # min plot value
	ggplot(d, aes(x=ReportedAge, y=PredictedAge))+
		geom_point(size=3)+
		stat_smooth(method = lm, formula = y ~ x,se=T)+
		ggtitle("")+
		xlab(xlab)+
		ylab(ylab)+
		xlim(minVal, maxVal)+
		ylim(minVal, maxVal)+
		theme_minimal()+
		theme(
			axis.text=element_text(size=19),
			axis.title=element_text(size=22),
			plot.title = element_text(size=23),
			plot.margin = margin(0, 0.5, 0, 0, "cm")
		)
}

pdf(paste0(plotPath, "bulkFetal_ReportedvsPredicted_FetalClock.pdf"), width=7, height=7)
RepvsPred(df, predCol='PredictedAge')
dev.off()



