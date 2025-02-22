

## Apply epigenetic clock functions to late fetal, child and adult samples aged between 26 pcw - 104 years.
# Hovath Clock function (https://www.rdocumentation.org/packages/wateRmelon/versions/1.16.0/topics/agep)
# Cortical Clock function (https://github.com/gemmashireby/CorticalClock)

library(ggplot2)
library(wateRmelon)
library(grid)

PathToBetas <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/2_normalised/"
corticalScript <- "/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/Bulk/3_analysis/Scripts/"
plotPath <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/Figures/"
setwd("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/3_analysis/")

#1. Load data ===================================================================================================================
load(paste0(PathToBetas,"EPICBrainLifecourse.rdat"))
pheno <- epic.pheno
betas <- epic.betas
betas <- betas[,which(pheno$AgeBin=='Late' | pheno$Phenotype %in% c('Child','Adult'))]
pheno <- pheno[which(pheno$AgeBin=='Late' | pheno$Phenotype %in% c('Child','Adult')),] #677
identical(rownames(pheno),colnames(betas))

#2. Run Horvath clock ===========================================================================================================
pheno$PredAge <- agep(betas)
write.csv(pheno,file="predAgePostnatal_Horvath.csv")

#3. Run Cortical clock ==========================================================================================================
dir <- "/gpfs/mrc0/projects/Research_Project-MRC190311/references/CortexClock/" #reference coeffs
source(paste0(corticalScript,"corticalClock_GShireby.r"))

#4. Correlate reported vs predicted ages ========================================================================================
horvath <- read.csv("predAgePostnatal_Horvath.csv", row.names=1)
shireby <- read.csv("CorticalPred.csv", row.names=1)
identical(horvath$Basename, shireby$Basename)
df <- data.frame(ReportedAge=horvath$Age.years, Horvath=horvath$PredAge, Shireby=shireby$brainpred)

RepvsPred <- function(df, predCol, xlab="Reported age (years)", ylab="Predicted age (years)"){
	d <- data.frame(ReportedAge=df$ReportedAge, PredictedAge=df[,predCol])
	cor.age <- signif(cor.test(d$ReportedAge, d$PredictedAge)$estimate, digits=3)
	print(cor.age)
	#grob <- grobTree(textGrob(paste("Cor =",cor.age), x=0.2,  y=0.95, gp=gpar(col="black", fontsize=20)))
	maxVal <- max( max(d$ReportedAge), max(d$PredictedAge) ) #max plot value
	minVal <- min( min(d$ReportedAge), min(d$PredictedAge) ) #min plot value
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
			plot.margin = margin(0.1, 0.5, 0.1, 0.1, "cm")
		)#+
		#annotation_custom(grob)
}

# horvath clock
pdf(paste0(plotPath, "bulkPostnatal_ReportedvsPredicted_HorvathClock2.pdf"), width=7, height=7)
RepvsPred(df, predCol='Horvath')
dev.off()

# cortex clock
pdf(paste0(plotPath, "bulkPostnatal_ReportedvsPredicted_CorticalClock2.pdf"), width=7, height=7)
RepvsPred(df, predCol='Shireby')
dev.off()




