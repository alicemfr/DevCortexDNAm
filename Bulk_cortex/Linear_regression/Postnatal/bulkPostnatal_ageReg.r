

## Run Age linear regression on postnatal bulk samples at the bulk fetal dDMPs ##


library(data.table)


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"EPICBrainLifecourse.rdat"))

betas <- epic.betas
pheno <- epic.pheno


#2. Load results ================================================================================================================
res.fetal <- data.frame(readRDS(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))


#3. Filter samples and sites ====================================================================================================
postnat <- which(pheno$AgeCat=='Postnatal')
pheno <- pheno[postnat,]
betas <- betas[,postnat]
pheno <- droplevels(pheno)
identical(rownames(pheno),colnames(betas))

# load cell-type proportion predictions
cellpred <- read.csv("postnatalCellProps.csv", row.names=1)
pheno <- cbind(pheno, cellpred)

# limit postnatal betas to fetal dDMPs
sig <- res.fetal[which(res.fetal$P.Age<9e-8),] #50913
betas.sig <- betas[which(rownames(betas) %in% rownames(sig)),]
identical(rownames(pheno),colnames(betas.sig))


#4. Define regression function ==================================================================================================
# include 3 of the 4 cell-types to avoid collinearity
agereg.function <- function(row, pheno){
	lmod <- lm(row ~ Age.years + Sex + Cohort + Sox10 + IRF8 + TN, data=pheno)
	age_stats <- coef(summary(lmod))['Age.years',c('Estimate','Std. Error','Pr(>|t|)')]
	sex_stats <- coef(summary(lmod))['SexM',c('Estimate','Std. Error','Pr(>|t|)')]
	stats <- c(age_stats, sex_stats)
	return(stats)
}


#5. Run EWAS ====================================================================================================================
mat <- matrix(NA, nrow=nrow(betas.sig), ncol=6)
prg <- txtProgressBar()
for(i in 1:nrow(betas.sig)){
	setTxtProgressBar(prg, i/nrow(betas.sig))
	mat[i,] <- agereg.function(as.numeric(betas.sig[i,]), pheno)
}

rownames(mat) <- rownames(betas.sig)
colnames(mat) <- c("Age.ES","Age.SE","Age.P", "Sex.ES","Sex.SE","Sex.P")

saveRDS(mat, file=paste0(resultsPath, "postnatalAgeReg_bulkFetalDMPs.rds"))