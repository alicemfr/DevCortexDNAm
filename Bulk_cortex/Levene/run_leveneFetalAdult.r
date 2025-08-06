

## Run Levene's test per DNAm site, comparing variance of early- and mid-fetal samples vs adult samples ##


library(data.table)
library(car)
library(pbapply)
'%ni%' <- Negate('%in%')


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"EPICBrainLifecourse.rdat")) # epic.betas epic.pheno
betas <- epic.betas
pheno <- epic.pheno


#2. Limit samples ===============================================================================================================
samples <- which(pheno$Phenotype %in% c('Fetal','Adult') & pheno$AgeBin!='Late')
pheno <- pheno[samples,]
betas <- betas[,samples]
identical(rownames(pheno),colnames(betas))

pheno$Phenotype <- as.factor(pheno$Phenotype)


#3. Run Levene's test ===========================================================================================================
levtest <- function(x){
	x <- as.numeric(x)                              # make DNAm numeric
	fetal <- which(pheno$Phenotype=='Fetal')        # fetal samples
	adult <- which(pheno$Phenotype=='Adult')        # adult samples
	fetal.var <- var(x[fetal])                      # variance of DNAm for fetal samples
	adult.var <- var(x[adult])                      # variance of DNAm for adult samples
	mod <- leveneTest(x ~ Phenotype, data=pheno)    # run levene's test comparing DNAm variance between fetal and adult
	p <- mod['group','Pr(>F)']                      # extract p-value for levene's test
	return( c(fetal.var, adult.var, p) )            # return fetal variance, adult variance and levene's test p-value
}

lev <- t(pbapply(betas, 1, levtest))                # run levtest per row of DNAm. Transpose output so rows = sites, cols = stats
colnames(lev) <- c('Var.fetal','Var.adult','Levene.P')

saveRDS(lev, file=paste0(AnalysisPath, "BulkLifecourse_levene.rds"))