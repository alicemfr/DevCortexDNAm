

## EWAS functions for FANS dataset ##


# Cell-type EWAS
Cell <- function(row, pheno){

	if(age.group=='Fetal'){
		full.model <- lmer(row ~ NewCellType + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
		null.model <- lmer(row ~ Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
	}
	else{
		full.model <- lmer(row ~ NewCellType + Sox10 + IRF8 + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
		null.model <- lmer(row ~ Sox10 + IRF8 + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
	}
	
	ANOVA.P <- anova(full.model,null.model)["full.model","Pr(>Chisq)"]
	stats <- coef(summary(full.model))['NewCellTypeNon-neuronal',c('Estimate','Std. Error','Pr(>|t|)')]
	stats.combo <- c(ANOVA.P, stats) #combine stats from ANOVA and lmer

	return(stats.combo)
}


# Age*Cell EWAS
AgeCell <- function(row, pheno){
	if(age.group=='Fetal'){
		full.model <- lmer(row ~ Age*NewCellType + NewCellType + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
		null.model <- lmer(row ~ NewCellType + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
	}
	else{
		full.model <- lmer(row ~ Age*NewCellType + NewCellType + Sox10 + IRF8 + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
		null.model <- lmer(row ~ NewCellType + Sox10 + IRF8 + Age + Sex + Plate + (1|Individual_ID), data=pheno, REML = FALSE)
	}
	
	ANOVA.P <- anova(full.model,null.model)["full.model","Pr(>Chisq)"]
	stats <- c(coef(summary(full.model))['Age:NewCellTypeNon-neuronal',c('Estimate','Std. Error','Pr(>|t|)')], coef(summary(full.model))['Age',c('Estimate','Std. Error','Pr(>|t|)')], coef(summary(full.model))['NewCellTypeNon-neuronal',c('Estimate','Std. Error','Pr(>|t|)')])
	stats.combo <- c(ANOVA.P, stats) #combine stats from ANOVA and lmer

	return(stats.combo)
}


# Age EWAS - separate cell types
AgeCellSpecific <- function(row, pheno){
	if(age.group=='Fetal'){
		full.model <- lm(row ~ Age + Sex + Plate, data=pheno)
	}
	else{
		full.model <- lm(row ~ Age + Sex + Plate, data=pheno)
	}
	
	Age <- coef(summary(full.model))['Age',c('Estimate','Pr(>|t|)')]
	Sex <- coef(summary(full.model))['SexM',c('Estimate','Pr(>|t|)')]
	stats <- c(Age, Sex)
	return(stats)
}