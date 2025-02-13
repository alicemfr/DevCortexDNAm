

## Calculate neuron (NeuN), oligodendrocyte (Sox10), microglia (IRF8) and astrocyte-enriched (TN) proportions in bulk postnatal samples ##
# https://github.com/ejh243/CETYGO


library(CETYGO)


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"EPICBrainLifecourse.rdat"))
betas <- epic.betas
pheno <- epic.pheno


#2. Limit to postnatal samples ==================================================================================================
postnat <- which(pheno$AgeCat=='Postnatal')
pheno <- pheno[postnat,]
betas <- betas[,postnat]
pheno <- droplevels(pheno)
identical(rownames(pheno),colnames(betas))


#3. Run CETYGO ==================================================================================================================
predProp <- projectCellTypeWithError(betas, modelBrainCoef[['ANOVA']][[2]])
predProp <- predProp[,!colnames(predProp)=='nCGmissing']
colnames(predProp) <- c('IRF8','NeuN','Sox10','TN','CETYGO')

write.csv(predProp, "postnatalCellProps.csv")