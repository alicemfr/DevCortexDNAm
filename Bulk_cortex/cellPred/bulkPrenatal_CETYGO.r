

## Deconvolute early/mid-fetal bulk cortex samples ##

library(CETYGO)
library(ggplot2)
library(reshape2)


#1. Load testing data ===========================================================================================================

# bulk fetal
load(paste0(PathToBetas,"fetalBulk_EX3_23pcw_n91.rdat"))
betas.bulk <- betas
pheno.bulk <- pheno
identical(rownames(pheno.bulk),colnames(betas.bulk))
pheno.bulk$Age <- pheno.bulk$PCW

# FANS lifecourse
load(paste0(PathToBetas,"FANSlifecourse_N_NN.rdat"))
# extract fetal and child samples > 20pcw
betas.fans.satb2 <- betas[,SampleSheet$Phenotype %in% c('Fetal','Child')]
pheno.fans.satb2 <- SampleSheet[SampleSheet$Phenotype %in% c('Fetal','Child'),]
pheno.fans.satb2 <- droplevels(pheno.fans.satb2)
betas.fans.satb2 <- betas.fans.satb2[,-which(pheno.fans.satb2$Phenotype=='Fetal' & pheno.fans.satb2$Age<20)]
pheno.fans.satb2 <- pheno.fans.satb2[-which(pheno.fans.satb2$Phenotype=='Fetal' & pheno.fans.satb2$Age<20),]
identical(rownames(pheno.fans.satb2),colnames(betas.fans.satb2))
pheno.fans.satb2$NewCellType <- as.character(pheno.fans.satb2$NewCellType)
# extract adult samples
betas.fans.neun <- betas[,SampleSheet$Phenotype=='Adult']
pheno.fans.neun <- SampleSheet[SampleSheet$Phenotype=='Adult',]
identical(rownames(pheno.fans.neun),colnames(betas.fans.neun))
pheno.fans.neun$NewCellType <- as.character(pheno.fans.neun$NewCellType)


#2. Load training data ==========================================================================================================
load("GSE38214_betas_pheno.rdat") # betas pheno
# betas and pheno for 'GSE38214' donwloaded from GEO

betas <- betas[,pheno$Cell_Type %in% c('Embryonic stem (ES) cells H9','ES-derived neural precursor cells')]
pheno <- pheno[pheno$Cell_Type %in% c('Embryonic stem (ES) cells H9','ES-derived neural precursor cells'),]
identical(rownames(pheno),colnames(betas))
pheno$CT <- pheno$Cell_Type
pheno$CT[pheno$CT=='Embryonic stem (ES) cells H9'] <- 'ES'
pheno$CT[pheno$CT=='ES-derived neural precursor cells'] <- 'ES-NPC'

# add late-fetal and child SATB2+ to training data
pheno.fans.satb2 <- pheno.fans.satb2[,c('Sex','NewCellType')]
pheno.fans.satb2 <- pheno.fans.satb2[pheno.fans.satb2$NewCellType=='Neuronal',] # limit to SATB2+ (neuronal)
pheno.fans.satb2$CT <- as.character(pheno.fans.satb2$NewCellType)
pheno.fans.satb2$CT[pheno.fans.satb2$CT=='Neuronal'] <- 'SATB2pos' #SATB2+
colnames(pheno.fans.satb2) <- c('Sex','Cell_Type','CT')
betas.fans.satb2 <- betas.fans.satb2[,colnames(betas.fans.satb2) %in% rownames(pheno.fans.satb2)]

# add adult NeuN+ to training data
pheno.fans.neun <- pheno.fans.neun[,c('Sex','NewCellType')]
pheno.fans.neun$CT <- as.character(pheno.fans.neun$NewCellType)
pheno.fans.neun$CT[pheno.fans.neun$CT=='Neuronal'] <- 'NeuNpos' # NeuN+
pheno.fans.neun$CT[pheno.fans.neun$CT=='Non-neuronal'] <- 'NeuNneg' # NeuN-
colnames(pheno.fans.neun) <- c('Sex','Cell_Type','CT')

pheno <- rbind(pheno, pheno.fans.satb2, pheno.fans.neun)
betas.fans.satb2.2 <- betas.fans.satb2[intersect(rownames(betas.fans.satb2),rownames(betas)),]
betas <- betas[intersect(rownames(betas),rownames(betas.fans.satb2.2)),]
betas.fans.neun.2 <- betas.fans.neun[intersect(rownames(betas.fans.neun),rownames(betas)),]
betas <- betas[intersect(rownames(betas),rownames(betas.fans.neun.2)),]
all(rownames(betas) %in% rownames(betas.fans.satb2.2))
all(rownames(betas) %in% rownames(betas.fans.neun.2))
betas <- cbind(betas, betas.fans.satb2.2, betas.fans.neun.2)
identical(rownames(pheno),colnames(betas))
betas <- as.matrix(betas)


#3. Create custom model =========================================================================================================
customModel <- pickCompProbesMatrix(rawbetas = betas,
                                    cellTypes = unique(pheno$CT),
                                    cellInd = pheno$CT,
                                    numProbes = 100,
                                    probeSelect = "any")
									
									
#4. Deconvolute testing data ====================================================================================================

rInd <- rownames(betas.bulk)[rownames(betas.bulk) %in% rownames(customModel$coefEsts)]
predProp <- projectCellTypeWithError(betas.bulk, customModel$coefEsts[rInd,])
predProp <- data.frame(predProp)


#5. Plot ========================================================================================================================
df.plot <- melt(predProp)
df.plot <- df.plot[-which(df.plot$variable=='nCGmissing'),]
colnames(df.plot)[colnames(df.plot)=='variable'] <- 'Cell-type'

age <- pheno.bulk$Age
df <- cbind(df.plot, Age=rep(age, length(unique(df.plot$`Cell-type`))))

# Scatter: predicted proportion vs sample age
df.plot <- df[-which(df$`Cell-type`=='CETYGO'),]
ggplot(df.plot, aes(x=Age, y=value, group=`Cell-type`, colour=`Cell-type`))+
	geom_smooth(aes(fill=`Cell-type`))+
	ylab("Predicted proportion")+
	xlab("Age (pcw)")+
	theme_minimal()+
	theme(axis.text=element_text(size=20), axis.title=element_text(size=20))