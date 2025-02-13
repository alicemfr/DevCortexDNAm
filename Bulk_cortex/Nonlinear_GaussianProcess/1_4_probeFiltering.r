

## Filter non-variable probes and probes with constant DNAm >0.9 or <0.1 ##


#1. Load data ===================================================================================================================
load(paste0(PathToBetas,"fetalBulk_EX3_23pcw_n91.rdat"))                                            #pheno file and original betas
betas <- readRDS(paste0(outputPath,"FetalBrainBetasEX3_23pcw_LOOZ_newmethod_ethsnp.rdat"))          #filtered betas
identical(rownames(pheno), colnames(betas))                                                         #check order of samples in pheno file and filtered betas match


#2. Resave pheno file ===========================================================================================================
pheno <- pheno[,c('Individual_ID','Sample_ID','PCW','Sex')]
pheno$Sample <- rownames(pheno)
# GPMethylation requires 'Age','Sex','Sample', where Sample is the ID that matches colnames of betas
colnames(pheno) <- c('Individual_ID','Sample_ID','Age','Sex','Sample')
write.csv(pheno, file=paste0(outputPath,"FetalBrain_Pheno_23pcw_EX3.csv"))


#3. Load non-variable probe info from '1_3_bulkFetal_nonVarProbes.r' ===================================================================
# non-variable defined as a probe with DNAm change of <5% in middle 80% of samples
nonvar <- readRDS(paste0(outputPath,"nonVarprobes.rds"))                                            #252,504 non-variable probes
betas  <- betas[-which(rownames(betas) %in% nonvar[,'Probe']),]                                     #543,800 variable probes remaining
 

#4. Remove probes where all samples have mean DNAm - 3sd >0.9 or mean DNAm + 3sd <0.1 ===========================================
sds <- apply(betas, 1, sd, na.rm=T)
mns <- apply(betas, 1, mean, na.rm=T)
stats <- data.frame(Probe=rownames(betas), Mean=mns, Std=sds)

# First filter: find which probes have mean DNAm >0.9 or <0.1
mn0.9 <- stats[which(stats$Mean>0.9),] # 1379 probes 
mn0.1 <- stats[which(stats$Mean<0.1),] # 9796 probes 

# Second filter: expand range to look for probes where all samples within 3sd from mean are consistently >0.9 or <0.1
# Reasoning: A probe could have a mean DNAm of 0.9 but this could be due to a nonlinear trend oscillating about a mean of 0.9, 
# which would be biologically interesting so we'd want to keep this probe. Filtering on just mean would remove this probe. The
# mean +/- 3sd approach requires 95% of samples to be >0.9 or <0.1, the regions where technical variation is likely to be driving the effects.

MeanStd0.9 <- c()
for(i in 1:nrow(mn0.9)){
	statement = mn0.9$Mean[i] - 3*mn0.9$Std[i] > 0.9
	MeanStd0.9 <- c(MeanStd0.9, statement)
}

MeanStd0.1 <- c()
for(i in 1:nrow(mn0.1)){
	statement = mn0.1$Mean[i] + 3*mn0.1$Std[i] < 0.1
	MeanStd0.1 <- c(MeanStd0.1, statement)
}

mn0.9Probes <- rownames(mn0.9)[which(MeanStd0.9==TRUE)] #3 probes
mn0.1Probes <- rownames(mn0.1)[which(MeanStd0.1==TRUE)] #1 probes

Constant0.90.1probes <- c(mn0.9Probes, mn0.1Probes)
stats.constant <- stats[match(Constant0.90.1probes, rownames(stats)),]

# Plot a few out to check it's doing what we expect
# par(mfrow=c(3,3))
# for(i in 1:9){	probe <- mn0.9Probes[i]; 	plot(pheno$PCW, betas[probe,]) }

betas.noConst <- betas[-which(rownames(betas) %in% rownames(stats.constant)),] # 543,796 probes remaining


#5. Save betas ==================================================================================================================
# GPMethylation requires betas object to be named betas.dasen, and the betas file to have extension .RData
betas.dasen <- betas.noConst
save(betas.dasen, file=paste0(outputPath,"FetalBrain_Betas_noHiLowConst.RData"))