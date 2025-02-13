
## Runs FANS EWASs ##

args <- commandArgs(trailingOnly=TRUE)

mod <- args[1]
age.group <- args[2]
output.file <- args[3]


#1. Load data ===================================================================================================================
print("Loading data")

print(paste0("Loading betas: ", MethylationPath, "FANSlifecourse_N_NN.rdat"))
load(paste0(MethylationPath, "FANSlifecourse_N_NN.rdat"))
print(paste0("File will be saved to ...", AnalysisPath, output.file))

columns <- c('Estimate.','SE.','P.')
output.colnames <- c("Anova.P",paste0(columns, mod))
if(mod=="AgeCell"){ output.colnames <- c('Anova.P', paste0(columns,'Age:CellType'), paste0(columns,'Age'), paste0(columns,'CellType')) }
if(mod=="AgeSex"){ output.colnames <- c('Anova.P', paste0(columns,'Age:Sex'), paste0(columns,'Age'), paste0(columns,'CellType')) }
if(mod=="SexCell"){ output.colnames <- c('Anova.P', paste0(columns,'Sex:CellType'), paste0(columns,'Sex'), paste0(columns,'CellType')) }
if(mod=='AgeCellSpecific'){ output.colnames <- c("ES.Age","P.Age","ES.Sex","P.Sex") }

print("Output colnames will be ...")
print(output.colnames)


#2. Separate datasets ===========================================================================================================
print(paste0("Separating into ", age.group, " subset."))

if(age.group=='Fetal'){
	samples <- which(SampleSheet$Phenotype==age.group & SampleSheet$AgeBin!='Late') # only include early and mid fetal samples
}else{
	samples <- which(SampleSheet$Phenotype==age.group)
}
betas.samples <- betas[,samples]                 
pheno.samples <- SampleSheet[samples,]
pheno.samples <- droplevels(pheno.samples)
identical(rownames(pheno.samples),colnames(betas.samples))

if(age.group=='Adult'){
	# create dummy variables for 2 (Sox10, IRF8) out of 3 (Sox10, IRF8, DoubleNeg) non-neuronal cell populations
	pheno.samples$Sox10 <- rep(0,nrow(pheno.samples))
	pheno.samples$Sox10[pheno.samples$Cell_Type=='Sox10+'] <- 1
	pheno.samples$IRF8 <- rep(0,nrow(pheno.samples))
	pheno.samples$IRF8[pheno.samples$Cell_Type=='IRF8+'] <- 1
}

#3. Attach packages =============================================================================================================
library(lme4)
library("lmerTest")
library(stats)
library(doParallel)

cl<-makeCluster(10)
registerDoParallel(cl)


#4. Load functions ==============================================================================================================
print("Loading functions")
source(paste0(ScriptsPath, "FANS_EWAS_functions.r"))

print(paste0("Model to be used ...", mod))

# If running the within cell-type model, separate betas into neuronal and non-neuronal
if(mod=='AgeCellSpecific'){
	neur <- which(pheno.samples$NewCellType=='Neuronal')
	nonneur <- which(pheno.samples$NewCellType=='Non-neuronal')

	betas.neur <- betas.samples[,neur]
	pheno.neur <- pheno.samples[neur,]
	pheno.neur <- droplevels(pheno.neur)
	if(identical(rownames(pheno.neur),colnames(betas.neur))!=TRUE){stop("rownames(pheno.neur) != colnames(betas.neur)")}
	
	betas.nonneur <- betas.samples[,nonneur]
	pheno.nonneur <- pheno.samples[nonneur,]
	pheno.nonneur <- droplevels(pheno.nonneur)
	if(identical(rownames(pheno.nonneur),colnames(betas.nonneur))!=TRUE){stop("rownames(pheno.nonneur) != colnames(betas.nonneur)")}
}

ANOVA.model <- get(mod)


#5. Run EWAS ====================================================================================================================
print(paste("Starting EWAS using", mod, "function"))
print(paste("Starting EWAS on", nrow(pheno.samples), "samples."))

clusterExport(cl, list(mod, "lmer"))

if(mod=='AgeCellSpecific'){

	# neuronal
	print("Running EWAS 1")
	res.neur <- foreach(i=1:nrow(betas.neur), .combine=rbind, .packages=c('lme4','stats'), .verbose=F) %dopar%{
				ANOVA.model(row=as.numeric(betas.neur[i,]), pheno.neur)}
	print("Finished EWAS 1")
	
	# non-neuronal
	print("Running EWAS 2")
	res.nonneur <- foreach(i=1:nrow(betas.nonneur), .combine=rbind, .packages=c('lme4','stats'), .verbose=F) %dopar%{
				ANOVA.model(row=as.numeric(betas.nonneur[i,]), pheno.nonneur)}
	print("Finished EWAS 2")
	
}else{
	print("Starting EWAS")
	res <- foreach(i=1:nrow(betas.samples), .combine=rbind, .packages=c('lme4','stats'), .verbose=F) %dopar%{
				ANOVA.model(row=as.numeric(betas.samples[i,]), pheno.samples)}
	print("Finished EWAS.")
}

print("Finalising results table.")
if(mod=='AgeCellSpecific'){
	rownames(res.neur)<-rownames(betas); colnames(res.neur) <- paste0('Neuronal.',output.colnames)
	rownames(res.nonneur)<-rownames(betas); colnames(res.nonneur) <- paste0('Non.neuronal.',output.colnames)
	res <- cbind(res.neur, res.nonneur)
}else{
	rownames(res)<-rownames(betas); colnames(res) <- output.colnames
}


#6. Save results ================================================================================================================
print(paste0("Saving results to ... ", AnalysisPath, output.file))
saveRDS(res, file=paste0(AnalysisPath, output.file))