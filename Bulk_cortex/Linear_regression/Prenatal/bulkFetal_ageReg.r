

## Linear regression to test association of DNA methylation with Age and Sex in bulk fetal cortex ##


library(data.table)
library(doParallel)


#1. Load data ===================================================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))


#2. Pre-processing ==============================================================================================================
pheno$PCW <- as.numeric(as.character(pheno$PCW)) # PCW = age in post-conceptual weeks (numeric 6-23)
pheno$Sex <- as.factor(pheno$Sex) # Sex (factor M/F)
pheno$Plate <- as.factor(pheno$Plate) # Plate = experimental array plate number (factor 1/2)


#3. Regression function =========================================================================================================
agereg.function <- function(row, pheno){
  lmod <- lm(row ~ PCW + Sex + Plate, data=pheno)
  age <- summary(lmod)$coefficients['PCW', c(1,2,4)] # Extract Age Estimate, Standard.Error and P.value
  sex <- summary(lmod)$coefficients['SexM', c(1,2,4)] # Extract Sex Estimate, Standard.Error and P.value
  stats <- c(age, sex)
  return(stats)
}


#4. Run Age EWAS ================================================================================================================
print(paste0("Starting regression on ", nrow(betas), " sites."))
cl <- makeCluster(10)
registerDoParallel(cl)
res <- foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
	agereg.function(row=betas[i,], pheno)
}
rownames(res) <- rownames(betas)
colnames(res) <- c("Beta.Age", "SE.Age", "P.Age", "Beta.Sex", "SE.Sex", "P.Sex")


#4. Filter and save results =====================================================================================================

# annotate with Illumina EPIC manifest (MethylationEPIC_v-1-0_B4.csv) with extra annotations for sites in known autism, developmental delay and schizophrenia genes
epicManifest <- fread(paste0(EPICmanPath, "epicMan_genelist_AUT.DDD.SCZ.csv"), data.table=F)
epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c('CHR','MAPINFO','UCSC_RefGene_Name','UCSC_RefGene_Group','Relation_to_UCSC_CpG_Island','GeneList','TargetGene')]
res <- cbind(res, as.data.frame(epicMan))

# remove Y chromosome as DNAm values uninformative for females
res <- res[-which(res$CHR=='Y'),]

# extract unique mentions of genes
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))
saveRDS(res, file=paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))

