

## Compare mean DNA methylation at the start of early fetal, the end of mid-fetal, and adult ##


library(data.table)
library(pbapply)
library(ggplot2)
library(viridis)
'%ni%' <- Negate('%in%')


avDNAm <- function(x){
	mn <- mean(as.numeric(x))
	return(mn)
}


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"EPICBrainLifecourse.rdat"))

betas <- epic.betas
pheno <- epic.pheno

# limit to top 10,000 fetal dDMPs
res <- readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))
sig <- res[which(res$P.Age<9e-8),]
ord <- sig[order(sig$P.Age),]
betas <- betas[which(rownames(betas) %in% rownames(ord)),] #41518
betas <- betas[match(rownames(ord)[which(rownames(ord) %in% rownames(betas))],rownames(betas)),]
identical(rownames(betas),rownames(ord)[which(rownames(ord) %in% rownames(betas))])
betas <- betas[1:10000,]


#2. Limit samples ===============================================================================================================
samples <- which(pheno$Phenotype %in% c('Fetal','Adult') & pheno$AgeBin!='Late')
pheno <- pheno[samples,]
betas <- betas[,samples]
identical(rownames(pheno),colnames(betas))
  

#3. Average DNAm for last 3 fetal ages ==========================================================================================
fetal <- pheno[which(pheno$Phenotype=='Fetal'),]
ages <- as.numeric(unique(fetal$Age))
ord <- ages[order(ages,decreasing=T)] # order oldest to youngest
top <- ord[1:3]
oldest <- fetal[which(fetal$Age %in% top),]
betas.oldest <- betas[,which(colnames(betas) %in% rownames(oldest))]
identical(rownames(oldest),colnames(betas.oldest))
av.fetal.oldest <- pbapply(betas.oldest, 1, avDNAm)


#4. Average DNAm for first 3 fetal ages =========================================================================================
ord <- ages[order(ages,decreasing=F)] # order youngest to oldest
top <- ord[1:3]
youngest <- fetal[which(fetal$Age %in% top),]
betas.youngest <- betas[,which(colnames(betas) %in% rownames(youngest))]
identical(rownames(youngest),colnames(betas.youngest))
av.fetal.youngest <- pbapply(betas.youngest, 1, avDNAm)


#5. Average DNAm for the adult samples ==========================================================================================
adult <- pheno[which(pheno$Phenotype=='Adult'),]
betas.adult <- betas[,which(colnames(betas) %in% rownames(adult))]
identical(rownames(adult),colnames(betas.adult))
av.adult <- pbapply(betas.adult, 1, avDNAm)


#6. Combine results =============================================================================================================
df <- data.frame(Fetal.youngest.mean=av.fetal.youngest, Fetal.oldest.mean=av.fetal.oldest, Adult.mean=av.adult)
df <- df*100


#7. Plot results ================================================================================================================
# Early-fetal vs Adult
ggplot(df, aes(x=Fetal.youngest.mean, y=Adult.mean) ) +
  geom_bin2d(bins = nrow(df)/100) +
  scale_fill_continuous(type = "viridis") +
  xlab("Average DNAm (%)\nStart of early-fetal") + ylab("Average DNAm (%)\nAdult") +
  ylim(0,100)+ xlim(0,100)+
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))

# Mid-fetal vs Adult
ggplot(df, aes(x=Fetal.oldest.mean, y=Adult.mean) ) +
  geom_bin2d(bins = nrow(df)/100) +
  scale_fill_continuous(type = "viridis") +
  xlab("Average DNAm (%)\nEnd of mid-fetal") + ylab("Average DNAm (%)\nAdult") +
  ylim(0,100)+ xlim(0,100)+
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))

# Early-fetal vs Mid-fetal
ggplot(df, aes(x=Fetal.youngest.mean, y=Fetal.oldest.mean) ) +
  geom_bin2d(bins = nrow(df)/100) +
  scale_fill_continuous(type = "viridis") +
  xlab("Average DNAm (%)\nStart of early-fetal") + ylab("Average DNAm (%)\nEnd of mid-fetal") +
  ylim(0,100)+ xlim(0,100)+
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))