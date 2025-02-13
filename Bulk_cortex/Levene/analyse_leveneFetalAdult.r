

## Compare Levene's test stats (early- and mid-fetal samples vs adult samples) ##


library(data.table)
library(scales)
library(ggplot2)
'%ni%' <- Negate('%in%')


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


#2. Load results ================================================================================================================
lev <- data.frame(readRDS(paste0(AnalysisPath, "BulkLifecourse_levene.rds")))
lev <- lev[which(rownames(lev) %in% rownames(betas)),]


#3. Compare stats ===============================================================================================================
# T-test of fetal and adult mean variances
t.test(lev$Var.fetal, lev$Var.adult)$p.value


#4. Plot ========================================================================================================================
df <- data.frame(Var=c(lev[,'Var.fetal'],lev[,'Var.adult']), Group=c(rep('Fetal',nrow(lev)),rep('Adult',nrow(lev))))
df$Var <- df$Var*(100**2) # variance calculated with respect to DNAm measured as a proportion (0-1). Converting to variance of DNAm as a percentage (0-100)
pdf(paste0(plotPath, "BulkLifecourse_levene_top10000DMPs_densityVariance_percDNAm.pdf"), width=12, height=8)
ggplot(df, aes(x=Var, fill=Group)) +
	geom_density(alpha=0.3)+
	ylab("Density")+
	xlab("Variance")+
	theme_minimal()+
	theme(axis.text=element_text(size=19), 
	axis.title=element_text(size=22), 
	plot.title = element_text(size=23), 
	legend.title = element_text(size=18), 
	legend.text = element_text(size=16))
dev.off()
