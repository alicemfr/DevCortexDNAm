

## Test assocation of global DNA methylation with Age in bulk fetal cortex ##


library(data.table)

epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F) # Illumina EPIC manifest

#1. Load data ===================================================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))


#2. Pre-processing ==============================================================================================================
pheno$PCW <- as.numeric(as.character(pheno$PCW))
pheno$Sex <- as.factor(pheno$Sex)
pheno$Plate <- as.factor(pheno$Plate)


#3. Limit to autosomal probes ===================================================================================================
probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicMan <- epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
XY <- rownames(probes)[probes$CHR %in% c('X','Y')]
betas.noXY <- betas[-which(rownames(betas) %in% XY),] #790180


#4. Load bulk fetal Age EWAS results ============================================================================================
res <- readRDS(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annot_sig.rds")
betas.noXY.nonsig <- betas.noXY[-which(rownames(betas.noXY) %in% rownames(res)),] # remove significantly linear probes. 739839


#5. Calculate average DNAm per sample ===========================================================================================
# *100 to convert to % DNAm
avBetas1 <- as.numeric(colMeans(betas.noXY*100)) # all autosomal sites
avBetas2 <- as.numeric(colMeans(betas.noXY.nonsig*100)) # autosomal sites excluding dDMPs


#6. Regress average DNAm on age =================================================================================================
mod1 <- lm(avBetas1 ~ PCW + Sex + Plate, data=pheno) # ES = -0.0194, p = 3.33e-10
mod2 <- lm(avBetas2 ~ PCW + Sex + Plate, data=pheno) # ES = 0.007936, p = 0.00974


#7. Plot ========================================================================================================================
pdf(paste0(plotPath, "globalDNAm.pdf"))

par(mfrow=c(2,2))
par(mar=c(5.1, 4.5, 4.1, 2.1))

# all autosomal sites
plot(pheno$PCW, avBetas1, pch=16, xlab='Age (pcw)', ylab='Mean DNA methylation (%)', main='', cex.lab=1.2)
abline(mod1)
mtext("                                                             All autosomal DNA methylation sites", side=3, line=2, cex=1.2)
# as above but with fixed y-axis
plot(pheno$PCW, avBetas1, pch=16, xlab='Age (pcw)', ylab='Mean DNA methylation (%)', ylim=c(0,100), cex.lab=1.2)
abline(mod1)
# autosomal sites excluding dDMPs
plot(pheno$PCW, avBetas2, pch=16, xlab='Age (pcw)', ylab='Mean DNA methylation (%)', cex.lab=1.2)
abline(mod2)
mtext("                                                            Excluding dDMPs", side=3, line=2, cex=1.2)
# as above but with fixed y-axis
plot(pheno$PCW, avBetas2, pch=16, xlab='Age (pcw)', ylab='Mean DNA methylation (%)', ylim=c(0,100), cex.lab=1.2)
abline(mod2)

dev.off()