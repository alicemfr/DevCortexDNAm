

## Principal Component Analysis (PCA) on FANS lifecourse DNAm dataset ##

library(ggplot2)
library(data.table) # for fread
library(RColorBrewer)
library(viridis)
library(scales) # for rescale
library(pscl) # for pR2

celltype_cols <- c(plasma(4)[2],viridis(4)[3]) # plasma(4)[2] = #9C179EFF = purple; viridis(4)[3] = #35B779FF = green


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"FANSlifecourse_N_NN.rdat"))
betas <- betas #693964    600
pheno <- SampleSheet

probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
epicMan <- epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
XY <- which(probes$CHR %in% c('X','Y'))
betas.noXY <- betas[-XY,] # limit to autosomal sites (679917)


#2. Organise pheno columns for plotting =========================================================================================
pheno$Sex <- as.factor(pheno$Sex)
pheno$Cell_Type[pheno$Cell_Type=='SATB2 +'] <- 'SATB2+'
pheno$Cell_Type[pheno$Cell_Type=='SATB2 -'] <- 'SATB2-'
pheno <- droplevels(pheno)
pheno$Cell_Type <- factor(pheno$Cell_Type, levels=c('SATB2+','SATB2-','NeuN+','Sox10+','Double-','IRF8+'))


#3. Run PCA & plot ==============================================================================================================
source(paste0(scriptsPath, "pcaFunctions.r"))

# PCA on all samples ------------------------------------------------------------------------------------------------------------
pca.all <- pca.Gene(pheno, betas.noXY, varProbes=TRUE, nVarprobes=10000, scale=TRUE)

# correlation heatmap for PC1-4
tmp <- pheno
tmp$Phenotype <- as.numeric(tmp$Phenotype)
tmp$NewCellType <- as.numeric(tmp$NewCellType)
tmp$Sex <- as.numeric(tmp$Sex)

pdf(paste0(plotPath,"FANS_PCAcorHeatmap_4PCs_AgeCellSex_narrow.pdf"), width=3, height=4)
corPs <- pcaCorPlot(samplesheet=tmp, columns=c('Phenotype','NewCellType','Sex'), nPCs=4, pca.res=pca.all, colLabels=c('Age','Cell-type','Sex'), returnP=T)
dev.off()


# PCA on adult samples ----------------------------------------------------------------------------------------------------------
pheno.adult <- pheno[which(pheno$Phenotype=='Adult'),]
betas.adult <- betas.noXY[,which(pheno$Phenotype=='Adult')]
identical(rownames(pheno.adult),colnames(betas.adult))
pca.adult <- pca.Gene(pheno.adult, betas.adult, varProbes=TRUE, nVarprobes=10000)

# use adult PCA loadings to project fetal samples into same space
pheno.f <- pheno[which(pheno$Phenotype=='Fetal'),]
f.betas.sub <- betas[which(rownames(betas) %in% names(pca.adult$center)),which(pheno$Phenotype=='Fetal')]
f.betas.sub <- f.betas.sub[match(names(pca.adult$center), rownames(f.betas.sub)),]
identical(rownames(f.betas.sub), names(pca.adult$center))
pca.f <- PCAproject(pca.adult, f.betas.sub)

# plot PC1 vs PC2
pdf(paste0(plotPath,"FANS_adult_projF_colAbLabel_pchPheno_PC1_PC2.pdf"), width=7, height=7)
par(mar=c(5.1, 5, 4.1, 2.1))
ab_colours <- c('#35B779FF','#9C179EFF','#247c52','#ABE2BB','#79C478','#c378c4') #order: Double-, NeuN+, Sox10+, IRF8+, SATB2-, SATB2+
pca.adult.fc <- plotPCA(pca.res=pca.adult, pheno=pheno.adult, title='',colourBy='Cell_Type', pcA=1,pcB=2, legend=T, legendPos='bottomright', projection=TRUE, projectData=pca.f, projectPheno=pheno.f, colourBy.project = 'Cell_Type', pointBy='Phenotype', legendOrder=c('Fetal','Adult','NeuN+','SATB2+','SATB2-','Sox10+','Double-','IRF8+'), cex.legend=1, cex.lab=1.5, cex.main=1.5, cex.axis=1.5, colours=c(ab_colours,rep('black',2)), returnPCA=TRUE)
dev.off()