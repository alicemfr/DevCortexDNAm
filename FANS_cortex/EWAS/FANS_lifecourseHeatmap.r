

## FANS heatmap ##

library(pheatmap)
library(data.table)
library(viridis)
library(RColorBrewer)

epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)

celltype_cols <- c(plasma(4)[2],viridis(4)[3])


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"FANSlifecourse_N_NN.rdat"))
plot_samp <- SampleSheet


#2. Remove X/Y chr ==============================================================================================================
probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicMan <- epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
XY <- which(probes$CHR  %in% c('X','Y'))
plot_betas <- betas[-XY,]


#3. Heatmap of adult samples at fetal cell-type DMPs ============================================================================
plot_samp <- SampleSheet
plot_betas <- betas[-XY,]
# keep only adult
samples <- which(plot_samp$Phenotype=='Adult')
plot_samp <- plot_samp[samples,]
plot_betas <- plot_betas[,samples]
identical(rownames(plot_samp),colnames(plot_betas))

# subset to Fetal Cell EWAS DMPs
dmps <- read.csv(paste0(AnalysisPath, "FACS_Cell_EWAS_fetal_anno_sig.csv"), row.names=1)
plot_betas.dmp <- plot_betas[rownames(plot_betas) %in% rownames(dmps),]
m <- as.matrix(plot_betas.dmp)

# create annotation matrix
nCuts <- nlevels(plot_samp$NewCellType)
annotation_col = data.frame(CellType=plot_samp$NewCellType, Antibody=plot_samp$Cell_Type)
rownames(annotation_col) <- rownames(plot_samp)

# run heatmap
pheatmap(m, annotation_col=annotation_col,  show_colnames=FALSE, show_rownames=FALSE, cutree_cols=nCuts, main = "", filename=paste0(plotPath,"FANS_heatmap_FetalCellDMPs_Adult.pdf"), width=10, height=8)