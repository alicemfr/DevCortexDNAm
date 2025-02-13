

## Effect size scatter plots for bulk/neuronal/non-neuronal Age effect sizes ##

library(ggplot2)
library(viridis)
library(grid)
library(VennDiagram)

'%ni%' <- Negate('%in%')

ES.plot <- function(res.df, col1, col2, xlim, ylim, xlab, ylab, scale.factor=100){
	X <- res.df[,col1]
	Y <- res.df[,col2]
	df <- data.frame(X,Y)
	cor.xy <- signif(cor.test(X,Y)$estimate, digits=3)

	ggplot(df, aes(x=X, y=Y) ) +
	  geom_vline(xintercept = 0)+ geom_hline(yintercept = 0)+
	  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth=0.5, fullrange=TRUE)+
	  geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=0.5)+
	  geom_bin2d(bins = nrow(df)/scale.factor) +
	  scale_fill_continuous(type = "viridis") +
	  xlab(xlab) + ylab(ylab) +
	  xlim(xlim)+ ylim(ylim)+
	  theme_minimal() +
	  theme(axis.text=element_text(size=19), axis.title=element_text(size=22), plot.title = element_text(size=23))
}


#1. Load results ================================================================================================================
bulk <- as.data.frame(readRDS(paste0(bulkRes,"ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")))
fans <- as.data.frame(readRDS(paste0(fansRes,"FANS_AgeCellSpecific_EWAS_Fetal.rds")))


#2. Bulk DMPs: bulk vs neuronal & non-neuronal ==================================================================================
# limit to overlapping probes
bulk.dmps <- bulk[which(bulk$P.Age<9e-8),] #50913 - total number of bulk dmps
fans.bulk.dmps <- fans[which(rownames(fans) %in% rownames(bulk.dmps)),] #42114 - number of bulk dmps we have data for in the fans dataset
bulk.dmps <- bulk.dmps[match(intersect(rownames(fans.bulk.dmps), rownames(bulk.dmps)), rownames(bulk.dmps)),]
identical(rownames(bulk.dmps),rownames(fans.bulk.dmps))

bulk.ageES <- bulk.dmps$Beta.Age*100 # multiply by 100 to convert to % DNAm change
neur.ageES <- fans.bulk.dmps$Neuronal.ES.Age*100
non_neur.ageES  <- fans.bulk.dmps$Non.neuronal.ES.Age*100

df <- data.frame(Bulk=bulk.ageES, Neur=neur.ageES, Non_neur=non_neur.ageES)
rownames(df) <- rownames(bulk.dmps)

cor.bulk.neur <- cor(df$Bulk, df$Neur) #0.9233224
cor.bulk.non.neur <- cor(df$Bulk, df$Non_neur) #0.8695539
cor.neur.non.neur <- cor(df$Neur, df$Non_neur) #0.7734067

pdf(paste0(AnalysisPath,"bulkDMPs_AgeES_bulk_vs_neuronal.pdf"), width=7, height=7)
ES.plot(df, col1='Bulk', col2='Neur', xlim=c(-8,8), ylim=c(-8,8), xlab="Bulk cortex\nAge effect size (%)", ylab="Neuronal cell-type\nAge effect size (%)")
dev.off()

pdf(paste0(AnalysisPath,"bulkDMPs_AgeES_bulk_vs_non.neuronal.pdf"), width=7, height=7)
ES.plot(df, col1='Bulk', col2='Non_neur', xlim=c(-8,8), ylim=c(-8,8), xlab="Bulk cortex\nAge effect size (%)", ylab="Non-neuronal cell-types\nAge effect size (%)")
dev.off()

pdf(paste0(AnalysisPath,"bulkDMPs_AgeES_neuronal_vs_non.neuronal.pdf"), width=7, height=7)
ES.plot(df, col1='Neur', col2='Non_neur', xlim=c(-8,8), ylim=c(-8,8), xlab="Neuronal cell-type\nAge effect size (%)", ylab="Non-neuronal cell-types\nAge effect size (%)")
dev.off()

# examples from 'bulkDMPs_AgeES_neuronal_vs_non.neuronal' that show large effect in one cell-type and not other
# selecting thresholds based on figure
Nspec <- df[abs(df$Neur)>5 & abs(df$Non_neur)<1,] #52. Neuronal-specific: large ES in neuronal, small ES in non-neuronal
NNspec <- df[abs(df$Non_neur)>3 & abs(df$Neur)<1,] #28. Non-neuronal-specific: large ES in non-neuronal, small ES in neuronal
Nneg_NNpos <- df[df$Neur<(-1) & df$Non_neur>1,] #15. Neuronal -ve ES, non-neuronal +ve ES
Npos_NNneg <- df[df$Neur>1 & df$Non_neur<(-1),] #1. Neuronal +ve ES, non-neuronal -ve ES


#3. FANS DMPs: neuronal & non-neuronal ==========================================================================================
neur.dmps <- fans[which(fans$Neuronal.P.Age<9e-8),] #1872
non_neur.dmps <- fans[which(fans$Non.neuronal.P.Age<9e-8),] #820

df.neur <- data.frame(Neur=neur.dmps$Neuronal.ES.Age*100, Non_neur=neur.dmps$Non.neuronal.ES.Age*100)
df.non_neur <- data.frame(Neur=non_neur.dmps$Neuronal.ES.Age*100, Non_neur=non_neur.dmps$Non.neuronal.ES.Age*100)

cor.neur.dmps <- cor(df.neur$Neur, df.neur$Non_neur) #0.6335795
cor.non_neur.dmps <- cor(df.non_neur$Neur, df.non_neur$Non_neur) #0.8732221

pdf(paste0(fansPlots,"neuronalDMPs_AgeES_neuronal_vs_non.neuronal.pdf"), width=7, height=7)
ES.plot(df.neur, col1='Neur', col2='Non_neur', xlim=c(-8,8), ylim=c(-8,8), xlab="Neuronal cell-type\nAge effect size (%)", ylab="Non-neuronal cell-type\nAge effect size (%)", scale.factor=20)
dev.off()

pdf(paste0(fansPlots,"non.neuronalDMPs_AgeES_neuronal_vs_non.neuronal.pdf"), width=7, height=7)
ES.plot(df.non_neur, col1='Non_neur', col2='Neur', xlim=c(-8,8), ylim=c(-8,8), xlab="Non-neuronal cell-type\nAge effect size (%)", ylab="Neuronal cell-type\nAge effect size (%)", scale.factor=8)
dev.off()


#4. Neuronal vs non-neuronal Venn diagram =======================================================================================
setwd(fansPlots)
venn.diagram(x=list(rownames(neur.dmps), rownames(non_neur.dmps)), category.names = c("Neuronal DMPs" , "Non-neuronal DMPs"), filename = 'FACS_AgeCellSpecific_EWAS_fetal_Venn.png', output=TRUE, cat.pos=-1)


#5. Bulk / neuronal / non-neuronal Venn diagram =================================================================================
setwd(AnalysisPath)

# re-define neuronal and non-neuronal DMPs based on number of bulk DMPs. p-value threshold = bonferroni of 0.05/nbulkDMPs
pthresh <- 0.05/nrow(bulk.dmps) #1.187254e-06
neur.dmps <- fans[which(fans$Neuronal.P.Age<pthresh),] #4088
non_neur.dmps <- fans[which(fans$Non.neuronal.P.Age<pthresh),] #1820

venn.diagram(x=list(rownames(bulk.dmps), rownames(neur.dmps), rownames(non_neur.dmps)), category.names = c("Bulk DMPs", "Neuronal DMPs" , "Non-neuronal DMPs"), filename = 'bulkDMPs_AgeCellSpecificDMPs_Venn.png', output=TRUE, cat.pos=c(-30,30,180))

NinBulk <- length(which(rownames(neur.dmps) %in% rownames(bulk.dmps))) #3574
NNinBulk <- length(which(rownames(non_neur.dmps) %in% rownames(bulk.dmps))) #1766
percNinBulk <- (NinBulk/nrow(neur.dmps))*100 #87.42661
percNNinBulk <- (NNinBulk/nrow(non_neur.dmps))*100 #97.03297

BulkinNandNN <- length(which(rownames(bulk.dmps) %in% rownames(neur.dmps) & rownames(bulk.dmps) %in% rownames(non_neur.dmps))) #774
percBulkinNandNN <- (BulkinNandNN/nrow(bulk.dmps))*100 #1.837869


#6. Bulk / neuronal / non-neuronal specific DMPs ================================================================================

# make combined phenotype file and betas matrix
load(paste0(bulkMethylation, "fetalBulk_EX3_23pcw_n91.rdat")) # bulk fetal DNAm
bulk.betas <- betas
bulk.pheno <- pheno
bulk.pheno$CellType <- rep('Bulk',nrow(bulk.pheno))

load(paste0(fansMethylation,"FANSlifecourse_N_NN.rdat")) # FANS lifecourse DNAm
fans.betas <- betas[,SampleSheet$Phenotype=='Fetal']
fans.pheno <- SampleSheet[SampleSheet$Phenotype=='Fetal',]
identical(rownames(fans.pheno),colnames(fans.betas))

# match betas
bulk.betas.m <- bulk.betas[rownames(bulk.betas) %in% rownames(fans.betas),]
fans.betas.m <- fans.betas[rownames(fans.betas) %in% rownames(bulk.betas.m),]
identical(rownames(bulk.betas.m),rownames(fans.betas.m))

# cbind fans samples to bulk betas
betas.combo <- cbind(bulk.betas.m, fans.betas.m)
b <- bulk.pheno[,c('Sample_ID','PCW','CellType')]
colnames(b) <- c('Sample_ID','Age','CellType')
f <- fans.pheno[,c('Sample_ID','Age','NewCellType')]
colnames(f) <- c('Sample_ID','Age','CellType')
pheno.combo <- rbind(b,f)
pheno.combo$Age <- as.numeric(pheno.combo$Age)
pheno.combo$CellType <- as.factor(pheno.combo$CellType)
identical(rownames(pheno.combo),colnames(betas.combo))


# extract example probes to plot
bulk.dmps <- bulk[which(bulk$P.Age<9e-8),] #50913 - total number of bulk dmps
fans.bulk.dmps <- fans[which(rownames(fans) %in% rownames(bulk.dmps)),] #42114 - number of bulk dmps we have data for in the fans dataset
bulk.dmps <- bulk.dmps[match(intersect(rownames(fans.bulk.dmps), rownames(bulk.dmps)), rownames(bulk.dmps)),]
identical(rownames(bulk.dmps),rownames(fans.bulk.dmps))

# re-calc p-thresh
pthresh <- 0.05/nrow(bulk.dmps) #1.187254e-06
neur.dmps <- fans[which(fans$Neuronal.P.Age<pthresh),] #4088
non_neur.dmps <- fans[which(fans$Non.neuronal.P.Age<pthresh),] #1820

# combine data
df <- data.frame(bulk[rownames(bulk) %in% rownames(fans),c('Beta.Age','P.Age','CHR','MAPINFO','Gene')], fans[rownames(fans) %in% rownames(bulk),c('Neuronal.ES.Age','Neuronal.P.Age','Non.neuronal.ES.Age','Non.neuronal.P.Age')]) #692927
df$Bulk <- df$P.Age<9e-8 #42114
df$Neur <- df$Neuronal.P.Age<pthresh #4086 - 2 probes fewer than before because 2 probes not present in bulk dataset
df$Non.neur <- df$Non.neuronal.P.Age<pthresh #1818 - 2 probes fewer than before because 2 probes not present in bulk dataset

colnames(df)[which(colnames(df) %in% c('P.Age','Neuronal.P.Age','Non.neuronal.P.Age'))] <- c('P.B','P.N','P.NN')
colnames(df)[which(colnames(df) %in% c('Beta.Age','Neuronal.ES.Age','Non.neuronal.ES.Age'))] <- c('ES.B','ES.N','ES.NN')

df$BminusN <- df$ES.B - df$ES.N
df$BminusNN <- df$ES.B - df$ES.NN
df$NminusNN <- df$ES.N - df$ES.NN

all.dmps <- df[df$Neur==T & df$Non.neur==T & df$Bulk==T,] #774
A <- all.dmps[order(abs(all.dmps$ES.B),decreasing=T),] # order by bulk age ES: largest to smallest 

neur.spec.dmps <- df[df$Neur==T & df$Non.neur==F & df$Bulk==F,] #511
bigN <- neur.spec.dmps[abs(neur.spec.dmps$BminusN)>0.025 & abs(neur.spec.dmps$NminusNN)>0.025,]

non_neur.spec.dmps <- df[df$Neur==F & df$Non.neur==T & df$Bulk==F,] #51
bigNN <- non_neur.spec.dmps[abs(non_neur.spec.dmps$BminusNN)>0.02 & abs(non_neur.spec.dmps$NminusNN)>0.02,]

bulk.spec.dmps <- df[df$Neur==F & df$Non.neur==F & df$Bulk==T,] #37548
bigB <- bulk.spec.dmps[abs(bulk.spec.dmps$BminusN)>0.02 & abs(bulk.spec.dmps$BminusNN)>0.02,]

largeN <- df[rownames(df) %in% rownames(Nspec),]
largeN <- largeN[order(abs(largeN$NminusNN),decreasing=T),]

largeNN <- df[rownames(df) %in% rownames(NNspec),]
largeNN <- largeNN[order(abs(largeNN$NminusNN),decreasing=T),]

Nneg_NNpos <- df[rownames(df) %in% rownames(Nneg_NNpos),]
Npos_NNneg <- df[rownames(df) %in% rownames(Npos_NNneg),]

# plot examples
source(paste0(scriptsPath, "fetal_plotFunctions.r"))

colours <- c('black',plasma(4)[2],viridis(4)[3])
toPlot <- 'Npos_NNneg' # options = 'largeN', 'largeNN', 'Nneg_NNpos', 'Npos_NNneg'
res.plot <- get(toPlot)
p.col <- 'P.NN'

myplots <- list()
for(i in 1:20){
	p <- ageByGroup_scatter(res.plot, pheno.combo,betas.combo, i, age.column='Age', p.column=p.col, group.column='CellType', group.plot.name='Cell-type', xlab="Age (pcw)", ylab="DNA methylation (%)", colours, Ylim=c(-45,100))
	myplots[[i]] <- p
}
pdf(paste0(AnalysisPath, "plot_B_N_NN.", toPlot, "_fullAge.pdf"), width=8, height=8)
myplots
dev.off()