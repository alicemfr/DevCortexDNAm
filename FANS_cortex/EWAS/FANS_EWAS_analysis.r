

## Analyses the FANS lifecourse EWAS results ##


#1. Load libraries & define functions ===========================================================================================
library(scales)
library(stringr)
library(viridis)
library(data.table)
library(gridExtra)

'%ni%' <- Negate('%in%')

orderRes <- function(res, column){
	if(nrow(res)!=0){
		res <- res[order(res[,paste0(column)]),]
	}
	return(res)
}

#2. Load data ===================================================================================================================
load(paste0(MethylationPath,"FANSlifecourse_N_NN.rdat"))


#3. Load results ================================================================================================================
study <- 'Cell' #change to the EWAS study of interest
if(study=='AgeCell'){
	fl <- paste0(resPath, "FACS_", study, "_EWAS_fetal_anno")
}else{
	fl <- paste0(resPath, "FACS_", study, "_EWAS_fetal_adult_anno")
}
res <- data.frame(fread(paste0(fl,".csv"), data.table=F), row.names=1)
res.sig <- data.frame(fread(paste0(fl,"_sig.csv"), data.table=F), row.names=1)


#4. Create Sig columns ==========================================================================================================
# Fetal vs adult significance within cell type
if(study=='AgeCellSpecific'){
	celltype <- 'Non.neuronal' # options: c('Neuronal','Non.neuronal')
	res$Sig <- rep('Not significant',nrow(res))
	res$Sig[which(res[,paste0('Sig.Fetal.',celltype)]==TRUE)] <- 'Fetal'
	res$Sig[which(res[,paste0('Sig.Adult.',celltype)]==TRUE)] <- 'Adult'
	res$Sig[which(res[,paste0('Sig.Fetal.',celltype)]==TRUE & res[,paste0('Sig.Adult.',celltype)]==TRUE)] <- 'Both'
	res$Sig <- factor(res$Sig, levels=c("Not significant", "Fetal", "Adult", "Both"))
}else{
	res$Sig <- rep('Not significant',nrow(res))
	res$Sig[which(res[,'Fetal.Anova.P']<9e-8)] <- 'Fetal'
	res$Sig[which(res[,'Adult.Anova.P']<9e-8)] <- 'Adult'
	res$Sig[which(res[,'Fetal.Anova.P']<9e-8 & res[,'Adult.Anova.P']<9e-8)] <- 'Both'
	res$Sig <- factor(res$Sig, levels=c("Not significant", "Fetal", "Adult", "Both"))
}

# Neuronal vs non-neuronal significance within age group
if(study=='AgeCellSpecific'){
	res$Fetal.Sig <- rep('Not significant')
	res$Fetal.Sig[res[,paste0('Sig.Fetal.Neuronal')]==TRUE] <- 'Neuronal'
	res$Fetal.Sig[res[,paste0('Sig.Fetal.Non.neuronal')]==TRUE] <- 'Non-neuronal'
	res$Fetal.Sig[res[,paste0('Sig.Fetal.Neuronal')]==TRUE & res[,paste0('Sig.Fetal.Non.neuronal')]==TRUE] <- 'Both'
	res$Fetal.Sig <- factor(res$Fetal.Sig, levels=c('Not significant','Neuronal','Non-neuronal','Both'))
	
	res$Adult.Sig <- rep('Not significant')
	res$Adult.Sig[res[,paste0('Sig.Adult.Neuronal')]==TRUE] <- 'Neuronal'
	res$Adult.Sig[res[,paste0('Sig.Adult.Non.neuronal')]==TRUE] <- 'Non-neuronal'
	res$Adult.Sig[res[,paste0('Sig.Adult.Neuronal')]==TRUE & res[,paste0('Sig.Adult.Non.neuronal')]==TRUE] <- 'Both'
	res$Adult.Sig <- factor(res$Adult.Sig, levels=c('Not significant','Neuronal','Non-neuronal','Both'))
}


#5. Effect size scatter plots ===================================================================================================
source(paste0(ScriptsPath, "FANS_EWAS_plottingFunctions.r"))

# plot title & filename
if(study=='AgeCellSpecific'){
	ttl <- paste0(celltype, ' Age EWAS')
	flnm <- paste0(celltype,'.Age')
	
}else{	
	studies <- c('Age','Sex','Cell','AgeCell','AgeSex','SexCell')
	study_titles <- c('Age','Sex','CellType','Age*CellType','Age*Sex','Sex*CellType')
	sdy <- study_titles[which(studies==study)]
	ttl <- paste0(sdy, " EWAS")
	flnm <- study
}

# plot columns
if(study=='AgeCellSpecific'){
	col.fetal <- paste0('Fetal.',celltype,'.ES.Age')
	col.adult <- paste0('Adult.',celltype,'.ES.Age')
	sig.col.fetal <- paste0('Sig.Fetal.',celltype)
	sig.col.adult <- paste0('Sig.Adult.',celltype)
}else{
	col.fetal <- paste0('Fetal.Estimate.',study)
	col.adult <- paste0('Adult.Estimate.',study)
	sig.col.fetal <- 'Sig.Fetal'
	sig.col.adult <- 'Sig.Adult'
}


# all probes
df <- data.frame(Fetal=res[,col.fetal]*100, Adult=res[,col.adult]*100, Significance=res$Sig)
cor.all <- cor(df$Fetal, df$Adult)

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_colourSig.pdf"), width=8, height=8)
ESscatter_fetalVSadult(df, age.group='All', type='Point')
dev.off()

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_colourDensity2.pdf"), width=8, height=8)
ttl <- ''
ESscatter_fetalVSadult(df, age.group='All', type='2D', bins=240)
dev.off()


# fetal DMPs
df <- data.frame(Fetal=(res[which(res[,sig.col.fetal]==TRUE),col.fetal]), Adult=res[which(res[,sig.col.fetal]==TRUE),col.adult], Significance=res$Sig[which(res[,sig.col.fetal]==TRUE)])

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_sigFetal_colourSig.pdf"), width=8, height=8)
ESscatter_fetalVSadult(df, age.group='Fetal', type='Point')
dev.off()

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_sigFetal_colourSig.pdf"), width=8, height=8)
ESscatter_fetalVSadult(df, age.group='Fetal', type='2D')
dev.off()


# adult DMPs
df <- data.frame(Fetal=(res[which(res[,sig.col.adult]==TRUE),col.fetal]), Adult=res[which(res[,sig.col.adult]==TRUE),col.adult], Significance=res$Sig[which(res[,sig.col.adult]==TRUE)])

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_sigAdult_colourSig.pdf"), width=8, height=8)
ESscatter_fetalVSadult(df, age.group='Adult', type='Point')
dev.off()

pdf(paste0(plotPath,"EffectSize_", flnm, "_fetalVSadult_sigAdult_colourSig.pdf"), width=8, height=8)
ESscatter_fetalVSadult(df, age.group='Adult', type='2D')
dev.off()


# Effect size: Neuronal vs non-neuronal within age.group
if(study=='AgeCellSpecific'){
	
	age.group <- 'Fetal'
	Xlim <- c(-0.08,0.08)*100 #fetal: -0.08,0.08, adult: -0.005,0.005
	Ylim <- c(-0.08,0.08)*100
	
	eitherSig <- res[-which(res[,paste0(age.group,'.Sig')]=='Not significant'),] #2420
	
	# all DMPs
	df <- data.frame(Neuronal=eitherSig[,paste0(age.group,'.Neuronal.ES.Age')]*100, Non.neuronal=eitherSig[,paste0(age.group,'.Non.neuronal.ES.Age')]*100, Significance=eitherSig[,paste0(age.group,'.Sig')]) #2420
	Cor <- cor(df$Neuronal, df$Non.neuronal) #0.5642321
	
	pdf(paste0(plotPath,"EffectSize_", study, "_NvsNN_", age.group, "_anySig_colourSig.pdf"), width=8, height=8)
	ESscatter_AgeCellSpecific(df, age.group, Xlim, Ylim)
	dev.off()
	
	# cell-type specific DMPs
	celltype <- 'Non-neuronal'
	df2 <- df[df$Significance %in% c(celltype,'Both'),]
	Cor <- cor(df2$Neuronal, df2$Non.neuronal) #Neuronal=0.634, Non-neuronal=0.873
	
	pdf(paste0(plotPath,"EffectSize_", study, "_NvsNN_", age.group, "_", celltype, "Sig_colourSig.pdf"), width=8, height=8)
	ESscatter_AgeCellSpecific(df2, age.group, Xlim, Ylim)
	dev.off()
	
}


#6. Rescale ages for x-axis =====================================================================================================
fetal.ages.scale <- rescale(SampleSheet[which(SampleSheet$Phenotype=='Fetal'),'Age'], to = c(0, 33))
child.ages.scale <- rescale(SampleSheet[which(SampleSheet$Phenotype=='Child'),'Age'], to = c(41, 60))
adult.ages.scale <- rescale(SampleSheet[which(SampleSheet$Phenotype=='Adult'),'Age'], to = c(75, 100))
SampleSheet$age.rescale <- rep(NA, nrow(SampleSheet))
SampleSheet$age.rescale[which(SampleSheet$Phenotype=='Fetal')] <- fetal.ages.scale
SampleSheet$age.rescale[which(SampleSheet$Phenotype=='Child')] <- child.ages.scale
SampleSheet$age.rescale[which(SampleSheet$Phenotype=='Adult')] <- adult.ages.scale



#7. Plot ========================================================================================================================
#restore scientific notation for plotting p-value in titles
options(scipen=0)

source("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/plotAges.r") # scatter plot across life-course with loess fit line
source("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/plotBoxes.r") # boxplot of any phenotypic factor
source("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/fetal_plotFunctions.r") # fetal plotting functions


# Fetal vs Adult ----------------------------------------------------------------------------------------------------------------

age.group <- 'Fetal' # main age.group of interest

if(study=='AgeCellSpecific'){
	celltype <- 'Neuronal' # options: 'Neuronal', 'Non.neuronal'
	Pcol <- paste(age.group,celltype,'P.Age',sep='.')
	EScol <- paste(celltype,'ES.Age',sep='.')
	sigCol <- paste('Sig',age.group,celltype,sep='.')
}else{
	Pcol <- paste(age.group,'P',study,sep='.')
	EScol <- paste('Estimate',study,sep='.')
	sigCol <- paste('Sig',age.group,sep='.')
}

#1. DMPs
res.sig <- res[res[,sigCol]==TRUE,]

#2. DMPs ordered by p-value
res.ord <- orderRes(res.sig, Pcol)

#3. DMPs for one age.group only
res.sig.spec <- res.ord[which(res.ord$Sig==age.group),] # only probes where the age.group specified is significant, not both

#4. DMPs with opposite ES in fetal vs adult
both.sig <- res[which(res$Sig=='Both'),]
opp <- both.sig[which( sign(both.sig[,paste0('Fetal.',EScol)]) != sign(both.sig[,paste0('Adult.',EScol)]) ),] #Cell:1427

#5. Examples where ES is high in both fetal and adult
max.ES.f <- max(abs(opp[,paste0('Fetal.',EScol)]))
max.ES.a <- max(abs(opp[,paste0('Adult.',EScol)]))
opp.max.f <- opp[which(abs(opp[,paste0('Fetal.',EScol)])>=(max.ES.f-0.2*max.ES.f)),] #top 20% greatest fetal ES
opp.max.f.ord <- opp.max.f[order(abs(opp.max.f[,paste0('Adult.',EScol)]), decreasing=T),] #top adult ES of these
opp.max.a <- opp[which(abs(opp[,paste0('Adult.',EScol)])>=(max.ES.a-0.2*max.ES.a)),] #top 20% greatest adult ES
opp.max.a.ord <- opp.max.a[order(abs(opp.max.a[,paste0('Fetal.',EScol)]), decreasing=T),] #top fetal ES of these
opp.plot <- rbind(opp.max.f.ord[1:5,], opp.max.a.ord[1:5,])
opp.plot <- na.omit(opp.plot)


#6. Quadrants of ES scatters
both.sig <- res[which(res$Sig=='Both'),]
df <- both.sig
top.right <- df[which(df[,col.fetal]>0 & df[,col.adult]>0),]  # +ve fetal, +ve adult
bottom.right <- df[which(df[,col.fetal]>0 & df[,col.adult]<0),]  # +ve fetal, -ve adult
bottom.left <- df[which(df[,col.fetal]<0 & df[,col.adult]<0),]  # -ve fetal, -ve adult
top.left <- df[which(df[,col.fetal]<0 & df[,col.adult]>0),]  # -ve fetal, +ve adult
quad <- bottom.right # options: c('top.right','bottom.right','bottom.left','top.left')

# top examples where ES is high in both fetal and adult
max.ES.f <- max(abs(quad[,col.fetal]))
max.ES.a <- max(abs(quad[,col.adult]))
quad.max.f <- quad[which(abs(quad[,col.fetal])>=(max.ES.f-0.2*max.ES.f)),] #top 20% greatest fetal ES
quad.max.f.ord <- quad.max.f[order(abs(quad.max.f[,col.adult]), decreasing=T),] #top adult ES of these
quad.max.a <- quad[which(abs(quad[,col.adult])>=(max.ES.a-0.2*max.ES.a)),] #top 20% greatest adult ES
quad.max.a.ord <- quad.max.a[order(abs(quad.max.a[,col.fetal]), decreasing=T),] #top fetal ES of these
if(nrow(quad.max.f.ord)<5 | nrow(quad.max.a.ord)<5){ print('At least one of the dfs has <5 rows. Edit next line to account for this') }
F <- quad.max.f.ord[1:5,]
A <- quad.max.a.ord[1:5,]
if(any(rownames(F) %in% rownames(A))){ print('Some probes are in both lists. Remove duplicates before continuing') }
rownames(F)[which(rownames(F) %in% rownames(A))]
quad.plot <- rbind(F, A)

# Go to plots below...



# Neuronal vs non-neuronal ------------------------------------------------------------------------------------------------------
age.group <- 'Fetal'
Pcol.N <- paste0(age.group,'.Neuronal.P.Age')
EScol.N <- paste0(age.group,'.Neuronal.ES.Age')
sigCol.N <- paste0('Sig.',age.group,'.Neuronal')
Pcol.NN <- paste0(age.group,'.Non.neuronal.P.Age')
EScol.NN <- paste0(age.group,'.Non.neuronal.ES.Age')
sigCol.NN <- paste0('Sig.',age.group,'.Non.neuronal')

#1. Opposite ES
both.sig <- res[res[,sigCol.N]==TRUE & res[,sigCol.NN]==TRUE,]
opp <- both.sig[which( sign(both.sig[,EScol.N]) != sign(both.sig[,EScol.N]) ),]

#2. DMP in one cell-type only
neur <- res[which(res[,paste0(age.group,'.Sig')]=='Neuronal'),]
neur.sub <- neur[neur[,paste0(age.group,'.Non.neuronal.P.Age')]>0.05,] # maximise non-neuronal p-value
neur.ord <- neur.sub[order(neur.sub[,paste0(age.group,'.Neuronal.ES.Age')],decreasing=T),] # maximise neuronal age effect size

# Go to plots below...



# Life-course scatter plot ##----------------------------------------------------------------------------------------------------

colourBy <- 'NewCellType' # column indicating Neuronal/Non-neuronal
pchBy <- NULL
pch <- 16
colours <- c(plasma(4)[2],viridis(4)[3])
loessBy <- 'NewCellType'
p.age <- 'fetal'
res.plot <- neur.ord
res.column <- Pcol
if(nrow(res.plot)<30){ N <- nrow(res.plot) }else{ N <- 30 }

pdf(paste0(plotPath, "plot_",study,"_N_only.pdf"), height=8, width=12)
for(i in 1:N){
	plotAges(res=res.plot, pheno=SampleSheet, betas=betas, plot.column='age.rescale', res.column=res.column, p.age=p.age, i=i, gene=res.plot[i,'Gene'], pchBy=pchBy, pch=pch, chr=res.plot[i,'CHR'], colourBy=colourBy, colours=colours, xaxt='n', plot.loess=TRUE, loessSpan=NULL, loessBy=loessBy, pchByLegend=TRUE, loessByLegend=FALSE, colourByLegend=TRUE)
}
dev.off()



# Fetal-only scatter plot ##-----------------------------------------------------------------------------------------------------
colours <- c(plasma(4)[2],viridis(4)[3])
#res.sig.fetal <- res.sig[res$Sig.Fetal==TRUE,] # comment out if res.sig is already fetal-only
res.ord <- orderRes(res.sig, Pcol)
pheno.f <- SampleSheet[SampleSheet$Phenotype=='Fetal' & SampleSheet$AgeBin!='Late',]
betas.f <- betas[,colnames(betas) %in% rownames(pheno.f)]
identical(rownames(pheno.f),colnames(betas.f))

# extract types of interaction
Age.neg.Cell.neg <- res.ord[which(res.ord$Fetal.Estimate.CellType<0 & res.ord$Fetal.Estimate.Age<0),] # SYT1
Age.pos.Cell.neg <- res.ord[which(res.ord$Fetal.Estimate.CellType<0 & res.ord$Fetal.Estimate.Age>0),] # PEX14
Age.neg.Cell.pos <- res.ord[which(res.ord$Fetal.Estimate.CellType>0 & res.ord$Fetal.Estimate.Age<0),] # OAT
Age.pos.Cell.pos <- res.ord[which(res.ord$Fetal.Estimate.CellType>0 & res.ord$Fetal.Estimate.Age>0),] # GJA1
#poi <- which(rownames(res.ord) %in% rownames(Age.neg.Cell.neg)[1:50])

# For supp fig
goi <- c('SYT1', 'PEX14', 'OAT', 'GJA1') # gene of interest - examples from each group above
poi <- c() # position of interest
for(i in 1:length(goi)){
	poi <- c(poi, which(res.ord$Gene %in% goi[i])[1]) # first instance of gene
}

# For main fig
poi <- 1 #cg04831806 - SDK1

myplots <- list()
for(j in 1:length(poi)){
	i <- poi[j]
	p <- ageByGroup_scatter(res.ord, pheno.f, betas.f, i, age.column='Age', p.column=Pcol, group.column='NewCellType', group.plot.name='Cell-type', xlab="Age (pcw)", ylab="DNA methylation (%)", colours)
	myplots[[j]] <- p
}

pdf(paste0(plotPath, "plot_",study,"_Fetal_4examples.pdf"), width=7, height=7)
myplots
dev.off()

pdf(paste0(plotPath, "plot_",study,"_Fetal_topDMP_SDK1.pdf"), width=7, height=7)
myplots
dev.off()


# Life-course boxplot ##---------------------------------------------------------------------------------------------------------
colourBy <- 'NewCellType'
pchBy <- NULL
pch <- 16
colours <- c(viridis(4)[3],plasma(4)[2]) # needs to be reverse of above
loessBy <- 'NewCellType'
p.age <- 'fetal'
res.plot <- opp.plot
res.column <- paste0(age.group,'.P.',study)
if(nrow(res.plot)<30){ N <- nrow(res.plot) }else{ N <- 30 }

pdf(paste0(plotPath, "boxplot_",study,"_oppES.pdf"), height=8, width=12)
for(i in 1:N){
	plotBoxes(res=res.plot, pheno=SampleSheet, betas=betas, plot.column=colourBy, res.column=res.column, p.age=p.age, i=i, gene=res.plot[i,'Gene'], pchBy=pchBy, pch=pch, chr=res.plot[i,'CHR'], colours=colours, inter=TRUE, interCol='Phenotype', subsetTo=c('Fetal','Adult'), subsetToCol='Phenotype', boxNames=c('Neuronal','Non-neuronal','Neuronal','Non-neuronal'), labsPos=c(1.5,3.5))
}
dev.off()