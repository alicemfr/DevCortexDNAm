

## Analysis of adult (from the bulk lifecourse dataset) age regression results run on bulk fetal DMPs.


library(ggplot2)
library(data.table)
library(scales)
library(plyr) # for ddply
library(dplyr)
library(viridis)


#1. Load data ===================================================================================================================
load(paste0(MethylationPath,"EPICBrainLifecourse.rdat"))

betas <- epic.betas
pheno <- epic.pheno

fetal.ages.scale <- rescale(pheno[which(pheno$Phenotype=='Fetal'),'Age'], to = c(0, 33))        # max fetal age is 33pcw
child.ages.scale <- rescale(pheno[which(pheno$Phenotype=='Child'),'Age.years'], to = c(41, 55)) # fetal ends and child starts at 40pcw (0 years)
adult.ages.scale <- rescale(pheno[which(pheno$Phenotype=='Adult'),'Age.years'], to = c(65, 100))
pheno$age.rescale <- rep(NA, nrow(pheno))
pheno$age.rescale[which(pheno$Phenotype=='Fetal')] <- fetal.ages.scale
pheno$age.rescale[which(pheno$Phenotype=='Child')] <- child.ages.scale
pheno$age.rescale[which(pheno$Phenotype=='Adult')] <- adult.ages.scale


#2. Load prenatal results =======================================================================================================
res.fetal <- data.frame(readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")))


#3. Load postnatal results ======================================================================================================
res <- data.frame(readRDS(paste0(analysisPath, "postnatalAgeReg_bulkFetalDMPs.rds"))) #41518
colnames(res) <- c('Beta.Age','SE.Age','P.Age', 'Beta.Sex','SE.Sex','P.Sex')

# match to fetal results
res.fetal <- res.fetal[which(res.fetal$P.Age<9e-8),] #50913
res.fetal.sub <- res.fetal[which(rownames(res.fetal) %in% rownames(res)),] #41518
res <- res[which(rownames(res) %in% rownames(res.fetal.sub)),] #41518
identical(rownames(res.fetal.sub), rownames(res))
perc <- (nrow(res)/nrow(res.fetal))*100 #81.54695%

# annotate EWAS results table
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE)
epicMan<-epicManifest[match(rownames(res), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
res <- cbind(res, as.data.frame(epicMan))

# extract unique mentions of genes
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))


#4. Combine pre and postnatal results ===========================================================================================
res.sig <- res                                           # postnatal results table (all of which are significant in fetal)
res.fetal.sub.sig <- res.fetal.sub                       # fetal results table (all of which are significant in fetal)
identical(rownames(res.sig),rownames(res.fetal.sub.sig))

res.all.sig <- cbind(Post=res.sig[,c('Beta.Age','SE.Age','P.Age', 'Beta.Sex','SE.Sex','P.Sex')], Pre=res.fetal.sub.sig[,c('Beta.Age','SE.Age','P.Age', 'Beta.Sex','SE.Sex','P.Sex')], res.sig[,c('CHR','Gene')]) #41518


#5. Characterise level of postnatal significance ================================================================================
bonf <- 0.05/nrow(res) #Bonferroni P = 1.204297e-06

res.all.sig$bonfSig <- rep(0,nrow(res.all.sig))
res.all.sig$bonfSig[res.all.sig$Post.P.Age<bonf] <- 1 #1003

res.all.sig$Significance <- rep("> 0.05",nrow(res.all.sig))
res.all.sig$Significance[which(res.all.sig$Post.P.Age<0.05)] <- "< 0.05"
res.all.sig$Significance[which(res.all.sig$Post.P.Age<bonf)] <- "< bonferroni"
res.all.sig$Significance[which(res.all.sig$Post.P.Age<9e-8)] <- "< genome-wide"
res.all.sig$Significance <- as.factor(res.all.sig$Significance)

Cor <- cor(res.all.sig$Post.Beta.Age, res.all.sig$Pre.Beta.Age) #-0.009251818. Correlation of all prenatal DMPs. 
postsig <- which(res.all.sig$bonfSig==1) #1003
perc <- (length(postsig)/nrow(res.all.sig))*100 #2.41582%. Percentage of prenatal DMPs that are also significant in postnatal samples.
Cor.postsig <- cor(res.all.sig[postsig,'Post.Beta.Age'],res.all.sig[postsig,'Pre.Beta.Age']) #0.1126558. Correlation of probes that are significant in post and prenatal.
postsig.post.sign <- sign(res.all.sig[postsig,'Post.Beta.Age']) # within the postnatal DMPs, sign of beta for postnatal samples
postsig.pre.sign <- sign(res.all.sig[postsig,'Pre.Beta.Age']) # within the postnatal DMPs, sign of beta for prenatal samples
table(postsig.post.sign, postsig.pre.sign)
                 # postsig.pre.sign
# postsig.post.sign  -1   1
               # -1  64  79
               # 1  196 664
			   
hyper.hyper <- length(which(postsig.post.sign==1 & postsig.pre.sign==1)) #664
hyper.hyper.perc.ofHyperdDMPs <- (hyper.hyper/length(which(postsig.pre.sign==1)))*100 #89.36743%. hyper prenatal dDMPs that are also sig hyper postnatally

hypo.hypo <- length(which(postsig.post.sign==(-1) & postsig.pre.sign==(-1))) #64
hypo.hypo.perc.ofHypodDMPs <- (hypo.hypo/length(which(postsig.pre.sign==(-1))))*100 #24.61538%

hyper.hypo <- length(which(postsig.post.sign==1 & postsig.pre.sign==(-1))) #196
hypo.hyper <- length(which(postsig.post.sign==(-1) & postsig.pre.sign==1)) #79
opp <- which(postsig.post.sign != postsig.pre.sign) #275
perc <- (length(opp)/length(postsig))*100 #27.41775%

levels(res.all.sig$Significance) <- c("0.05 > P > bonferroni", "bonferroni > P > genome-wide", "P < genome-wide", "P > 0.05")
res.all.sig$Significance <- factor(res.all.sig$Significance, levels=c("P > 0.05", "0.05 > P > bonferroni", "bonferroni > P > genome-wide", "P < genome-wide"))


#6. Effect size scatter plot: pre vs postnatal for prenatal dDMPs, coloured by postnatal significance ===========================
ES.scatter <- function(res, plot.col1, plot.col2, xlab, ylab, main, pos.legend='topleft', cex.legend=1.2, cex.main=1.8, cex.lab=1.5, cex.axis=1.5, colBySigRes=TRUE, colBySigCol='bonfSig', colBySigThresh=9e-8, colNonSig='lightgrey', colSig=c("midnightblue", "skyblue"), sigPcol='Post.P.Age', colBySigLabels=NULL, xlim=NULL, ylim=NULL, colour='black', colourOrder=NULL, includeCor=FALSE, specificLevels=NULL, plotLegend=TRUE){


	df <- data.frame(x=res[,plot.col1]*100, y=res[,plot.col2]*100, significance=res[,colBySigCol], P=res[,sigPcol])
	df.nonsig <- df[which(df$sig==0),]
	df.sig <- df[which(df$sig==1),]
	df.sig <- df.sig[order(df.sig[,'P']),]

	# calc p-value quartiles for legend labels
	logP <- -log10(df.sig[,'P']) # convert to -log10(p) to calculate quartiles
	Q100 <- max(logP) # most significant
	Q0 <- min(logP) # least significant
	Q50 <- Q0 + ((Q100-Q0)/2)
	Q75 <- Q50 + ((Q100-Q50)/2)
	Q25 <- Q50 - ((Q100-Q50)/2)
	label.names <- c(Q0, Q25, Q50, Q75, Q100)
	label.names <- signif(10**(-label.names), digits=2) # return to original p-values

	colfunc <- colorRampPalette(colSig) # colours to range from
	legend_image <- as.raster(matrix(colfunc(nrow(df.sig)), ncol=1)) # establish colour gradient

	par(mar=c(6, 6, 4.1, 2.1))

	#1. plot non-significant sites
	plot(df.nonsig$x, df.nonsig$y, xlab='', ylab=ylab, main=main, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, col=colNonSig, xlim=xlim, ylim=ylim)
	mtext(text=xlab, side=1, line=4, cex=cex.axis) 

	#2. plot significant sites
	points(x=df.sig$x, y=df.sig$y, xlab='Effect Size - Prenatal', ylab='Effect Size - Postnatal', main='', col=colfunc(nrow(df.sig)), cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
	abline(v=0, h=0);abline(0,1, lty=3)

	#3. Plot legend
	# Edited colour gradient legend code snippet from: https://stackoverflow.com/questions/13355176/gradient-legend-in-base/70522655#70522655
	op <- par(fig=c(0.65,0.9, 0.5,0.96), new=TRUE)
	plot(c(0.5, 1), c(0, 1), type='n', axes=FALSE, xlab='', ylab='') # position of legend
	rasterImage(legend_image, 0, 0, 1, 1) # fill legend with colour gradient
	label.pos <- seq(0, 1, l=5)    
	axis(4, at=label.pos, pos=1, labels=F, col=0, col.ticks=1, tck=-.1)  ## axis ticks
	mtext(label.names, 4, 0.5, at=label.pos, las=2, cex=0.8)                    ## tick labels
	mtext('Postnatal\nage p-value', 3, -.125, cex=0.9, adj=.1, font=2)          ## legend title

	par(op) #reset par
	par(mfrow=c(1,1))

}

pdf(paste0(AnalysisPath, "adultDMPs_fetalDMPs_BULKlifecourse_allProbes_colourGradPval_blue_3.pdf"), width=8, height=8)
	ES.scatter(res=res.all.sig, plot.col1='Pre.Beta.Age', plot.col2='Post.Beta.Age', xlab='Prenatal\n% change per week', ylab='Postnatal\n% change per year', main='', xlim=c(-0.06,0.08)*100, ylim=c(-0.002,0.0025)*100)
dev.off()


#7. Examples from quadrants =====================================================================================================
source(paste0(scriptsPath, "plotAges.r"))

runPlotAges <- function(res.plot, p.age, filename, nPlots){
	pdf(paste0(AnalysisPath, filename), height=8, width=12)
	extraDetailsLabels=c(min(pheno[which(pheno$Phenotype=='Fetal'),'Age']),c('40 / 0 ',signif(max(pheno[which(pheno$Phenotype=='Child'),'Age.years']),digits=1),min(pheno[which(pheno$Phenotype=='Adult'),'Age.years'])),max(pheno[which(pheno$Phenotype=='Adult'),'Age.years']))
	for(i in 1:nPlots){
		plotAges(res=res.plot, pheno=pheno, betas=betas, plot.column='age.rescale', res.column='Post.P.Age', p.age=p.age, i=i, gene=res.plot[i,'Gene'], chr=res.plot[i,'CHR'], pchByLegend=FALSE, loessByLegend=FALSE, colours='black', colourBy=NULL, colourByLegend=FALSE, xaxt='n', extendXlim=FALSE, extraDetailsLabels=extraDetailsLabels, axisPos=c(40,55,65), cexLab=1.5, cexAxis=1.3)
	}
	dev.off()
}

df <- res.all.sig[which(res.all.sig$Significance=='P < genome-wide'),] # limit to dDMPs that are also significant postnatally
col.fetal <- 'Pre.Beta.Age'
col.adult <- 'Post.Beta.Age'

top.right <- df[which(df[,col.fetal]>0 & df[,col.adult]>0),]     # +ve fetal, +ve adult
bottom.right <- df[which(df[,col.fetal]>0 & df[,col.adult]<0),]  # +ve fetal, -ve adult
bottom.left <- df[which(df[,col.fetal]<0 & df[,col.adult]<0),]   # -ve fetal, -ve adult
top.left <- df[which(df[,col.fetal]<0 & df[,col.adult]>0),]      # -ve fetal, +ve adult

perc.top.right <- (nrow(top.right)/nrow(df))*100 #74.44444%
perc.top.left <- (nrow(top.left)/nrow(df))*100 #15.41667%
perc.bottom.right <- (nrow(bottom.right)/nrow(df))*100 #2.5%
perc.bottom.left <- (nrow(bottom.left)/nrow(df))*100 #7.638889%

quad <- top.left # options: c('top.right','bottom.right','bottom.left','top.left')

# examples where ES is large in both fetal and adult
max.ES.f <- max(abs(quad[,col.fetal]))
max.ES.a <- max(abs(quad[,col.adult]))
quad.max.f <- quad[which(abs(quad[,col.fetal])>=(max.ES.f-0.2*max.ES.f)),] #top 20% greatest fetal ES
quad.max.f.ord <- quad.max.f[order(abs(quad.max.f[,col.adult]), decreasing=T),] #top adult ES of these
quad.max.a <- quad[which(abs(quad[,col.adult])>=(max.ES.a-0.2*max.ES.a)),] #top 20% greatest adult ES
quad.max.a.ord <- quad.max.a[order(abs(quad.max.a[,col.fetal]), decreasing=T),] #top fetal ES of these

F <- quad.max.f.ord[1:5,]
A <- quad.max.a.ord[1:5,]
quad.plot <- rbind(F, A)

# plot
runPlotAges(res.plot=quad.plot, p.age='postnatal', filename="ageReg_PrevsPost_topleft_maxESmethod.pdf", nPlots=nrow(quad.plot))