

## Comparison of bulk fetal "linear" and "non-linear" probe lists ##

library(data.table)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(VennDiagram)


#1. Load results ================================================================================================================
reg <- readRDS(paste0(linearPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")) #807806. linear regression stats for all probes
gpstats <- fread(paste0(GPPath, "gpstats.tsv"), header=T, sep='\t', data.table=F) #543796. Gaussian process stats for all probes
nonlin.probes <- readRDS(paste0(nonlinearPath, "nonlinearMetrics.rds")) #73035. list of nonlinear probes as defined by Gaussian process model

#2. Subset EWAS stats for linear & nonlinear probes =============================================================================
# linear regression stats
lin <- reg[which(reg$P.Age<9e-8),] #50913. Linear regression stats for the linear probes
lin.non <- reg[which(rownames(reg) %in% rownames(nonlin.probes)),] #73035. Linear regression stats for the nonlinear probes

# nonlinear stats
non.lin <- gpstats[which(gpstats$Probe %in% rownames(lin)),] #49597. GP stats for linear probes
non <- gpstats[which(gpstats$Probe %in% rownames(nonlin.probes)),] #73035. GP stats for nonlinear probes

# combine stats
linear <- data.frame(lin[which(rownames(lin) %in% non.lin$Probe),c('Beta.Age','P.Age')], non.lin[,c('PreferedModelLR','mat52_ell')],Group=rep('Linear',nrow(non.lin))) #49597
nonlinear <- data.frame(lin.non[,c('Beta.Age','P.Age')], non[,c('PreferedModelLR','mat52_ell')],Group=rep('Nonlinear',nrow(non))) #73035
df <- rbind(linear, nonlinear)

# add linear effect size and p-value to nonlinear probes file
nonlinearProbes <- read.csv(paste0(nonlinearPath, "nonlinearProbes.csv"),row.names=1)
df.m <- df[rownames(df) %in% rownames(nonlinearProbes),]
df.m <- df.m[match(rownames(nonlinearProbes),rownames(df.m)),]
df.m$dDMP <- df.m$P.Age<9e-8
nonlinearProbes2 <- cbind(nonlinearProbes[,c('IlmnID','CHR','UCSC_RefGene_Group','Gene')], df.m[,c('Beta.Age','P.Age','dDMP')], nonlinearProbes[,c('Module','Timescale','logTimescale','MatvsConst','MatvsLinear','logSNR')])
colnames(nonlinearProbes2)[colnames(nonlinearProbes2) %in% c('Beta.Age','P.Age')] <- c('Linear.ES','Linear.P')
write.csv(nonlinearProbes2, file=paste0(nonlinearPath, "nonlinearProbes_plusLinearStats.csv"))

# percentage of significant linear & nonlinear probes present in each results list
overlap <- length(which(rownames(lin) %in% rownames(nonlin.probes))) #23252
perc.lin.in.non <- (overlap/nrow(nonlin.probes))*100 #31.83679% of nonlinear probes are also significantly linear
perc.non.in.lin <- (overlap/nrow(lin))*100 #45.67006% of significantly linear probes are also nonlinear

# percentage of non-significant linear & nonlinear probes present in each results list
nonsig.lin <- reg[which(reg$P.Age>=9e-8),] #756893
overlap <- length(which(rownames(nonsig.lin) %in% rownames(nonlin.probes))) #49783
perc.lin.in.non <- (overlap/nrow(nonlin.probes))*100 #68.16321% of nonlinear probes are also nonsignificantly linear
perc.non.in.lin <- (overlap/nrow(nonsig.lin))*100 #6.577284% of nonsignificantly linear probes are also nonlinear


#3. Compare nonlinear and linear stats ==========================================================================================
Beta.Age <- ggplot(df, aes(x=Group, y=Beta.Age))+ geom_jitter(aes(colour=Group))
P.Age <- ggplot(df, aes(x=Group, y=-log10(P.Age)))+ geom_jitter(aes(colour=Group))+geom_hline(yintercept=-log10(9e-8), linetype='dotted', size=1.2)
const_mll <- ggplot(df, aes(x=Group, y=const_mll))+ geom_jitter(aes(colour=Group))
linear_mll <- ggplot(df, aes(x=Group, y=linear_mll))+ geom_jitter(aes(colour=Group))
mat52_mll <- ggplot(df, aes(x=Group, y=mat52_mll))+ geom_jitter(aes(colour=Group))
mat52_ell <- ggplot(df, aes(x=Group, y=mat52_ell))+ geom_jitter(aes(colour=Group))
grid.arrange(Beta.Age,P.Age,const_mll,linear_mll,mat52_mll,mat52_ell, nrow=2, ncol=3)

tb1 <- table(df$Group,df$PreferedModelLR)
tb <- melt(tb1)
colnames(tb) <- c('Original','Predicted','Frequency')
ggplot(tb, aes(fill=Predicted, y=Frequency, x=Original))+ geom_bar(position="stack", stat="identity")
	
            # Linear Mat52
  # Linear      7316 42281
  # Nonlinear      0 73035
  
linpred <- tb1['Linear',]
mat52pred <- tb1[,'Mat52']
perc.linpred <- as.numeric((linpred[2]/sum(linpred)))*100 #85.24911%. percentage linear probes predicted nonlinear
perc.mat52pred <- as.numeric((mat52pred[1]/sum(mat52pred)))*100 #36.66534%. percentage of the matern52-assigned probes that are significantly linear

venn.diagram(x=list(Linear=rownames(lin), Nonlinear=rownames(nonlin.probes)), category.names=c('Linear dDMPs', 'Nonlinear sites'), filename=paste0(plotPath,"Venn_linear_nonlinear.png"), output=TRUE, cat.pos=181)


#4. Plot linear DMPs that are actually nonlinear ================================================================================
load(paste0(methylationPath, "fetalBulk_EX3_23pcw_n91.rdat"))
load(paste0(GPPath, "gpmeanvar.rdata")) 

plot.example <- function(probe, betas, timescale, linear_pval, gene, chr, res, ribbon.factor=1){
  
  i <- which(rownames(res)==probe)
  df <- data.frame(Age=pheno$PCW, Betas=(betas[probe,])*100)
  
  #1. prepare ribbon
  mn <- gp_mean[probe,]
  std <- sqrt(gp_var[probe,])
  upper <- mn + ribbon.factor*std
  lower <- mn - ribbon.factor*std
  df2 <- data.frame(Ages=as.numeric(colnames(gp_mean)),Mean=(mn*100), Upper=(upper*100), Lower=(lower*100))
  if(any(df2$Lower<0 | df2$Upper>100)){
	df2[which(df2$Lower<0),'Lower'] <- 0
	df2[which(df2$Upper>100),'Upper'] <- 100
  }
  
  #2. create the title using probe ID, gene name (if applicable) and chr
  if(gene!=''){                                        # if gene is present...
	ttl <- paste0(probe, ' - ', gene, ' - CHR ',chr)   # combine probe ID, gene name and chromosome
  }else{                                               # if there's no gene name...
	ttl <- paste0(probe, ' - CHR ', chr)               # just probe ID and chr
  }
  
  #3. plot
  ggplot()+
    geom_point(data=df, aes(x=Age, y=Betas), size=3)+
    geom_line(data=df2, aes(x=Ages, y=Mean), size=1.2, colour="#3366FF")+
    geom_ribbon(data=df2, aes(x=Ages, ymin=Lower, ymax=Upper), fill="grey60", alpha=0.4)+
    xlab('Age (pcw)')+
    ylab('DNA methylation (%)')+
    ylim(0,100)+
    ggtitle(paste0(ttl, '\nTimescale ', timescale, '       Linear p = ', linear_pval))+
    theme_minimal()+
    theme(axis.text=element_text(size=19), axis.title=element_text(size=20))+
    theme(plot.title = element_text(size=21))
}

df.sub <- df[df$Group=='Linear',] # rows are duplicated from the way df was created. Only keep one version of rownames
df.sub <- df.sub[df.sub$P.Age<1e-20,] # subset of highly significant linear DMPs #2714
df.sub <- df.sub[df.sub$mat52_ell>=10 & df.sub$mat52_ell<20,] # limit timescale to small but biologically meaningful range #214
df.ord <- df.sub[order(df.sub$mat52_ell),] # order by timescale: smallest to largest
df.hypo <- df.ord[df.ord$Beta.Age<0,]
df.hyper <- df.ord[df.ord$Beta.Age>0,]
df.plot <- rbind(df.hypo[1:10,], df.hyper[1:10,])

pdf(paste0(plotPath, "linearDMP_nonlinear_examples_HypoHyper.pdf"))

for(probe in rownames(df.plot)[1:20]){
	timescale <- signif(df.plot[probe, 'mat52_ell'], digits=3)
	linear_pval <- signif(df.plot[probe, 'P.Age'], digits=3)
	chr <- lin[probe,'CHR']
	gene <- lin[probe,'Gene']
	show(plot.example(probe=probe, betas=betas, timescale=timescale, linear_pval=linear_pval, gene=gene, chr=chr, res=df.plot, ribbon.factor=1))
}

dev.off()