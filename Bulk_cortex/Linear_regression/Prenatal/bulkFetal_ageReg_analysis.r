

## Analysis of Age linear regression in bulk fetal cortex ## 


library(data.table)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(plyr) # for ddply
'%ni%' <- Negate('%in%') # not in

source(paste0(scriptsPath, "fetal_plotFunctions.r"))


#1. Load betas ===============================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))


#2. Load results =============================================================================================
res <- data.frame(readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")))


#3. Significant results ======================================================================================
res.sig <- res[which(res$P.Age<9e-8),] #50913
res.ord <- res.sig[order(res.sig$P.Age),] # smallest p-value to largest
saveRDS(res.ord, file=paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annot_sig.rds")) # save EWAS results for dDMPs

# percentage sig of all probes
(nrow(res.sig)/nrow(res))*100 #6.302627 %

# number of sig genes
unq.gn <- unique(unlist(strsplit(as.character(res.sig$Gene), ";"))) #11688

# n.probes annot to autosomal/X
res.chr <- res[-which(is.na(res$CHR)),]                     # 807607. 199 probes with NA CHR
res.aut <- res.chr[which(res.chr$CHR != 'X'),]              # 789981
res.X <- res.chr[which(res.chr$CHR == 'X'),]                # 17626

# n.sig.probes annot to autosomal/X
res.sig.chr <- res.sig[-which(is.na(res.sig$CHR)),]         # 50901. 12 probes with NA CHR
res.sig.aut <- res.sig.chr[which(res.sig.chr$CHR != 'X'),]  # 50329
perc.sig.aut <- nrow(res.sig.aut)/nrow(res.sig.chr)*100     # 98.87625 %
res.sig.X <- res.sig.chr[which(res.sig.chr$CHR == 'X'),]    # 572
perc.sig.X <- nrow(res.sig.X)/nrow(res.sig.chr)*100		    # 1.12375 %


#4. Split into hyper- and hypo-methylated ====================================================================
res.up <- res.ord[which(res.ord$Beta.Age>0),]			    # 22,133
res.down <- res.ord[which(res.ord$Beta.Age<0),]			    # 28,780
perc.res.up <- nrow(res.up)/nrow(res.ord)*100			    # 43.4722 %
perc.res.down <- nrow(res.down)/nrow(res.ord)*100		    # 56.5278 %

# sign test to see if enrichment of hypomethylated probes compared to hypermethylated
res.dir <- sign(res.sig$Beta.Age)
   # -1     1
# 28780 22133
binom.test(length(which(res.dir==-1)), length(res.dir), alternative='greater')$p.value # p = 1.612832e-191. Hypomethylation DMPs are more common than expected.


#5. Plot =====================================================================================================

# DNAm vs Age scatter plot
poi <- which(rownames(res.ord) %in% c(rownames(res.up)[1], rownames(res.down)[1])) # poi=position of interest. The index of the top res.up and top res.down in res.ord

myplots <- list()
for(j in 1:length(poi)){
	i <- poi[j]
	p <- age_scatter(res.ord, pheno, betas, i, age.column='PCW')
	myplots[[j]] <- p
}

pdf(paste0(plotPath, "ageReg_fetalHypoHyper.pdf"), width=7, height=7)
myplots
dev.off()


# Density of fetal DMP effect sizes
res.sig$Sign <- rep('Hypermethylated',nrow(res.sig))
res.sig$Sign[which(res.sig$Beta.Age<0)] <- 'Hypomethylated'
res.sig$Sign <- as.factor(res.sig$Sign)
df <- data.frame(Beta=abs(res.sig$Beta.Age)*100, Sign=res.sig$Sign) # convert DNAm to %
mu <- ddply(df, "Sign", summarise, grp.mean=mean(Beta)) # Hypermethylated 1.025026; Hypomethylated 1.532173
p <- t.test(Beta.Age ~ Sign, data=res.sig)$p.value # p = 0

pdf(paste0(plotPath, "Density_ageReg_dDMPDir.pdf"), width=8, height=8)
ggplot(df, aes(x=Beta, fill=Sign))+
	geom_density(alpha=0.3)+
	ylab("Density")+
	xlab("% change in DNA methylation per week")+
	geom_vline(data=mu, aes(xintercept=grp.mean, color=Sign), linetype="dashed")+
	theme_minimal()+
	theme(	axis.text=element_text(size=19),
			axis.title=element_text(size=22),
			plot.title = element_text(size=23),
			legend.position=c(0.75,0.9),
			legend.title = element_blank(),
			legend.text = element_text(size=16)
	)
dev.off()






















