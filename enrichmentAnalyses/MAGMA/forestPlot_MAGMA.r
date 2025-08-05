

## Plot MAGMA results ##

library(ggplot2)

# bulk fetal MAGMA results
bulkRes_ASD <- read.table(paste0(bulkPath, "ASD2019_ageReg_fetalBrain_MAGMA.gsa.out"), header=T)
bulkRes_SCZ <- read.table(paste0(bulkPath, "SCZ_ageReg_fetalBrain_MAGMA.gsa.out"), header=T)

# FANS fetal MAGMA results
fansRes_ASD_N <- read.table(paste0(fansPath, "ASD2019_fetalNeuronal_MAGMA.gsa.out"), header=T)
fansRes_ASD_NN <- read.table(paste0(fansPath, "ASD2019_fetalNon.neuronal_MAGMA.gsa.out"), header=T)
fansRes_SCZ_N <- read.table(paste0(fansPath, "SCZ_fetalNeuronal_MAGMA.gsa.out"), header=T)
fansRes_SCZ_NN <- read.table(paste0(fansPath, "SCZ_fetalNon.neuronal_MAGMA.gsa.out"), header=T)

# combine results for plotting
ASD <- rbind(bulkRes_ASD, fansRes_ASD_N, fansRes_ASD_NN)
SCZ <- rbind(bulkRes_SCZ, fansRes_SCZ_N, fansRes_SCZ_NN)

plot_forest <- function(res, ylab){
	p <- res$P
	b <- res$BETA
	se <- res$SE	
	
	df <- data.frame(P=p,Beta=b,SE=se)
	df$Error <- 1.96 * df$SE # 95% confidence interval
	df$P <- signif(as.numeric(df$P), digits=3)
	df$Upper <- df$Beta + df$Error
	df$Lower <- df$Beta - df$Error
	df$Label <- ylab
	df$Label <- factor(df$Label, levels=rev(unique(df$Label))) # reverses the factor level ordering for labels after coord_flip()
	
	# make xlim (originally ylim) symmetrical
	vals <- c(min(df$Lower),max(df$Upper))
	xlim <- max(abs(vals))
	
	# positioning of p-value text
	ps <- df$P
	y.p.pos <- nrow(df):1
	x.p.pos <- df$Beta

	forest <- ggplot(data=df, aes(x=Label, y=Beta, ymin=Lower, ymax=Upper)) +
			geom_pointrange()+			
			geom_hline(yintercept=0, lty=2)+  # add a dotted line at x=1 after flip
			coord_flip()+  # flip coordinates (puts labels on y axis)
			xlab("")+ ylab("Effect size")+
			ylim(-xlim,xlim)+
			theme_bw()+  # white background
			annotate("text", label=paste('p =', ps), y=x.p.pos, x=y.p.pos+0.4, size=4, colour="black")+
			theme(axis.text=element_text(size=15), axis.title=element_text(size=18))+
			theme(legend.position = "none")+
			theme(aspect.ratio = 1/3)	
	show(forest)
}

ylab <- c('Bulk', 'Neuronal', 'Non-neuronal')

setwd(plotPath)

pdf("MAGMA_fetalAgeReg_ASD2019.pdf", width=8, height=3)
plot_forest(ASD, ylab)
dev.off()

pdf("MAGMA_fetalAgeReg_SCZ2022.pdf", width=8, height=3)
plot_forest(SCZ, ylab)
dev.off()