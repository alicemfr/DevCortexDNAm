
## Heatmap of ATAC-seq peak enrichment for each nonlinear module ##

library(reshape2) # for melt()
library(ggplot2)
library(RColorBrewer)


#1. Load enrichment results - nonlinear =========================================================================================
modules <- c('All', 'turquoise', 'blue', 'brown', 'green', 'red', 'yellow')

stats <- read.csv(paste0("Enrichment_Peaks_nonlinear_",modules[1],".csv"), row.names=1)
cells <- rownames(stats)
# extract significance
ps <- stats[,'P.adj']
ps <- as.matrix(ps)
colnames(ps) <- paste(modules[1],'P.adj',sep='.')
rownames(ps) <- cells
# extract effect sizes
stats <- stats[,c('Estimate')]
stats <- as.matrix(stats)
colnames(stats) <- paste(modules[1],'Estimate',sep='.')
rownames(stats) <- cells

for(m in modules[-1]){
	print(m)
	stats.m <- read.csv(paste0("Enrichment_Peaks_nonlinear_",m,".csv"), row.names=1)
	# extract significance
	p <- stats.m[,'P.adj']
	p <- as.matrix(p)
	colnames(p) <- paste(m,'P.adj',sep='.')
	ps <- cbind(ps, p)
	# extract effect sizes
	stats.m <- stats.m[,c('Estimate')]
	stats.m <- as.matrix(stats.m)
	colnames(stats.m) <- paste(m,'Estimate',sep='.')
	stats <- cbind(stats, stats.m)
}

# assign effect size of zero if non-significant
sig <- ps<5e-2
stats[!sig] <- 0

stats <- data.frame(stats)
stats$Cell_type <- rownames(stats)



#2. Load enrichment results - linear ============================================================================================
stats.l <- read.csv("fetalBulk_peakEnrichment_logRegStats.csv", row.names=1)
cells <- rownames(stats.l)
# extract significance
ps <- stats.l[,'P.adj']
ps <- as.matrix(ps)
colnames(ps) <- paste('linear','P.adj',sep='.')
rownames(ps) <- cells
# extract effect sizes
stats.l <- stats.l[,c('Estimate')]
stats.l <- as.matrix(stats.l)
colnames(stats.l) <- paste('linear','Estimate',sep='.')
rownames(stats.l) <- cells
# assign effect size of zero if non-significant
sig <- ps<5e-2
stats.l[!sig] <- 0

stats.l <- data.frame(stats.l)
stats.l$Cell_type <- rownames(stats.l)


#3. Combined linear and nonlinear results =======================================================================================
identical(stats.l$Cell_type, stats$Cell_type)
stats.l$Cell_type <- NULL
stats <- cbind(stats.l,stats)

stat.names <- c('linear',modules)
stat.labels <- c('All linear sites','All nonlinear sites', modules[-which(modules=='All')])


#4. Heatmap =====================================================================================================================
m <- melt(stats)
levels(m$variable) <- stat.names
m$Cell_type <- factor(m$Cell_type, levels=rev(unique(m$Cell_type))) # reverses the factor level ordering for labels after coord_fixed()
p <- ggplot(m, aes(y=Cell_type, x=variable, fill=value))+
	geom_tile(color = "white",lwd = 1.5,linetype = 1)+
	ylab("")+
	xlab("")+
	scale_fill_gradient2(low="steelblue3", mid='white', high="firebrick")+#, limits=c(2,9.5))+#, limits=c(0.75,7))+
	geom_vline(xintercept=c(0.5,2.5), size=0.2, colour='black')+
	guides(fill=guide_colourbar(title="Effect\nsize"))+
	theme_minimal()+
	scale_x_discrete(breaks=stat.names,labels=stat.labels)+
	theme(	plot.margin=unit(c(0.5,0.5,1.5,0.5),"cm"),
			axis.text.x=element_text(angle=45, vjust=0.5, hjust=1, size=12, margin=margin(-30,0,0,0)),
			axis.text=element_text(size=10),
			axis.title=element_text(size=18),
			legend.justification = "top"
		)
pdf(paste0("Heatmap_Peaks_nonlinear.pdf"), width=8, height=10)
show(p)
dev.off()