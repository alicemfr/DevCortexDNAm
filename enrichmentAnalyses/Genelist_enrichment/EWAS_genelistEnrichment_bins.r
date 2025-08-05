

## Test enrichment of gene lists within EWAS results ##


epicAnnotGeneList <- read.csv(paste0(refPath, "EPIC_annot_SFARI_SCHEMA.csv"), row.names=1)
sfari <- unique(epicAnnotGeneList$SFARI.Gene[-which(is.na(epicAnnotGeneList$SFARI.Gene)|epicAnnotGeneList$SFARI.Gene=='')])
schema <- unique(epicAnnotGeneList$SCHEMA.Gene[-which(is.na(epicAnnotGeneList$SCHEMA.Gene)|epicAnnotGeneList$SCHEMA.Gene=='')])
glistFile <- data.frame(Gene=c(sfari, schema), GeneList=c(rep('SFARI',length(sfari)),rep('SCHEMA',length(schema))))
glistFile <- rbind(glistFile, data.frame(Gene=unique(glistFile$Gene), GeneList=rep('Combined', length(unique(glistFile$Gene)))))

#===================================================================================================================================

apply_tests<-function(glist, pathway, gene_test, gene_size){

	model <- glm(pathway ~ gene_test + gene_size, family="binomial")
	p1 <- summary(model)$coefficients['gene_test','Pr(>|t|)']
	p2 <- summary(model)$coefficients['gene_size','Pr(>|t|)']
	b1 <- summary(model)$coefficients[c("gene_test"),'Estimate']
	b2 <- summary(model)$coefficients[c("gene_size"),'Estimate']
	or <- exp(b1)
	
	# divide the p-value by 2 to convert two-sided to one-sided test
	if(or > 1){															                        
		p1 <- p1/2
	}else{
		p1 <- 1-(p1/2)
	}
	
	return(c(
	glist,                                                    # the genelist we've tested
	sum(gene_size[which(pathway == 1)]),                      # total number of probes annotated to genelist genes
	length(which(pathway == 1)),                              # number of genes in genelist
	sum(gene_test[which(pathway == 1  & gene_test >= 1)]),    # total number of DMPs annotated to genelist genes
	length(which(pathway == 1 & gene_test >= 1)),             # number of genelist genes that contain a DMP
	b1,                                                       # effect size for gene_test
	or,                                                       # odds ratio for gene_test
	p1,                                                       # p-value for gene_test
	b2,                                                       # effect size for gene_size
	p2                                                        # p-value for gene_size
	))
	
}

geneEnrichment <- function(inputFile, resultCol, epicAnnotGeneList, outputFile, glists, glistFile, comboList=TRUE, bin.width=10000){

	print(paste('Input file:', inputFile))
	
	if(grepl(".rds", inputFile)){
	  res <- data.frame(readRDS(inputFile))
	}
	if(grepl(".csv", inputFile)){
	  library(data.table)
	  res <- data.frame(fread(inputFile, data.table=F), row.names=1)
	}
	
	print(dim(res))

	# annotate EWAS results with EPIC manifest that has been pre-annotated with whether probe is in genelist
	epicMan <- epicAnnotGeneList[match(rownames(res), epicAnnotGeneList$IlmnID),]
	res <- cbind(res, as.data.frame(epicMan))
	if(any(res$CHR %in% 'Y')){ res <- res[-which(res$CHR=='Y'),] }
	print(dim(res))
	
	if(comboList){
		res_glist_cols <- res[,colnames(res) %in% glists] # glist cols are 0 or 1 for whether probe is in gene list gene
		glist_sum <- apply(res_glist_cols, 1, sum) # sum across glist cols for each probe
		res$Combined <- as.numeric(glist_sum>0) # probes where sum>0 indicates probe is present in at least one gene list
		glists <- c(glists, "Combined") # add Combined column as a gene list to test
	}

	#1. Number of probes per gene -----------------------------------------------------------------------------------------------
	gene_size <- table(unlist(strsplit(res$Gene, "\\;")))
	print(length(gene_size))
	
	#2. Number of DMPs annotated to gene ----------------------------------------------------------------------------------------
	res$DMP <- as.numeric(res[,resultCol]<9e-8)
	res <- res[order(res[,resultCol]),] # order p-values smallest to largest
	
	# Identify how many bins to use
	all.dmp.indx <- rownames(res)[which(res$DMP==1)] # DMP probe IDs
	n.bins <- ceiling(length(which(res$DMP==1))/bin.width) # round up to next integer
	w <- bin.width
	
	stats <- matrix(NA, nrow=1,ncol=11)
	colnames(stats) <- c("Bin","GeneList", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "Estimate.GenesinTestList", "OR.GenesinTestList", "P.GenesinTestList", "Estimate.GenesSize", "P.GeneSize")
	
	for(n in 1:n.bins){
		print(n)
		
		# Determine start and end indeces for bins
		# e.g. n=2, start=1, end=20000
		start <- 1
		if(n==n.bins){
			end <- length(all.dmp.indx) # stop final bin where dmps end
		}else{
			end <- n*w
		}
		print(start)
		print(end)
		
		# DMPs within bin
		keep.dmp.indx <- all.dmp.indx[start:end]
		if(n!=n.bins & length(keep.dmp.indx)!=(n*w)){
			stop(paste("keep.dmp.indx length =", length(keep.dmp.indx), "but expected", (n*w)))
		}
		
		# DMPs not in bin
		rem.dmp.indx <- all.dmp.indx[all.dmp.indx %ni% keep.dmp.indx]
		if(n!=n.bins & length(rem.dmp.indx)!=(length(all.dmp.indx)-(n*w))){
			stop(paste("rem.dmp.indx length =", length(rem.dmp.indx), "but expected", (length(all.dmp.indx)-(n*w))))
		}
		
		# filter results for DMPs to keep
		if(length(rem.dmp.indx!=0)){ # zero when n==n.bins
			res.bin <- res[-which(rownames(res) %in% rem.dmp.indx),]
		}else{
			res.bin <- res
		}
		
		gene_dmp <- table(unlist(strsplit(res.bin$Gene[res.bin$DMP==1], "\\;")))
		gene_test <- gene_dmp[match(names(gene_size), names(gene_dmp))]
		names(gene_test) <- names(gene_size)
		gene_test[is.na(gene_test)]<-0
	
		mat <- matrix(NA,nrow=length(glists),ncol=11)
		mat[,1] <- rep(n,length(glists))
	
		for(i in 1:length(glists)){
		  glist <- glists[i]
		  genes <- glistFile[glistFile$GeneList==glist,'Gene']
		  
		  #3. Is gene in genelist: 0/1 ------------------------------------------------------------------------------------------
		  pathway <- as.numeric(names(gene_size) %in% genes)
		  
		  mat[i,2:ncol(mat)] <- apply_tests(glist, pathway, gene_test, gene_size)
		}
		
		stats <- rbind(stats, mat)
	
	}
	
	stats <- stats[-1,] # remove first NA row
	write.csv(stats, outputFile)

}

# Run enrichment ================================================================================================================

# bulk
inputFile <- "ageReg_fetalBrain_EX3_23pcw.rds"
resultCol <- 'P.Age'
geneEnrichment(inputFile, resultCol, epicAnnotGeneList, outputFilename_bulk, glists=c('SFARI','SCHEMA'), glistFile, bin.width=10000)

# FANS - neuronal
inputFile <- "FACS_AgeCellSpecific_EWAS_fetal_adult_anno.csv"
resultCol <- 'Fetal.Neuronal.P.Age'
geneEnrichment(inputFile, resultCol, epicAnnotGeneList, outputFilename_SATB2, glists=c('SFARI','SCHEMA'), glistFile, bin.width=400)


# Plot ==========================================================================================================================

combo <- read.csv(fullRes) # original enrichment results with full DMP list (i.e. not binned)
bins <- read.csv(binnedRes), row.names=1) # binned enrichment results from above

n.tests <- nrow(combo) # 3: Combined, SFARI, SCHEMA

plotCombovsBins <- function(genelist, bin.breaks, bin.labels){
	df <- data.frame(X=1:(max(bins$Bin)+1), OR=c(combo[combo$GeneList==genelist,'OR.GenesinTestList'],bins[bins$GeneList==genelist,'OR.GenesinTestList']), Group=c('Combined',rep('Binned',max(bins$Bin))))

	p <- ggplot(df, aes(x=X, y=OR))+
		geom_line(aes(color=Group), linewidth=1.2)+
		geom_point(size=2)+
		xlab("")+
		ylab("Odds ratio")+
		#ylim(1,1.2)+ # for bulk
		ylim(1,3.5)+ # for FANS
		theme_minimal()+
		theme(axis.text=element_text(size=19),
			axis.title=element_text(size=19),
			plot.title=element_text(size=23),
			legend.position='none',
			plot.margin=unit(c(0.5,0.5,1,0.5),"cm"),
			axis.text.x=element_text(angle=45, vjust=0.5, hjust=1, size=12, margin=margin(-20,0,0,0)),
			axis.text.y=element_text(size=15),
			axis.title.y=element_text(margin=margin(0,10,0,0)))+
		scale_color_manual(values=c("lightgrey","black"))+
		scale_x_continuous(breaks=bin.breaks,labels=bin.labels)+
		geom_hline(yintercept=1, linetype="dashed", color="black")
	show(p)
}

# bulk
for(g in c('SFARI','SCHEMA','Combined')){
	bin.labels <- c("all dDMPs","1 - 10,000","1 - 20,000","1 - 30,000","1 - 40,000", "1 - 50,000","1 - 50,913")
	bin.breaks <- 1:length(bin.labels)
	pdf(paste0("ageReg_fetalBulk_", g, "_logP_bins.pdf"), width=5,height=5)
	plotCombovsBins(g, bin.breaks, bin.labels)
	dev.off()
}

# FANS - neuronal
for(g in c('SFARI','SCHEMA','Combined')){
	bin.labels <- c("all dDMPs","1 - 200","1 - 400","1 - 600","1 - 800", "1 - 1,000","1 - 1,200","1 - 1,400","1 - 1,600","1 - 1,800","1 - 1,872")
	bin.breaks <- 1:length(bin.labels)
	pdf(paste0("FANS_Neuronal_", g, "_logP_bins.pdf"), width=5,height=5)
	plotCombovsBins(g, bin.breaks, bin.labels)
	dev.off()
}