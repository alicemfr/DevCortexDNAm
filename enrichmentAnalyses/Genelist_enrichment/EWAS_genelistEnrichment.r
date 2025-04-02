

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

geneEnrichment <- function(inputFile, resultCol, epicAnnotGeneList, outputFile, glists, glistFile, comboList=TRUE){

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

	#1. Number of probes per gene
	gene_size <- table(unlist(strsplit(res$Gene, "\\;")))
	print(length(gene_size))
	
	#2. Number of DMPs annotated to gene
	res$DMP <- as.numeric(res[,resultCol]<9e-8)
	gene_dmp <- table(unlist(strsplit(res$Gene[res$DMP==1], "\\;")))
	gene_test <- gene_dmp[match(names(gene_size), names(gene_dmp))]
	names(gene_test) <- names(gene_size)
	gene_test[is.na(gene_test)]<-0
	print(length(gene_test))
	
	mat <- matrix(NA,nrow=length(glists),ncol=10)
	colnames(mat) <- c("GeneList", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "Estimate.GenesinTestList", "OR.GenesinTestList", "P.GenesinTestList", "Estimate.GenesSize", "P.GeneSize")
	
	for(i in 1:length(glists)){
	  glist <- glists[i]
	  print(glist)
	  genes <- glistFile[glistFile$GeneList==glist,'Gene']
	  
	  #3. Is gene in genelist: 0/1
	  pathway <- as.numeric(names(gene_size) %in% genes)
	  print(length(pathway))
	  
	  mat[i,] <- apply_tests(glist, pathway, gene_test, gene_size)
	  
	}
	
	write.csv(mat, outputFile)

}

geneEnrichment(inputFile, resultCol, epicAnnotGeneList, outputFile, glists=c('SFARI','SCHEMA'), glistFile)