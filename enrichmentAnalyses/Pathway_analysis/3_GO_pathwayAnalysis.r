

## Logistic regression GO pathway analysis ##
 # Created by Dr Eilis Hannon
 # Edited by Alice Franklin

args <- commandArgs(trailingOnly=TRUE)

inputFile <- as.character(args[1])
resultCol <- as.character(args[2])
outputPath <- as.character(args[3])
test <- as.character(args[4])
top10 <- as.factor(args[5])
if(top10==TRUE){
	test <- paste(test, 'top10perc', sep='_')
}
splitDirection <- as.factor(args[6])
if(splitDirection==TRUE){
	directionCol <- as.character(args[7])
	direction <- as.character(args[8])
	test <- paste(test, direction, sep='_')
}
outputFile1 <- paste0(outputPath,test,"_pathwayAnalysis.csv")
outputFile2 <- paste0(outputPath,test,"_pathwayAnalysis_reducedTerms.csv")

print(paste0('inputFile=',inputFile))
print(paste0('resultCol=',resultCol))
print(paste0('outputPath=',outputPath))
print(paste0('top10=',top10))
print(paste0('splitDirection=',splitDirection))
print(paste0('outputFile (all terms) =',outputFile1))
print(paste0('outputFile (reduced terms)=',outputFile2))


#===================================================================================================================================

uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
	

apply_tests<-function(each){

	## logistic regres.metasion approach
	pathway<-vector(length = nrow(bg_gene_go))
	pathway[grep(each, bg_gene_go[,"GOTerms"])]<-1                      # positions where that GO term comes up
	model<-glm(pathway ~ gene_test + gene_size)                         # testing the relationship between the genes that share that GO term and the presence of significant probes across genes, taking into account gene size
	p1<-summary(model)$coefficients[c("gene_test"),c(4)]                # gene_test p-value
	p2<-summary(model)$coefficients[c("gene_size"),c(4)]                # gene_size p-value
	or<-exp(summary(model)$coefficients[c("gene_test"),c(1)])           # odds ratio: the Estimate for gene_test is the log odds. Holding gene_size constant, the odds of a gene containing a sig probe is exp(Estimate.gene_test).
	beta<-summary(model)$coefficients[c("gene_size"),c(1)]              # gene_size Estimate
	if(or > 1){                                                         # divide the p-value by 2 to convert two-sided to one-sided test
		p1<-p1/2
	} else {
		p1<-1-(p1/2)
	}
		
	model<-glm(pathway ~ gene_test_ind + gene_size)                     # repeat for the binary variable gene_test_ind, where it doesn't matter how many sig probes there are
	p1.ind<-summary(model)$coefficients[c("gene_test_ind"),c(4)]
	p2.ind<-summary(model)$coefficients[c("gene_size"),c(4)]
	or.ind<-exp(summary(model)$coefficients[c("gene_test_ind"),c(1)])
	beta.ind<-summary(model)$coefficients[c("gene_size"),c(1)]
	if(or.ind > 1){
		p1.ind<-p1.ind/2
	} else {
		p1.ind<-1-(p1.ind/2)
	}
		
	return(unlist(c(
	
	each, 																		# the GO term we've tested
	names[match(each, names[,"ID"]),c("Name","NameSpace")], 					# the GO info about that term: Name and NameSpace
	sum(gene_size[which(pathway == 1)]), 										# the total number of probes across all genes that share this GO term
	length(which(pathway == 1)), 												# the number of genes  that share this GP term
	sum(gene_test[which(pathway == 1  & gene_test >= 1)]), 						# the total number of probes across all *significant* genes that share this GO term
	length(which(pathway == 1 & gene_test >= 1)), 								# the number of *significant* genes that share this GO term
	p1,																			# the p-value for gene_test
	or,																			# the odds ratio for gene_test
	p2,																			# the p-value for gene_size
	beta,																		# the estimate for gene_size
	p1.ind,																		# the p-value for gene_test, when using the binary gene_test_ind variable
	or.ind,																		# the odds ratio for gene_test, when using the binary gene_test_ind variable
	p2.ind,																		# the p-value for gene_size, when using the binary gene_test_ind variable
	beta.ind,  																	# the estimate for gene_size, when using the binary gene_test_ind variable
	paste(bg_gene_go[which(pathway == 1 & gene_test >= 1),1], collapse = "|")	# the genes that have significant probes and share this GO term
	
	)))
}


#===================================================================================================================================

if(grepl(".rds",inputFile)){
	res <- data.frame(readRDS(inputFile))
}

if(grepl(".csv|.tsv",inputFile)){
	library(data.table)
	res <- data.frame(fread(inputFile,data.table=F),row.names=1)
}

print(dim(res))

if(splitDirection==TRUE){
print(paste('splitting into', direction, 'probes'))
	if(direction=='Hyper'){
		res <- res[which(res[,directionCol] > 0),]
	}
	if(direction=='Hypo'){
		res <- res[which(res[,directionCol] < 0),]
	}
print(dim(res))
}

#6. Load gene ontology (GO) reference files and match significant probes to genes ===============================================
if(!any(grepl('UCSC_RefGene_Name',colnames(res)))){
	if(!any(grepl('epicManifest', ls()))){
		library(data.table)
		epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
	}
	epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c("CHR","UCSC_RefGene_Name")]
	res <- cbind(res, as.data.frame(epicMan))
	res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))
}


gene_size<-table(unlist(strsplit(res$Gene, "\\;")))										# total number of probes annotated to each gene
bg_genes<-names(gene_size)																# names of all the unique genes in the results

if(top10==TRUE){
	res.sig <- res[which(res[,resultCol]<9e-8),]										# all significant probes
	res.ord <- res.sig[order(res.sig[,resultCol]),]										# order most to least significant
	totsig <- nrow(res.ord) 															# total number of sig probes
	perc10 <- 0.1*totsig																# 10% of number of sig
	top10prb <- res.ord[1:perc10,]
	gene_test<-table(unlist(strsplit(top10prb$Gene, "\\;")))
}else{
	gene_test<-table(unlist(strsplit(res$Gene[which(res[,resultCol]<9e-8)], "\\;")))    # for significant probes, total number of probes annotated to each gene
}


gene_go<-read.csv(paste0(refPath, "GOTerms_annot_genes_collapsed.csv"), stringsAsFactors = FALSE, row.names=1)
gene_go<-gene_go[which(gene_go[,"GOTerms"] != ""),]  #19731

bg_gene_go<-gene_go[match(intersect(bg_genes, gene_go[,"SYMBOL"]), gene_go[,"SYMBOL"]),]# match the gene go annotation df to the genes in our results 
termCount<-table(unlist(strsplit(as.character(bg_gene_go[,"GOTerms"]), "\\|")))			# total number of GO terms across all genes in results
terms<-names(termCount)[which(termCount > 9 & termCount <= 2000)]

names<-read.csv(paste0(refPath, "unique_GOTerm_names.csv"), stringsAsFactors = FALSE, row.names=1)

gene_size<-gene_size[bg_gene_go[,"SYMBOL"]]												# limit gene size vector to only the genes that have GO annotation
gene_test<-gene_test[bg_gene_go[,"SYMBOL"]]												# limit significant gene list vector to only the genes that have GO annotation
gene_test[is.na(gene_test)]<-0															# because we matched bg_gene_go genes to gene_test, there will be genes that are not in both, so it will create NAs

gene_test_ind<-gene_test
gene_test_ind[which(gene_test_ind > 0)]<-1												# 1 = this gene contains any number of significant probes, 0 = contains no sig probes

table(gene_test_ind)


#7. Run pathway analysis ========================================================================================================
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)

r<-foreach(i=1:length(terms), .combine=rbind, .export = c("gene_size", "gene_test_ind")) %dopar%{
	apply_tests(terms[i])
}
colnames(r)<-c("ID", "Name", "Type", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "P.GenesinTestList", "OR", "P.GeneSize", "Beta.GeneSize","P.GenesinTestList.1", "OR.1", "P.GeneSize.1", "Beta.GeneSize.1", "GenesinTestListAndPathway")

# filter to terms with between 10 and 2000 genes
r<-as.data.frame(r, stringsAsFactors = FALSE)
r<-r[which(as.numeric(r[,"nGenesinPathway"]) < 2000 & as.numeric(r[,"nGenesinPathway"]) > 9),]

r<-r[order(as.numeric(r$P.GenesinTestList)),]
r<-r[which(as.numeric(r$P.GenesinTestList) < 0.05/nrow(r)),-which(colnames(r) %in% c('P.GenesinTestList.1','OR.1','P.GeneSize.1','Beta.GeneSize.1'))]

print(paste0('Saving results to: ', outputFile1))
write.csv(r, outputFile1)


#8. Reduce GO terms =============================================================================================================
# Some GO terms will be significant only because they share genes involved in another pathway. Find the key pathways.

r <- read.csv(outputFile1, row.names=1)
r<-r[order(as.numeric(r$P.GenesinTestList)),]
r<-r[which(as.numeric(r$P.GenesinTestList) < 0.05/nrow(r)),] #60
print(dim(r))
print(head(r))

if(all(grepl('GO',rownames(r)))){
	r$GOTerm <- rownames(r)
}else{
	colnames(r)[which(colnames(r)=='ID')] <- 'GOTerm'
}
r$GOTerm <- as.character(r$GOTerm)
r$Name <- as.character(r$Name)
r$Type <- as.character(r$Type)

output<-c()
while(nrow(r)>0){

	if(class(r) != "character"){

	# for all terms repeat analysis controlling for most significant terms
	best_term<-vector(length =  nrow(bg_gene_go))
	best_term[grep(r[1,'GOTerm'], bg_gene_go[,"GOTerms"])]<-1 # find how many other genes share the top GO term
	merge.id<-c()
	merge.name<-c()
	remove<-c()
	
	for(j in 2:nrow(r)){
		print(j)
		pathway<-vector(length = nrow(bg_gene_go))
		pathway[grep(r[j,'GOTerm'], bg_gene_go[,"GOTerms"])]<-1
		model<-glm(pathway ~ gene_test_ind + gene_size + best_term)
		if(!is.na(summary(model)$coefficients["gene_test_ind", "Pr(>|t|)"]) & summary(model)$coefficients["gene_test_ind", "Pr(>|t|)"] > 0.05){
			merge.id<-append(merge.id, r[j,'GOTerm'])
			merge.name<-append(merge.name, r[j,'Name'])
			remove<-append(remove, j)
		} 
	
	}
	merge.id<-paste(unlist(as.character(merge.id)), collapse = "|")
	merge.name<-paste(unlist(as.character(merge.name)), collapse = "|")

	output<-rbind(output, c(r[1,], MergeID=merge.id, MergeName=as.character(merge.name)))
	
	r<-r[-c(1, remove),]
	} else {
	
	output<-rbind(output, c(r, "", ""))
	r<-NULL
}
}

write.csv(output, outputFile2)