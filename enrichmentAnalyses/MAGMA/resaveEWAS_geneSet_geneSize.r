

## Create gene set file and gene size file from EWAS results table for use with MAGMA ##


library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(data.table)


#1. Load results ================================================================================================================
if(grepl(".rds", resFile)){
	res <- readRDS(resFile)
}
if(grepl(".csv", resFile)){
	res <- data.frame(fread(resFile, data.table=F), row.names=1)
}


#2. (optional) Remove non-variable probes =======================================================================================
if(filterNonVarProbes){
	nonvar <- readRDS(nonVarFile)
	res <- res[-which(rownames(res) %in% rownames(nonvar)),]
}


#3. Remove probes not annotated to a gene =======================================================================================
if(!any(grepl('Gene',colnames(res)))){
	if(!any(grepl('epicManifest', ls()))){
		library(data.table)
		epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
	}
	epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c("CHR","UCSC_RefGene_Name")]
	res <- cbind(res, as.data.frame(epicMan))
	uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
	res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))
}
res <- res[-which(res$Gene=='' | is.na(res$Gene)),]	


#4. Create background file (all genes) ==========================================================================================
geneSize_bg <- function(res){
	gene_size <- table(unlist(strsplit(res$Gene, "\\;")))
	gene_name <- names(gene_size)
	gene_size <- as.numeric(gene_size)
	gene_set <- rep('Background',length(gene_size))
	
	res.background <- data.frame(GeneSet=gene_set, GeneName=gene_name, GeneSize=gene_size)
	
	gene_entrez <- select(org.Hs.eg.db, keys=as.character(res.background$GeneName), columns="ENTREZID", keytype="SYMBOL")
	identical(unique(as.character(res.background$GeneName)), unique(as.character(gene_entrez$SYMBOL)))
	gene_entrez_match <- gene_entrez[match(as.character(res.background$GeneName),gene_entrez$SYMBOL),]
	identical(as.character(res.background$GeneName), as.character(gene_entrez_match$SYMBOL))
	res.background$EntrezID <- gene_entrez_match$ENTREZID
	
	return(res.background)
}
res.geneSize_bg <- geneSize_bg(res)


#5. Filter for significant probes ===============================================================================================
sigProbes <- function(res, Pcol, sigLevel=9e-8, EScol=NULL, Dir=NULL){
	sig <- res[which(res[,Pcol]<sigLevel),]
	m <- sig
	
	if(!is.null(EScol)){
		if(Dir=='hypo'){
			m <- sig[which(sig[,EScol]<0),]
		}
		if(Dir=='hyper'){
			m <- sig[which(sig[,EScol]>0),]
		}
	}
	
	return(m)
}

res.sig <- sigProbes(res, Pcol=Pcol)

if(splitDir){
	res.sig.hypo <- sigProbes(res, Pcol=Pcol, EScol=EScol, Dir='hypo')
	res.sig.hyper <- sigProbes(res, Pcol=Pcol, EScol=EScol, Dir='hyper')
}


#6. Create Gene Set and Gene Size files =========================================================================================
sigRes.geneSize <- function(sig, geneCol='Gene', geneSetName, background){

	gene_name <- unique(unlist(strsplit(sig[,geneCol], "\\;")))
	gene_set <- rep(geneSetName,length(gene_name))
	df <- data.frame(GeneSet=gene_set, GeneName=gene_name)

	# use the gene sizes calculated from the background set
	bg <- background[which(as.character(background$GeneName) %in% as.character(df$GeneName)),]
	bg <- bg[match(as.character(df$GeneName), as.character(bg$GeneName)),]
	identical(as.character(bg$GeneName), as.character(df$GeneName))
	
	df$GeneSize <- bg$GeneSize
	df$EntrezID <- bg$EntrezID

	# gene set
	gene.set <- df[,c('GeneSet','EntrezID')]
	gene.set <- gene.set[-which(is.na(gene.set$EntrezID)),]
	
	# gene size
	gene.size <- df[,c('EntrezID','GeneSize')]
	gene.size <- gene.size %>% distinct(EntrezID, .keep_all = TRUE)
	gene.size <- gene.size[-which(is.na(gene.size$EntrezID)),]

	return(list(GeneSet=gene.set, GeneSize=gene.size))

}

res.sig.df <- sigRes.geneSize(sig=res.sig, geneSetName=geneSetName, background=res.geneSize_bg)

if(splitDir){
	res.sig.hypo.df <- sigRes.geneSize(sig=res.sig.hypo, geneSetName=paste0(geneSetName,".hypo"), background=res.geneSize_bg)
	res.sig.hyper.df  <- sigRes.geneSize(sig=res.sig.hyper, geneSetName=paste0(geneSetName,".hyper"), background=res.geneSize_bg)
}


#7. Combine results =============================================================================================================
if(splitDir){
	res.gene.set <- rbind(res.sig.df$GeneSet, res.sig.hypo.df$GeneSet, res.sig.hyper.df$GeneSet)
}else{
	res.gene.set <- res.sig.df$GeneSet
}

# add background
geneSet_bg <- function(res.full, geneCol='Gene', geneSetName){
	
	gene_name <- unique(unlist(strsplit(res.full[,geneCol], "\\;")))
	gene_set <- rep(geneSetName,length(gene_name))
	df <- data.frame(GeneSet=gene_set, GeneName=gene_name)
	
	gene_entrez <- select(org.Hs.eg.db, keys=as.character(df$GeneName), columns="ENTREZID", keytype="SYMBOL")
	identical(unique(as.character(df$GeneName)), unique(as.character(gene_entrez$SYMBOL)))
	gene_entrez_match <- gene_entrez[match(as.character(df$GeneName),gene_entrez$SYMBOL),]
	identical(as.character(df$GeneName), as.character(gene_entrez_match$SYMBOL))
	df$EntrezID <- gene_entrez_match$ENTREZID
	if(length(which(is.na(df$EntrezID)))!=0){
		df <- df[-which(is.na(df$EntrezID)),]
	}
	
	return(df)
}

res_bg_geneSet <- geneSet_bg(res.full=res, geneSetName="Background")

# checking that Background set contains all genes across each gene set
'%ni%' <- Negate('%in%')
length(which(res.gene.set$EntrezID %ni% res_bg_geneSet$EntrezID)) # should be 0
res.gene.set.all <- rbind(res.gene.set, res_bg_geneSet[,c('GeneSet','EntrezID')])


#8. Save results ================================================================================================================

# Save gene set files
write.table(res.gene.set.all, file=paste0(outPath,outFilePrefix,"_geneSet.txt"), row.names=F, col.names=F, sep='\t', quote=FALSE)

# Save gene size files
res_bg_geneSize <- res.geneSize_bg[-which(is.na(res.geneSize_bg$EntrezID)),]
res_bg_geneSize <- res_bg_geneSize %>% distinct(EntrezID, .keep_all = TRUE)

write.table(res_bg_geneSize[,c('EntrezID','GeneSize')], file=paste0(outPath,outFilePrefix,"_geneSize.txt"), row.names=F, sep='\t', quote=FALSE)