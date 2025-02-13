

## Barplots of DMP distribution within genelist genes ##


library(viridis)

epicAnnotGeneList <- read.csv(paste0(refPath, "EPIC_annot_SFARI_SCHEMA.csv"), row.names=1)
glistFile <- read.csv(paste0(refPath, "GeneList_SFARI_SCHEMA.csv"), row.names=1)


#1. Load EWAS results ===========================================================================================================
if(grepl(".rds", inputFile)){
  res <- data.frame(readRDS(inputFile))
}
if(grepl(".csv", inputFile)){
  library(data.table)
  res <- data.frame(fread(inputFile, data.table=F), row.names=1)
}


#2. Calculate background number of significant genes ============================================================================
epicMan <- epicAnnotGeneList[match(rownames(res), epicAnnotGeneList$IlmnID),]
res <- cbind(res, as.data.frame(epicMan))
if(any(res$CHR %in% 'Y')){ res <- res[-which(res$CHR=='Y'),] }
res <- res[-which(res$Gene==''),] 

res.sig <- res[which(res[,resultCol]<9e-8),]
res.nonsig <- res[-which(rownames(res) %in% rownames(res.sig)),]

gene_size <- table(unlist(strsplit(res$Gene, "\\;")))
gene_dmp <- table(unlist(strsplit(res$Gene[which(res[,resultCol]<9e-8)], "\\;")))  # for significant probes, total number of probes annotated to each gene
background <- (length(gene_dmp)/length(gene_size))*100                             # percentage of genes with at least one sig probe vs all genes. 

gene_test <- gene_dmp[match(names(gene_size), names(gene_dmp))]
names(gene_test) <- names(gene_size)
gene_test[is.na(gene_test)]<-0

# Percentage DMPs
perc.sig <- (nrow(res.sig)/nrow(res))*100


#3. Calculate significant genes for gene lists ==================================================================================
glistFile <- rbind(glistFile, data.frame(Gene=glistFile$Gene, GeneList=rep('Combined', nrow(glistFile))))
res$Combined <- as.numeric(res$SFARI==1 | res$SCHEMA==1)
glists=c('SFARI','SCHEMA','Combined')

# Percentage DMPs
percDMPs <- perc.sig
for(i in 1:length(glists)){
  glist <- glists[i]
  res.pathway <- res[res[,glist]==1,]
  res.sig.pathway <- res.pathway[which(res.pathway[,resultCol]<9e-8),]
  perc.sig.pathway <- (nrow(res.sig.pathway)/nrow(res.pathway))*100
  percDMPs <- c(percDMPs,perc.sig.pathway)
}

pdf(plotFile_1, width=5, height=7)
par(mar=c(8, 6, 4.1, 2.1))
barplot(percDMPs, main="",ylab='Significant probes (%)', names.arg=c('All genes', glists), col=c('darkgray', mako(8)[c(4,7,3)]), ylim=c(0,(max(percDMPs)+1)), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2)
abline(h=percDMPs[1])
dev.off()

# Percentage of genes with at least one DMP
percs <- background
for(i in 1:length(glists)){
  glist <- glists[i]
  print(glist)
  genes <- glistFile[glistFile$GeneList==glist,'Gene']
  
  pathway <- as.numeric(names(gene_size) %in% genes)
  
  gene_size_pathway <- gene_size[pathway==1] # total number of probes annot to gene list genes
  gene_test_pathway <- gene_test[pathway==1] # number of DMPs annot to gene list genes
  
  pathway_sig_genes <- gene_test_pathway[gene_test_pathway!=0] # gene list genes that contain at least 1 DMP
  
  perc <- (length(pathway_sig_genes)/length(gene_test_pathway))*100
  percs <- c(percs, perc)
}

pdf(plotFile, width=5, height=7)
par(mar=c(8, 6, 4.1, 2.1))
barplot(percs, main="",ylab='% genes with at least\none significant probe', names.arg=c('All genes', glists), col=c('darkgray', mako(8)[c(4,7,3)]), ylim=c(0,max(percs)+2), las=2, cex.lab=1.5, cex.names=1.2, cex.axis=1.2)
abline(h=percs[1])
dev.off()