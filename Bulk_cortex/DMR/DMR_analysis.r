

## Perform DMR analysis on bulk fetal Age EWAS results ##


#1. Load data ===================================================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))


#2. Run DMR script ==============================================================================================================
source("DMR_functions.r")

res <- readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds"))
run.dmrff(betas=betas, res=res, ES.col='Beta.Age', SE.col='SE.Age', P.col='P.Age', sig.thresh=9e-8, 
		dmr.filename="ageReg_fetalBrain_EX3_23pcw_DMR", 
		stats.filename="ageReg_fetalBrain_EX3_23pcw_statsForDMR")
		

#3. Annotate DMR results ========================================================================================================
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Homo.sapiens)
require(AnnotationHub)
library(pbapply) # progress bar

pad <- 30000 # base pairs to plot either side of DMR on x-axis

hub <- AnnotationHub()
query(hub, c("cpg","hg19"))
cpgs <- as.data.frame(hub[["AH5086"]])
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
mat <- as.data.frame(transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL")))


grange.anno <- function(x){
	
	# extract chr, start and end position of DMR
	chr <- as.character(unlist(x[which(colnames(dmrs)=='chr')]))
	start <- as.numeric(x[which(colnames(dmrs)=='start')])
	end <- as.numeric(x[which(colnames(dmrs)=='end')])
	range <- c(start,end)
	
	# identify CpGs within this region
	chr2 <- paste("chr", chr, sep="")
	cpg <- cpgs[which(cpgs$seqnames == chr2 & cpgs$start >=(min(range)-pad) & 
	cpgs$start <= max(range)+pad),]
	
	# match to reference
	gr <- GRanges(seqnames = chr2, ranges = IRanges(start = min(range)-pad, 
	end = max(range)+pad))
	trans <- subsetByOverlaps(transcripts(txdb), gr)

	# add gene info
	m <- as.matrix(mat[mat$TXNAME %in% as.data.frame(trans)[,"tx_name"],"SYMBOL"])
	colnames(m) <- "SYMBOL"
	trans <- cbind(as.data.frame(trans),m)
	trans <- trans[is.na(trans$SYMBOL) == "FALSE",]
	
	# collapse genes
	genes <- as.character(unlist(trans$SYMBOL))
	gene <- paste(unique(unlist(strsplit(genes, "\\;"))), collapse = ";")
	return(gene)
}

dmr <- readRDS("ageReg_fetalBrain_EX3_23pcw_DMR.rds") 	#43884
dmrs <- dmr[which(dmr$p.adjust < 0.05 & dmr$n > 2),]	# 1356
dmrs <- dmrs[order(dmrs$p.adjust),]

dmr.genes <- pbapply(dmrs, 1, grange.anno)
dmrs$Gene <- dmr.genes
saveRDS(dmrs, file="ageReg_dmrGenes_GRangesAnnot.rds")

# calculate stats
dmrs <- readRDS("ageReg_dmrGenes_GRangesAnnot.rds") #1356
dmrs$range <- dmrs$end - dmrs$start
mean(dmrs$range) # 368.6475
res.dir <- sign(dmrs$estimate) # -1: 531, 1: 825
binom.test(length(which(res.dir==1)), length(res.dir), alternative='greater')$p.value #p-value = 6.869432e-16. Hypermethylated DMRs are enriched.
perc.hyper <- (length(which(res.dir==1))/length(res.dir))*100 #60.84071
length(unique(unlist(strsplit(as.character(dmrs$Gene), ";")))) #2292