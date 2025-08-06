

## Annotate EWAS results with ATAC-seq peaks ##

library(data.table)


# Run either...

# bulk fetal
EWASres <- paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")
outputFile <- "ageReg_fetalBrain_annot_allTissuePeaks.rds"
celltypes <- NULL
Pcol <- 'P.Age'

# FANS fetal
EWASres <- paste0(resultsPath, "FACS_AgeCellSpecific_EWAS_fetal_adult_anno.csv")
outputFile <- "FANS_AgeCellSpecific_annot_allTissuePeaks.rds"
celltypes <- c('Neuronal','Non.neuronal')
Pcols <- c()
for(ct in celltypes){ 
	Pcols <- c(Pcol, paste0('Fetal.',ct,'.P.Age'))
}


#1. Read EWAS results file ======================================================================================================
if(grepl('.rds',EWASres,)){ res <- readRDS(EWASres) }
if(grepl('.csv',EWASres,)){ 
	library(data.table)
	res <- data.frame(fread(EWASres, data.table=F), row.names=1)
}


#2. Filter probes ===============================================================================================================
if(length(which(is.na(res$CHR)))!=0){
	res <- res[-which(is.na(res$CHR)),]
}
res <- res[-which(res$CHR %in% c('X','Y')),]


#3. Add column to indicate DMPs =================================================================================================
if(is.null(celltypes)){
	# i.e. for bulk
	res$DMP <- res[,Pcol]<9e-8
}else{
	# i.e. for FANS
	for(i in 1:length(celltypes)){ # a separate DMP column per cell type
		ct <- celltypes[i]
		res[,paste('DMP',ct,sep='.')] <- res[,Pcols[i]]<9e-8
	}
}
res2 <- res

# load cell-type names
cells <- readRDS(paste0(refPath,"tissueNames_peakEnrichment.rds"))

for(cell in cells){

	print(cell)
	inputFile <- paste0(refPath,"top10000_",cell,"_peaks.csv")
	
	peaks <- data.frame(fread(inputFile, data.table=F),row.names=1)
	peaks <- peaks[-which(peaks$chr %in% c('X','Y')),]

	cpg_in_peak <- function(x){  # annotate whether CpG is within peak
		cpg <- as.numeric(x[which(colnames(res_chr)=='MAPINFO')])
		in_peak <- any(cpg>peaks_chr$start & cpg<peaks_chr$end)
		return(in_peak)
	}

	res_annot <- data.frame()
	for(i in 1:22){  # loop over chr 1-22
		print(i)
		peaks_chr <- peaks[peaks$chr==i,]
		res_chr <- res[res$CHR==i,]
		res_chr$Peak <- apply(res_chr,1,cpg_in_peak)
		res_annot <- rbind(res_annot,res_chr)
	}

	colnames(res_annot)[colnames(res_annot)=='Peak'] <- cell
	res2 <- res2[match(rownames(res_annot),rownames(res2)),]
	if(identical(rownames(res2),rownames(res_annot))){
		res2[,cell] <- res_annot[,cell]
	}else{
		stop("rownames(res2)!=rownames(res_annot)")
	}
}

saveRDS(res2, file=paste0(outputPath,outputFile))