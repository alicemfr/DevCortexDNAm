

## Annotate and resave EWAS results ##


#1. Load libraries & define functions ===========================================================================================
library(data.table)

# function to split UCSC_RefGene_Name into unique mentions of genes
uniqueAnno <- function(row){
	if(is.na(row)){
		row = ''
	}
	if(row != ""){
		return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
	}else{
		return(row)
	}
}

#2. Load manifest file ==========================================================================================================
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=T, data.table=F)


#3. Load results ================================================================================================================
study <- 'AgeCell' #change to the EWAS study of interest

# read in fetal and adult results and combine
res.fetal <- c()
res.adult <- c()
fetal.filename <- paste0(AnalysisPath, "July24/FANS_", study, "_EWAS_Fetal.rds")
adult.filename <- paste0(AnalysisPath, "July24/FANS_", study, "_EWAS_Adult.rds")

# load file if exists
if(file.exists(fetal.filename)){
	res.fetal <- data.frame(readRDS(fetal.filename)); colnames(res.fetal) <- paste0('Fetal.', colnames(res.fetal))
} else{ print(paste0("File does not exist: ", fetal.filename)) }
if(file.exists(adult.filename)){
	res.adult <- data.frame(readRDS(adult.filename)); colnames(res.adult) <- paste0('Adult.', colnames(res.adult))
} else{ print(paste0("File does not exist: ", adult.filename)) }

df <- data.frame(Probe=rownames(res.fetal))
f <- nrow(res.fetal); a <- nrow(res.adult) #nrow will be NULL if file does not exist

# combine dfs
if(!is.null(f)){ df <- cbind(df, res.fetal) }
if(!is.null(a)){ df <- cbind(df, res.adult) }

res <- df

# edit the column names in these cases so that they match the study name exactly
if(study=='AgeCell'){ colnames(res)[grep('Age.CellType',colnames(res))] <- gsub('Age.CellType','AgeCell',colnames(res)[grep('Age.CellType',colnames(res))]) }


#4. Annotate ====================================================================================================================
epicMan <- epicManifest[match(rownames(res), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")]
res <- cbind(as.data.frame(epicMan), res)
rownames(res) <- res$IlmnID

res <- res[-which(res$CHR=='Y'),]

if(grepl('Sex',study)){
	# Highlight probes found to be cross-hybridising to X/Y chr
	cross.df <- read.csv(paste0(refPath, "BLAT_crossXY_perc90.csv"), row.names=1)
	cross <- as.character(unique(cross.df$qName))
	res$Potential.Cross <- rep(0,nrow(res))
	res$Potential.Cross[which(rownames(res) %in% cross)] <- 1
}

# find unique mentions of genes
res$Gene <- unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))
# move 'Gene' column to after 'Relation_to_UCSC_CpG_Island'
res <- res[,c(1:which(colnames(res)=='Relation_to_UCSC_CpG_Island'), which(colnames(res)=='Gene'), which(colnames(res)=='Probe'):(ncol(res)-1))]


#5. Significant results =========================================================================================================
# Create columns indicating significance in fetal and adult.
# If analysing AgeCellSpecific model, split further by cell-type
if(study=='AgeCellSpecific'){
	if(!is.null(f)){
		res$Sig.Fetal.Neuronal <- res$Fetal.Neuronal.P.Age<9e-8
		res$Sig.Fetal.Non.neuronal <- res$Fetal.Non.neuronal.P.Age<9e-8
	}
	if(!is.null(a)){
		res$Sig.Adult.Neuronal <- res$Adult.Neuronal.P.Age<9e-8
		res$Sig.Adult.Non.neuronal <- res$Adult.Non.neuronal.P.Age<9e-8
	}
}else{
	if(!is.null(f)){ res$Sig.Fetal <- res$Fetal.Anova.P<9e-8 }
	if(!is.null(a)){ res$Sig.Adult <- res$Adult.Anova.P<9e-8 }
}

# age group to be used in filename
if(!is.null(f)){ grp <- "fetal" }
if(!is.null(a)){ grp <- paste(grp, "adult", sep='_') }

# save all results
write.csv(res, file=paste0(resPath, "FACS_", study, "_EWAS_", grp, "_anno.csv"))

# only keep probes that are significant in either fetal or adult
sigCols <- grep('Sig.',colnames(res))
allFALSE <- function(x){ f <- all(x[sigCols]==FALSE); return(f) } # probes not sig in any age group
toRemove <- apply(res, 1, allFALSE)
res.save <- res[!toRemove,] # remove probes not sig

# re-order columns
if(study=='AgeCellSpecific'){
	res.save <- res.save[,c(which(colnames(res.save)=='Fetal.Neuronal.ES.Age'):ncol(res.save), which(colnames(res.save)=='IlmnID'):which(colnames(res.save)=='Gene'))]
}else{
	res.save <- res.save[,c(which(colnames(res.save)=='Fetal.Anova.P'):ncol(res.save), which(colnames(res.save)=='IlmnID'):which(colnames(res.save)=='Gene'))]
}

# save significant results
write.csv(res.save, file=paste0(resPath, "FACS_", study, "_EWAS_", grp, "_anno_sig.csv"))