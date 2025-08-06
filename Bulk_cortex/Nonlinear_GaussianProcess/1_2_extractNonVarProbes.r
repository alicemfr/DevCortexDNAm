

## Filter probes on IQR ##

# path.betas = full path of DNAm matrix (betas), including filename. Can take .rdat or .rds. If .rdat, function assumes DNAm matrix object contains 'betas' in name.
# path.output = directory path of output file
# name.output = filename of output file. Will be appended to path.output
# nprobes = optional: if to run on subset of the data e.g. to test functionality, specify how many probes.


extractNonVarProbes <- function(path.betas, path.output, name.output, nprobes=NULL){

	print("Step 1 - Load betas")													
	if(grepl(".rds",path.betas)){                                   # If file has .rds suffix, uses readRDS() to load
		betas <- readRDS(path.betas)
	}else{                                                          # Else, assumes .rdat
		e <- try(load(path.betas))                                  # If only 1 object was saved in .rdat format, load() does not work
		if("try-error" %in% class(e)){                              # If error when using load(), use readRDS() instead
			print("Error when loading .rdat, probably because only 1 object present. Treating as .rds and using readRDS() instead.")
			betas <- readRDS(path.betas)
		}else{
			load(path.betas)
			betas <- get(ls()[grep('betas',ls())])                  # Assuming betas object contains 'betas' in name
		}
	}

	if(is.null(nprobes)){                                           # If nprobes not specified, nprobes equals all probes in betas
		nprobes <- nrow(betas)
	}
	
	var.filter <- function(x){
		probe <- betas[i,]
		perc90 <- as.numeric(sort(probe)[0.90*length(probe)])       # 90th percentile
		perc10 <- as.numeric(sort(probe)[0.10*length(probe)])       # 10th percentile
		mid80 <- perc90 - perc10                                    # Range of middle 80% values
			if(mid80<0.05){                                         # If range is < 5%, probe considered non-variable
				nonVarProbe <- rownames(betas)[i]				
				return(c(nonVarProbe, mid80))                       # Return the probe name and its middle 80% range
			}
	}
	
	print("Step 2 - Load doParallel and make clusters")             # Prepare clusters for running in parallel
		library(doParallel)
		cl <- makeCluster(10)
		registerDoParallel(cl)

	print("Step 3")
		print(paste0("Starting var.filter() on ", nprobes))
		clusterExport(cl, list("var.filter"), envir=environment())
		res <- foreach(i=1:nprobes, .combine=rbind, .verbose=F) %dopar%{
			var.filter(betas[i,])                                   # Run var.filter on each row (probe) from 1 to nprobes
		}

	print("Step 4")
		print("Finished var.filter()")
		colnames(res) <- c("Probe","mid80Range")
		rownames(res) <- res[,'Probe']
		
	print("Step 5")
		print(paste0("Saving results to ",paste0(path.output, name.output)))
		saveRDS(res, file=paste0(path.output, name.output))
	
}