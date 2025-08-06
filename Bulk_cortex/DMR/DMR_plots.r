

## Plot bulk fetal Age DMRs ##


#--------------------------------------------------------------------------------------------------------------------------------
# To fix bug with cache associated with AnnotationHub:
# From http://bioconductor.jp/packages/3.13/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html

moveFiles<-function(package){
	olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
	newdir <- tools::R_user_dir(package, which="cache")
	dir.create(path=newdir, recursive=TRUE)
	files <- list.files(olddir, full.names =TRUE)
	moveres <- vapply(files,
	FUN=function(fl){
	  filename = basename(fl)
	  newname = file.path(newdir, filename)
	  file.rename(fl, newname)
	},
	FUN.VALUE = logical(1))
	if(all(moveres)) unlink(olddir, recursive=TRUE)
}
moveFiles("AnnotationHub")
#--------------------------------------------------------------------------------------------------------------------------------


source("DMR_functions.r") # contains run.dmrff(), plotDMR(), plotDMR_wrap(). Sources miniman() within.


#1. Load results ================================================================================================================
res <- readRDS(paste0(resultsPath, "ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds")) # EWAS results
dmrs <- readRDS("ageReg_dmrGenes_GRangesAnnot.rds") # DMR results


#2. Plot DMRs ===================================================================================================================
# top DMR - PAX6
plotDMR_wrap(dmrs, chr=11, P.col='P.Age', ES.col='Beta.Age', pdfName="DMR_PAX6", manual.buffer=c(1000,1000,1000), ttl='PAX6', pad=100000)