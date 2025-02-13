

## Resave ASD GWAS summary statistics for input into MAGMA gene set analysis ##

library(data.table)

#1. Load ASD GWAS results =======================================================================================================
res <- fread(paste0(MAGMApath, "iPSYCH-PGC_ASD_Nov2017"), data.table=F)
dim(res)
#9112386       9


#2. Reformat file ===============================================================================================================
# add column for total number of samples
res$nTotal <- 18381+27969 #from paper: 18,381 ASD cases and 27,969 controls

# subset to common alleles
# Paper says "... we examined markers with minor allele frequency (MAF) ≥ 0.01, imputation INFO score ≥ 0.7, and supported by an effective sample size in > 70% of the total. This final meta-analysis included results for 9,112,387 autosomal markers ...". 
# Since res has 9112386 rows the results must contain only common alleles.

# subset to SNPs with high quality metric (INFO)
res <- res[which(res$INFO > 0.8),] #8479438


#3. Save full res ===============================================================================================================
write.table(res, file=paste0(MAGMApath, "ASD_GWAS_sumStats_pval.txt"), row.names=F, sep='\t', quote=FALSE)


#4. Subset to SNP, CHR and POS columns ==========================================================================================
res.sub <- res[,c('SNP','CHR','BP')]
colnames(res.sub) <- c('SNP','CHR','POS')


#5. Save sub res ================================================================================================================
write.table(res.sub, file=paste0(MAGMApath, "ASD_GWAS_sumStats.snploc"), row.names=F, col.names=F, sep='\t', quote=FALSE)