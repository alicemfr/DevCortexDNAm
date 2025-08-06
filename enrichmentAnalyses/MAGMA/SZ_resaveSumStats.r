

## Resave SCZ GWAS summary statistics for input into MAGMA gene set analysis ##

library(data.table)

#1. Load SZ GWAS results ========================================================================================================
res <- fread(paste0(MAGMApath, "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv"), data.table=F)
dim(res)
#7585077      16


#2. Reformat file ==========================================================================================================================
# add column for total number of samples
res$nTotal <- res$NCAS + res$NCON

# subset to common alleles
   # results already subset to common alleles.

# subset to SNPs with high quality metric (IMPINFO)
res <- res[which(res$IMPINFO > 0.8),] #6935162


#3. Save full res ==========================================================================================================================
write.table(res, file=paste0(MAGMApath, "SZ_GWAS_sumStats_pval.txt"), row.names=F, sep='\t', quote=FALSE)


#4. Subset to SNP, CHR and POS columns =====================================================================================================
res.sub <- res[,c('ID','CHROM','POS')]
colnames(res.sub) <- c('SNP','CHR','POS')
res.sub[,1] <- as.character(res.sub[,1])

#5. Save sub res ===========================================================================================================================
write.table(res.sub, file=paste0(MAGMApath, "SZ_GWAS_sumStats.snploc"), row.names=F, col.names=F, sep='\t', quote=FALSE)