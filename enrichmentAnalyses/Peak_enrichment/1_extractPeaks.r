

## Enrichment of ATAC-seq peaks within EWAS results ##

# Corresponding paper: Domcke et al (2020) https://www.science.org/doi/10.1126/science.aba7612
# sci-ATAC-seq data (85,261 cells) generated on 3 human fetal individuals (2 male, 1 female; 110-115 days post-conception, equivalent to 15.7-16.4pcw)
# Top 10,000 most specific peaks for each cell type: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149683


library(data.table)
library(stringr)


#1. Load data ===================================================================================================================
spec <- fread("GSE149683_File_S4.Specificity_scores_for_top_10000_sites_per_cell_type.txt", data.table=F) #790852


#2. Reformat cell-type names ====================================================================================================

# change spaces (' ') to underscores ('_') for these would-be-duplicate cell types:
exceptions <- c('ENS glia','ENS neurons','Retinal progenitors','Retinal pigment')
for(cell in exceptions){
	spec$cell_type <- gsub(cell, gsub(' ','_',cell), spec$cell_type)
}
all_tissues <- unique(spec$cell_type)
tissues <- word(all_tissues)              # extract first word of cell type
length(tissues)==length(unique(tissues))

for(cell in tissues){
	print(cell)
	sp <- spec[grep(cell, spec$cell_type),]
	sp <- sp[sp$zscore_filtered==FALSE,]  # limit to top 10,000 using zscore_filtered column
	peak_coords <- sp$feature             # extract coords
	peak_coords_split <- matrix(unlist(strsplit(peak_coords, '_')),ncol=3,byrow=T) #10000 3
	peak_coords_split <- as.data.frame(peak_coords_split)
	colnames(peak_coords_split) <- c('chr','start','end')
	peak_coords_split$chr <- gsub("chr","", peak_coords_split$chr)
	write.csv(peak_coords_split, file=paste0("top10000_",cell,"_peaks.csv"))
}

saveRDS(tissues, file="tissueNames_peakEnrichment.rds")

tissue_names <- data.frame(`Short name`=tissues, `Full name`=all_tissues)
write.csv(tissue_names, file="tissueNamesKey.csv")