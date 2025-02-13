

## Create gene ontology (GO) reference files for pathway analysis ## 

library(GO.db)
library(org.Hs.eg.db)
'%ni%' <- Negate('%in%') # not in


GO.db
# GODb object:
# | GOSOURCENAME: Gene Ontology
# | GOSOURCEURL: http://current.geneontology.org/ontology/go-basic.obo
# | GOSOURCEDATE: 2023-07-27
# | Db type: GODb
# | package: AnnotationDbi
# | DBSCHEMA: GO_DB
# | GOEGSOURCEDATE: 2023-Sep11
# | GOEGSOURCENAME: Entrez Gene
# | GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
# | DBSCHEMAVERSION: 2.1


tb <- toTable(org.Hs.egGO2ALLEGS) #3383814       4
  # gene_id      go_id Evidence Ontology
# 1       1 GO:0008150       ND       BP
# 2       2 GO:0000003      IEA       BP
# 3       2 GO:0001553      IEA       BP


# remove genes with sources of evidence most likely to be unreliable or biased
unrel <- c("IEA", "NAS", "RCA")
tb <- tb[tb$Evidence %ni% unrel,]


tb <- tb[,c('gene_id','go_id')] #2545401
tb <- unique(tb) #1615090 # remove duplicate genes for same GO term
  # gene_id      go_id
# 1       1 GO:0008150
# 4       2 GO:0001867
# 5       2 GO:0001868
# where gene_id = EntrezID
# Convert EntrezID to gene symbol
annots <- select(org.Hs.eg.db, keys=tb$gene_id, columns=c("SYMBOL"), keytype="ENTREZID")
identical(annots$ENTREZID, tb$gene_id)


# Merge tb and annots
df <- data.frame(SYMBOL=annots$SYMBOL, GOTerms=tb$go_id, EntrezID=annots$ENTREZID) #1615090
  # SYMBOL    GOTerms EntrezID
# 1   A1BG GO:0008150        1
# 2    A2M GO:0001867        2
# 3    A2M GO:0001868        2

saveRDS(df, "GOTerms_annot_genes.rds")


# Collapse all GOTerms for a gene
mat <- matrix(NA, ncol=3, nrow=length(unique(df$SYMBOL)))
colnames(mat) <- c('SYMBOL','GOTerms','EntrezID')
mat <- as.data.frame(mat)

for(i in 1:length(unique(df$SYMBOL))){ #19310
if(i%%1000==0){print(i)}
	gene <- unique(df$SYMBOL)[i] # extract unique genes
	entrez <- unique(df$EntrezID)[i] # entrez id for gene
	GOterms <- df[df$SYMBOL==gene,'GOTerms']
	clps <- paste(GOterms, collapse='|')
	
	mat[i,] <- c(gene, clps, entrez)
}

write.csv(mat, file="GOTerms_annot_genes_collapsed.csv")


# Extract descriptions for each GO term
desc <- select(GO.db, keys=unique(df$GOTerms), columns=c("TERM","ONTOLOGY"))
colnames(desc) <- c('ID','Name','NameSpace')
# remove commas in Name as causes a problem when saving results table as csv
desc$Name <- gsub(',','',desc$Name)
write.csv(desc, file="unique_GOTerm_names.csv")


#-------------------------------------------------------------------------------------------------------------------------------#
sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0

# locale:
# [1] C

# time zone: Europe/London
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods
# [8] base

# other attached packages:
# [1] org.Hs.eg.db_3.18.0  GO.db_3.18.0         AnnotationDbi_1.64.1
# [4] IRanges_2.36.0       S4Vectors_0.40.2     Biobase_2.62.0
# [7] BiocGenerics_0.48.1

# loaded via a namespace (and not attached):
 # [1] crayon_1.5.2            vctrs_0.6.4             httr_1.4.7
 # [4] cli_3.6.1               rlang_1.1.2             DBI_1.2.3
 # [7] png_0.1-8               bit_4.0.5               RCurl_1.98-1.16
# [10] Biostrings_2.70.3       KEGGREST_1.42.0         bitops_1.0-8
# [13] fastmap_1.1.1           GenomeInfoDb_1.38.8     memoise_2.0.1
# [16] BiocManager_1.30.25     compiler_4.3.2          RSQLite_2.3.7
# [19] blob_1.2.4              pkgconfig_2.0.3         XVector_0.42.0
# [22] R6_2.5.1                GenomeInfoDbData_1.2.11 tools_4.3.2
# [25] bit64_4.0.5             zlibbioc_1.48.2         cachem_1.0.8