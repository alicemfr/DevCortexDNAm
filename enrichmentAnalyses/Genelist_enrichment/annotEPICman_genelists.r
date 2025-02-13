

## Annotate EPIC manifest with gene lists ##

library(data.table)
library(tidyr)

uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }


#1. Load gene list files ========================================================================================================

# SFARI autism genes ------------------------------------------------------------------------------------------------------------
sfari <- read.csv("SFARI-Gene_genes_08-19-2024release_09-23-2024export.csv") #1191
sfari <- sfari[which(sfari$gene.score==1),] # limit to SFARI genes with scores S or 1. n=233
#table(sfari$syndromic) # 0=118; 1=115
sfari <- sfari[,c('gene.symbol','chromosome','gene.score','syndromic')]
sfari$syndromic2[sfari$syndromic==1] <- 'S'
sfari <- sfari %>% unite(., col = "Score", gene.score, syndromic2, na.rm=TRUE, sep = '')
sfari$syndromic <- NULL
colnames(sfari) <- c('Gene','CHR','Score')
sfari$altGeneName <- rep('',nrow(sfari)) # identified alternative gene names by searching on genecards.org
sfari$altGeneName[which(sfari$Gene=='KMT5B')] <- 'SUV420H1'
sfari$altGeneName[which(sfari$Gene=='SRPRA')] <- 'SRPR'
sfari$altGeneName[which(sfari$Gene=='NEXMIF')] <- 'KIAA2022'
write.csv(sfari, file="SFARI-Gene_genes_08-19-2024release_09-23-2024export_Score1.csv")

# Schizophrenia SCHEMA genes ----------------------------------------------------------------------------------------------------
scz <- read.csv("Singh2022_SCZ_GWAS.csv") # Supplementary Table 5 of Singh et al. (2022) https://doi.org/10.1038/s41586-022-04556-w
schema <- scz[which(scz$Q.meta<0.05),]
schema <- schema[,c('Gene.Symbol','Q.meta')]
colnames(schema) <- c('Gene','Q-value')
write.csv(schema, file="SCHEMA_genes.csv")

# Combine -----------------------------------------------------------------------------------------------------------------------
# All genes across SFARI and SCHEMA. Duplicate genes present
df <- data.frame(Gene=c(sfari$Gene, schema$Gene), GeneList=c(rep('SFARI',nrow(sfari)),rep('SCHEMA',nrow(schema))))
write.csv(df, file="GeneList_SFARI_SCHEMA.csv")


#2. Load EPIC manifest ==========================================================================================================
epicManifest <- fread(paste0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
epic <- epicManifest[,c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")] # limit to relevant columns
epic$Gene <- unlist(lapply(epic$UCSC_RefGene_Name, uniqueAnno)) # extract unique transcripts of genes


#3. Annotate EPIC manifest with gene lists ======================================================================================
annotEPIC_glist <- function(epic, genelist, glist_name){
  epic[,glist_name] <- rep(0,nrow(epic))
  epic[,paste(glist_name,'Gene',sep='.')] <- rep('',nrow(epic))
  
  prg <- txtProgressBar()
  for(i in 1:nrow(genelist)){
    setTxtProgressBar(prg, i/nrow(genelist))
    
    # first try original gene name
    g <- as.character(genelist$Gene[i])
    print(g)
    indx <- grep( paste( paste0("^",g,"$"), paste0("^",g,";"), paste0(";",g,"$"), paste0(";",g,";"), sep='|'), epic$Gene) #regex = ^gene$|^gene;|;gene$|;gene;
    
    if(length(indx)>0){
      epic[indx,paste(glist_name,'Gene',sep='.')] <- g
    }
    
    # if not found, try alternative gene name
    if(length(indx)==0){ 
      print(paste(g, "had no matches."))
      a <- as.character(genelist$altGeneName[i])
      if(g!=a){
        print(paste("Trying alternative gene name:", a))
        indx <- grep( paste( paste0("^",a,"$"), paste0("^",a,";"), paste0(";",a,"$"), paste0(";",a,";"), sep='|'), epic$Gene)
        if(length(indx)!=0){
          epic[indx,paste(glist_name,'Gene',sep='.')] <- a			
        }else{
          print(paste(a, "had no matches."))
        }
      }
    }
    
    # set genelist column to 1 for probes annotated to genelist genes
    epic[indx,glist_name] <- 1
    
  }
  return(epic)
}

epic.annot <- epic

# Annotate EPIC manifest with SFARI genes...
epic.annot <- annotEPIC_glist(epic.annot, genelist=sfari, glist_name='SFARI')

# ...and then with SCHEMA genes
epic.annot <- annotEPIC_glist(epic.annot, genelist=schema, glist_name='SCHEMA')

write.csv(epic.annot, file="EPIC_annot_SFARI_SCHEMA.csv")