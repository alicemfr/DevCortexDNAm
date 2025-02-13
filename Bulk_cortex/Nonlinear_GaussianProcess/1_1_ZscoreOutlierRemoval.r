

## Leave-one-out (LOO) Z-score outlier removal ## 

# Calculate Z-score for each sample by removing the sample from the dataset and calculating a Z-score for the remaining samples.
# Samples that create a large Z-score when they are removed from the dataset suggest they have a large influence on the distribution.
# Outlier defined as Z-score > 5sd. Usually use cut-off of 3sd, however using a more conservative threshold of 5 to account for the 
# dynamic nature of the data and prevent the removal of true dramatic biological trends seen at the extreme ages.

setwd(PathToBetas)

#1. Load data and filter samples ================================================================================================
load(paste0(PathToBetas, "fetalBulk_EX3_23pcw_n91.rdat"))
identical(rownames(pheno),colnames(betas)) 


#2. Exclude SNP probes from all ethnicities following McCartney et al. DOI: 10.1016/j.gdata.2016.05.012 =========================
snpProbes <- read.table(paste0(refPath,"SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
snpProbes <- snpProbes[which((snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95) | (snpProbes$AFR_AF >= 0.05 & snpProbes$AFR_AF <= 0.95) | (snpProbes$AMR_AF >= 0.05 & snpProbes$AMR_AF <= 0.95) | (snpProbes$SAS_AF >= 0.05 & snpProbes$SAS_AF <= 0.95) | (snpProbes$EAS_AF >= 0.05 & snpProbes$EAS_AF <= 0.95)),]
betas <- betas[-which(rownames(betas) %in% unique(snpProbes$IlmnID)), ]


#3. For each probe, calculate Z-score for all samples in a leave-one-out (LOO) fashion ==========================================
allOutliers <- matrix(NA, nrow=1, ncol=2) # keep note of which samples are failing
colnames(allOutliers) <- c('Probe','Sample')

for(i in 1:nrow(betas)){                                    # Loop over each probe (row)

	if(i%%10000==0){print(i)}                               # Prints every 10,000th iteration, to keep track of progress
		
		# Calculate Z-scores for each sample (column)
		Zscores <- c()
		for(j in 1:ncol(betas)){									
		betas.sample <- betas[i,j]                          # DNAm of sample of interest (SOI)
			if(!is.na(betas.sample)){                       # Skip if sample has NA DNAm
				betas.without.sample <- betas[i,-j]         # DNAm of all other samples (!SOI)
				mn <- mean(betas.without.sample,na.rm=T)    # Calculate mean (mean_!SOI)
				std <- sd(betas.without.sample,na.rm=T)     # Calculate sd (sd_!SOI)
				Z <- (betas.sample - mn)/std                # Z-score = (SOI - mean_!SOI) / sd_!SOI	
			}else{
				Z <- NA                                     # If the data point is NA, assign it NA Z-score
			}									
		Zscores <- c(Zscores, Z)                            # Aggregate all Z-scores for that probe (row)
		}
		
		# Outlier identification & removal
		outliers <- which(abs(Zscores)>5)                   # Define outlier as sample > 5 sd away from mean
		for(p in 1:length(outliers)){                       # Add probe and sample(s) to running tally of outliers
			allOutliers <- rbind(allOutliers, c(rownames(betas)[i], colnames(betas)[outliers][p]))
		}
		betas[i,outliers] <- NA                             # Replace DNAm with NA for the outlier sample(s)
		
}


#4. Save results ================================================================================================================
# save outlier sample information
saveRDS(allOutliers, file="fetalBrainOutlierList_newmethod_ethsnp_Feb23.rds")
# save betas with outlier data removed
saveRDS(betas, file="FetalBrainBetasEX3_23pcw_LOOZ_newmethod_ethsnp_Feb23.rdat")