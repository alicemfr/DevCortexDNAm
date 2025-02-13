

## 5. WGCNA - Combine pathway analyses and annotate with module name ##

modules <- c("turquoise", "blue", "brown", "yellow", "green", "red", "All")
df <- data.frame()

for(m in modules){
	res <- read.csv(paste0(outputPath,"nonlinear_",m,"_pathwayAnalysis_reducedTerms.csv"), row.names=1)
	res$Module <- rep(m,nrow(res))
	df <- rbind(df, res)
}

df <- df[,c(which(colnames(df)=='Module'), 1:(ncol(df)-1))]
write.csv(df, paste0(outputPath, "nonlinear_pathways_combined.csv"))