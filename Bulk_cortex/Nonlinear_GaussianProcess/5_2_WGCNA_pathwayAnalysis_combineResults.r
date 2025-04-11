

## 5. WGCNA - Combine pathway analyses and annotate with module name ##

modules <- c("turquoise", "blue", "brown", "yellow", "green", "red", "All")

df <- data.frame()
for(m in modules){
	f <- paste0(outputPath,"missMethylPathwayAnalysis_nonlinear_fetalBrain_",m,".csv")
	if(file.exists(f)){
		res <- read.csv(f, row.names=1)
		res$Module <- rep(m,nrow(res))
		df <- rbind(df, res)
	}
}

df <- df[,c(which(colnames(df)=='Module'), 1:(ncol(df)-1))]

write.csv(df, paste0(outputPath, "nonlinear_pathways_combined.csv"))