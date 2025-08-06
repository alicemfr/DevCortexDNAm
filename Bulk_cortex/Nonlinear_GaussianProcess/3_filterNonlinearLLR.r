

## Characterise and refine nonlinear probes identified with GPMethylation ##


library(data.table)
library(ggplot2)
library(gridExtra)


#1. Load data used in GPMethylation =============================================================================================
load(paste0(dataPath, "FetalBrain_Betas_LOOZ_EX3_noHiLowconst_newmethod_Feb23.RData"))
betas <- betas.dasen #543796


#2. Load GPMethylation results ==================================================================================================
gpstats <- fread(paste0(resultsPath, "gpstats.tsv"), header=T, sep='\t', data.table=F) #543796


#3. Filter to nonlinear probes ==================================================================================================
gpstats.nonlinear <- gpstats[which(gpstats$PreferedModelLR=='Mat52'),]      # 175886. nonlinear when preferred model is Matern5/2
gpstats.nonlinear <- gpstats.nonlinear[order(gpstats.nonlinear$mat52_ell),] # order on timescale (smallest = most nonlinear pattern)
rownames(gpstats.nonlinear) <- gpstats.nonlinear$Probe


#4. Calculate LLR and SNR =======================================================================================================
# Log-likelihood ratio (LLR)
matvsconst <- gpstats.nonlinear$mat52_mll - gpstats.nonlinear$const_mll # Matern 5/2 minus Constant maximum log likelihood (how much nonlinear over constant)
matvslin <- gpstats.nonlinear$mat52_mll - gpstats.nonlinear$linear_mll  # Matern 5/2 minus Linear maximum log likelihood (how much nonlinear over linear)

# Signal-to-noise ratio (SNR)
snr <- log10(gpstats.nonlinear$mat52_sigf/gpstats.nonlinear$mat52_sign)

# Mean methylation per probe
betas.match <- betas[which(rownames(betas) %in% rownames(gpstats.nonlinear)),]
betas.match <- betas.match[match(rownames(gpstats.nonlinear), rownames(betas.match)),]
identical(rownames(gpstats.nonlinear), rownames(betas.match))
meanMeth <- apply(betas.match, 1, mean, na.rm=T)

# Combine
df <- data.frame(Timescale=gpstats.nonlinear$mat52_ell, logTimescale=log10(gpstats.nonlinear$mat52_ell), MatvsConst=matvsconst, MatvsLinear=matvslin, logSNR=snr, meanMeth=meanMeth, logMeanMeth=log10(meanMeth))

nrow(df); range(df$Timescale)                               # 175886 probes. Timescale range: 3.118077 - 390.813925


#5. Filter probes on LLR ========================================================================================================
# remove any probes with LLR < 2 for either MatvsConst or MatvsLinear
df2 <- df[-which(df$MatvsConst<2 | df$MatvsLinear<2),]      # 73,531 remaining. Timescale range: 3.118077 - 282.707953

# An outlier probe with Timescale of 282 wasn't removed because its LLR is 2.744128 for both MatvsConst and MatvsLinear.
# Removing this probe manually as biologically unrealistic timescale.
df2 <- df2[-which(df2$Timescale>282),]                      # 73,530 remaining. Timescale range: 3.118077 - 93.283957

# Remove 495 probes with timescale < 10 as outlier-probe
df3 <- df2[-which(df2$Timescale<10),]                       # 73,035 remaining. Timescale range: 10.00414 - 93.28396

saveRDS(df3, file=paste0(postprocessPath,"nonlinearMetrics.rds"))


#6. Annotate nonlinear probes ===================================================================================================
d <- readRDS(paste0(postprocessPath,"nonlinearMetrics.rds"))

epicManifest <- fread(past0(refPath, "MethylationEPIC_v-1-0_B4.csv"), skip=7, fill=TRUE, data.table=F)
epicMan <- epicManifest[match(rownames(d), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
d <- cbind(d, as.data.frame(epicMan))
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }
d$Gene <- unlist(lapply(d$UCSC_RefGene_Name, uniqueAnno))
d <- d[,-which(colnames(d) %in% c('MAPINFO', 'UCSC_RefGene_Name','logMeanMeth'))]
write.csv(d, file=paste0(postprocessPath,"nonlinearProbes.csv"))


#7. Plot MLLs ===================================================================================================================

# Timescale: MatvsConstMLL
pdf(paste0(postprocessPath,"timescale_MatvsConstMLL.pdf"), width=8, height=7)
ggplot(df, aes(x=Timescale, y=MatvsConst) ) +
  geom_bin2d(bins = nrow(df)/1000) +
  scale_fill_continuous(type = "viridis") +
  xlab("Timescale") + ylab("LLR (Matern 5/2 - Constant)") +
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23)) +
  geom_hline(yintercept=c(2,5,10), linetype="dashed") +
  annotate("text", x=c(460,460,460), y=c(2,5,10), label=c("LLR = 2","LLR = 5","LLR=10"), size=5) + coord_cartesian(xlim=c(min(df$Timescale), 410), clip="off")
dev.off()


# Log timescale: MatvsConstMLL
pdf(paste0(postprocessPath,"log.timescale_MatvsConstMLL.pdf"), width=8, height=7)
ggplot(df, aes(x=logTimescale, y=MatvsConst) ) +
  geom_bin2d(bins = nrow(df)/1000) +
  scale_fill_continuous(type = "viridis") +
  xlab("log(Timescale)") + ylab("LLR (Matern 5/2 - Constant)") +
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))+
  geom_hline(yintercept=c(2,5,10), linetype="dashed") +
  annotate("text", x=c(3.3,3.3,3.3), y=c(2,5,10), label=c("LLR = 2","LLR = 5","LLR=10"), size=5) + coord_cartesian(xlim=c(min(df$logTimescale), 3), clip="off")
dev.off()


# Timescale: MatvsLinearMLL
pdf(paste0(postprocessPath,"timescale_MatvsLinearMLL.pdf"), width=8, height=7)
ggplot(df, aes(x=Timescale, y=MatvsLinear) ) +
  geom_bin2d(bins = nrow(df)/1000) +
  scale_fill_continuous(type = "viridis") +
  xlab("Timescale") + ylab("LLR (Matern 5/2 - Linear)") +
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))+
  geom_hline(yintercept=c(2,5,10), linetype="dashed") +
  annotate("text", x=c(460,460,460), y=c(2,5,10), label=c("LLR = 2","LLR = 5","LLR=10"), size=5) + coord_cartesian(xlim=c(min(df$Timescale), 410), clip="off")
dev.off()


# Log timescale: MatvsLinearMLL
pdf(paste0(postprocessPath,"log.timescale_MatvsLinearMLL.pdf"), width=8, height=7)
ggplot(df, aes(x=logTimescale, y=MatvsLinear) ) +
  geom_bin2d(bins = nrow(df)/1000) +
  scale_fill_continuous(type = "viridis") +
  xlab("log(Timescale)") + ylab("LLR (Matern 5/2 - Linear)") +
  theme_minimal() +
  theme(axis.text=element_text(size=19), axis.title=element_text(size=22)) +
  theme(plot.title = element_text(size=23))+
  geom_hline(yintercept=c(2,5,10), linetype="dashed") +
  annotate("text", x=c(3.3,3.3,3.3), y=c(2,5,10), label=c("LLR = 2","LLR = 5","LLR=10"), size=5) + coord_cartesian(xlim=c(min(df$logTimescale), 3), clip="off")
dev.off()