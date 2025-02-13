

## Identify non-variable probes from bulk fetal brain DNAm data ##

source(paste0(scriptsPath, "1_1_extractNonVarProbes.r"))
extractNonVarProbes(path.betas=paste0(PathToBetas,"FetalBrainBetasEX3_23pcw_LOOZ_newmethod_ethsnp.rdat"), path.output=outputPath, name.output="nonVarprobes.rds")