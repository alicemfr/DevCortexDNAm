

## Extract non variable probes from fetal brain DNAm data ##

source(paste0(scriptsPath, "1_2_extractNonVarProbes.r"))

# bulk fetal
extractNonVarProbes(path.betas=paste0(PathToBetas,"fetalBulk_EX3_23pcw_n91.rdat"), path.output=outputPath, name.output="fetalBulk_EX3_23pcw_n91.nonVarprobes.rds")


# fetal FANS
extractNonVarProbes(path.betas=paste0(PathToBetas,"fetalFANS_Neuronal.rdat"), path.output=outputPath, name.output="fetalFANS_Neuronal.nonVarprobes.rds")
extractNonVarProbes(path.betas=paste0(PathToBetas,"fetalFANS_Non.neuronal.rdat"), path.output=outputPath, name.output="fetalFANS_Non.neuronal.nonVarprobes.rds")