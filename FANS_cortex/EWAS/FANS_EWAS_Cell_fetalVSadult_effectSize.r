

## Compare fetal and adult cell-type effect sizes ##


#1. Load results ================================================================================================================
# probes significant to either fetal or adult in Cell-type EWAS
res <- read.csv(paste0(AnalysisPath, "FACS_Cell_EWAS_fetal_anno_sig.csv"), row.names=1) #454373


#2. Limit to significant sites in FETAL =========================================================================================
res.sig <- res[res$Fetal.Anova.P<9e-8,] #6531

# compare effect sizes between fetal and adult
ES.adult <- abs(res.sig$Adult.Estimate.Cell)*100
ES.fetal <- abs(res.sig$Fetal.Estimate.Cell)*100

mean(ES.adult) #22.61983
mean(ES.fetal) #17.81804

t.test(ES.adult, ES.fetal)$p.value # p = 8.652997e-98


#3. Limit to significant sites in ADULT =========================================================================================
res.sig <- res[res$Adult.Anova.P<9e-8,] #453675

# compare effect sizes between fetal and adult
ES.adult <- abs(res.sig$Adult.Estimate.Cell)*100
ES.fetal <- abs(res.sig$Fetal.Estimate.Cell)*100

mean(ES.adult) #14.42505
mean(ES.fetal) #1.609487

t.test(ES.adult, ES.fetal)$p.value # p=0 (i.e. p<1e-320)