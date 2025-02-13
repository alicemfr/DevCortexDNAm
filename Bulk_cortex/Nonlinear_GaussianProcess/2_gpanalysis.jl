

## Run GPMethylation on bulk fetal DNAm data ##
 # Julia v1.6.1

using GPMethylation

samplemetafile = "FetalBrain_Pheno_23pcw_EX3.csv"        # pheno file
betafile       = "FetalBrain_Betas_noHiLowConst.RData"   # betas RData file
outdir         = ""                                      # directory to save results
m              = 50                                      # number of points at which to calculate GP mean and var between min and max Age

run_gpregression(samplemetafile, betafile, outdir, m, verbose=false)