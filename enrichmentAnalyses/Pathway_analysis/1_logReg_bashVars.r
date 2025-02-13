

## Arguments to be passed to the Slurm submission script for pathway analysis ##


#===============================================================================================================================#
# bulk fetal

INPUT=ageReg_fetalBrain_EX3_23pcw_annotAllCols_filtered.rds
RESCOL=P.Age
TEST=ageReg_fetalBrain
TOP10=TRUE
SPLIT=TRUE
ESCOL=Beta.Age
DIRN=Hypo

JOB=$TEST

if [ $TOP10 = TRUE ] ; then
    JOB=${JOB}_top10perc
fi

if [ $SPLIT = TRUE ] ; then
    JOB=${JOB}_${DIRN}
fi

cd ${SCRIPTPATH}
LOG=${JOB}.log
ERR=${JOB}.err

echo $JOB

sbatch --job-name=$JOB --output=$LOG --error=$ERR 2_SubmissionScript_GO.pbs $INPUT $RESCOL $OUTPUT $TEST $TOP10 $SPLIT $ESCOL $DIRN


#===============================================================================================================================#
# FANS fetal

CELL='Non.neuronal'

INPUT=FACS_AgeCellSpecific_EWAS_fetal_adult_anno.csv
RESCOL=Fetal.${CELL}.P.Age
TEST=AgeCellSpecific_${CELL}
TOP10=FALSE
SPLIT=TRUE
ESCOL=Fetal.${CELL}.ES.Age
DIRN=Hyper

JOB=$TEST

if [ $TOP10 = TRUE ] ; then
    JOB=${JOB}_top10perc
fi

if [ $SPLIT = TRUE ] ; then
    JOB=${JOB}_${DIRN}
fi

cd ${SCRIPTPATH}
LOG=${JOB}.log
ERR=${JOB}.err

echo $JOB

sbatch --job-name=$JOB --output=$LOG --error=$ERR 2_SubmissionScript_GO.pbs $INPUT $RESCOL $OUTPUT $TEST $TOP10 $SPLIT $ESCOL $DIRN