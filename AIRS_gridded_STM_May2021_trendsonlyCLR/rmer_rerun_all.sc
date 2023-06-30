## this deleteds stuff from Output/QuantileXYZ where XYZ is 1-16; 
## so run it as        rmer_rerun_all.sc $1

clear
date

#echo "the argument is $1"
unset c
c=$(printf %02d $1)
echo "the argument is $c"

ls -lt Output/Quantile$c/test*.mat | wc -l
echo "removing latbins 01 .. 64 from Output/Quantile$c/"
/bin/rm Output/Quantile$c/test*.mat

/bin/rm slurm*.out
## loop forwards 1 -- 72
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 0
## loop backwards 72 -- 1
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 1
