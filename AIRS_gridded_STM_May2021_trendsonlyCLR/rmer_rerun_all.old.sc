clear
date

## echo "the argument is $1"

ls -lt Output/Quantile05/test*.mat | wc -l
echo "removing latbins 01 .. 64 from Output/Quantile05/"
/bin/rm Output/Quantile05/test*.mat

/bin/rm slurm*.out
## loop forwards 1 -- 72
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 0
## loop backwards 72 -- 1
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 1
