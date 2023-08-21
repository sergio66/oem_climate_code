## this deleteds stuff from Output/QuantileX  -- one lat
## so run it as        rmer_rerun_polar.sc $1 $2

clear
date

#echo "the argument is $1"
unset c
c=$(printf %02d $1)
echo "the first argument is $c"

unset d
c=$(printf %02d $2)
echo "the second argument is $d"

echo "I assume you removed the relevant files!!!"

ls -lt Output/Quantile$c/test*.mat | wc -l

/bin/rm slurm*.out
## loop forwards 1 -- 72
sbatch  -p high_mem   --array=$d sergio_matlab_jobB.sbatch 0
## loop backwards 72 -- 1
sbatch  -p high_mem   --array=$d sergio_matlab_jobB.sbatch 1
