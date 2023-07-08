## this deletes stuff from Output/QuantileXYZ where XYZ is 1-16; 
## so run it as        rmer_rerun_all.sc $1 $2
## where first arg == 1-16 for Quantile .... second arg is 1[default]/2 for cpu2021/highmem

clear
date

if [[ "$1" -eq "" ]]; then
  echo "no arguments, set $1 = 16"
  c=16
else
  #echo "the argument is $1"
  unset c
  c=$(printf %02d $1)
  echo "the first argument is $c"
fi

echo "the second argument is $2"

ls -lt Output_CAL/Quantile$c/test*.mat | wc -l
echo "removing latbins 01 .. 64 from Output/Quantile$c/"
/bin/rm Output_CAL/Quantile$c/test*.mat

/bin/rm slurm*.out

if [[ "$2" -eq "" ]]; then
  echo "second arg = [] : submitting 64+ jobs, 64- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch  -p cpu2021  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p cpu2021  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 1
elif [[ "$2" -eq "1" ]]; then
  echo "second arg = 1 : submitting 64+ jobs, 64- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch  -p cpu2021 --array=01-64 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p cpu2021 --array=01-64 sergio_matlab_jobB.sbatch 1
elif [[ "$2" -eq "2" ]]; then
  echo "second arg = 2 : submitting 64+ jobs, 64- jobs, to high_mem"
  ## loop forwards 1 -- 72
  sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-64 sergio_matlab_jobB.sbatch 1
else
  echo "rmer_rerun_all.sc $1 $2 : need second argument to be 1 (for high_mem) or 2 (for cpu2021) "
fi
