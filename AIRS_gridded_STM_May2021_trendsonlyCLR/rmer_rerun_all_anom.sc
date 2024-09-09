## this deletes stuff from Output/QuantileXYZ where XYZ is 1-16; 
## so run it as        rmer_rerun_all_anom.sc $1 $2 $3
## where first arg == 1-16 for Quantile .... second arg is 1[default]/2 for cpu2021/highmem ... third arg is how many processors (default 120)

clear
date

echo " "
echo "see how many are done by typing this ... "
echo "      watch \" ls -lt Output/Quantile03/*.mat | wc -l \" "
echo " "
echo "if not all finished, run these two commands"
echo "     sbatch  -p cpu2021  --exclude= --array=01-120 sergio_matlab_jobB.sbatch 0"
echo "     sbatch  -p cpu2021  --exclude= --array=01-120 sergio_matlab_jobB.sbatch 1"
echo " "


if [[ "$1" -eq "" ]]; then
  #echo "no arguments, set first arg to 16"
  #c=16
  echo "no arguments, set first arg to 3"
  c=3
  c=$(printf %02d 3)  
else
  #echo "the argument is $1"
  unset c
  c=$(printf %02d $1)
  echo "the first argument is $c"
fi

if [[ "$2" -eq "" ]]; then
  echo "no or one argument, set second arg to 1, send jobs to cpu2021"
else
  echo "the second argument is $2"
fi

if [[ "$3" -eq "" ]]; then
  echo "no or one or two arguments, set third arg to 120"
  p=120
else
  #echo "the argument is $3"
  unset p
  p=$(printf %02d $3)
  echo "the third argument is $p"
fi

ls -lt Output/Quantile$c/test*.mat | wc -l
echo "removing latbins 01 .. 120 from Output/Quantile$c/"
/bin/rm Output/Quantile$c/test*.mat

/bin/rm slurm*.out

if [[ "$2" -eq "" ]]; then
  echo "second arg = [] : submitting 120+ jobs, 120- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch  -p cpu2021  --exclude=cnode018,cnode019 --array=01-120 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p cpu2021  --exclude=cnode018,cnode019 --array=01-120 sergio_matlab_jobB.sbatch 1
elif [[ "$2" -eq "1" ]]; then
  echo "second arg = 1 : submitting 120+ jobs, 120- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch  -p cpu2021 --array=01-120 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p cpu2021 --array=01-120 sergio_matlab_jobB.sbatch 1
elif [[ "$2" -eq "2" ]]; then
  echo "second arg = 2 : submitting 120+ jobs, 120- jobs, to high_mem"
  ## loop forwards 1 -- 72
  sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-120 sergio_matlab_jobB.sbatch 0
  ## loop backwards 72 -- 1
  sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-120 sergio_matlab_jobB.sbatch 1
else
  echo "rmer_rerun_all.sc $1 $2 : need second argument to be 1 (for high_mem) or 2 (for cpu2021) "
fi
