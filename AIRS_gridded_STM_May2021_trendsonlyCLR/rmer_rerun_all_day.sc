## this deletes stuff from Output_Day/QuantileXYZ where XYZ is 1-16; 
## so run it as        rmer_rerun_all.sc $1 $2
## where first arg == 1-16 for Quantile .... second arg is 1[default]/2 for cpu2021/highmem

clear
date

echo " "
echo "see how many are done by typing this ... "
echo "      watch \"ls -lt Output_Day/Quantile03/*.mat | wc -l\" "
echo " "
echo "if not all finished, run these two commands"
echo "     sbatch  --exclude= --array=01-120 sergio_matlab_jobB.sbatch 0"
echo "     sbatch  --exclude= --array=01-120 sergio_matlab_jobB.sbatch 1"
echo " "

if [[ "$1" -eq "" ]]; then
  #echo "no arguments, set first arg to Q16"
  #c=16
  echo "no arguments, set first arg to Q03"
  c=3
  c=$(printf %02d 3)  
else
  #echo "the argument is $1"
  unset c
  c=$(printf %02d $1)
  echo "the first argument is $c"
fi

if [[ "$2" -eq "" ]]; then
  echo "no or one argument, set second arg to 1, send jobs to cpu2024"
else
  echo "the second argument is $2"
fi

ls -lt Output_Day/Quantile$c/test*.mat | wc -l
echo "removing latbins 01 .. 64 from Output_Day/Quantile$c/"
/bin/rm Output_Day/Quantile$c/test*.mat

/bin/rm slurm*.out

read -p "Press Enter to continue..."

if [[ "$2" -eq "" ]]; then
  echo "second arg = [] : submitting 64+ jobs, 64- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch   --array=01-64 sergio_matlab_chip.sbatch 0 -1
  ## loop backwards 72 -- 1
  sbatch   --array=01-64 sergio_matlab_chip.sbatch 1 -1
elif [[ "$2" -eq "1" ]]; then
  echo "second arg = 1 : submitting 64+ jobs, 64- jobs, to cpu2021"
  ## loop forwards 1 -- 72
  sbatch -p cpu2021  --array=01-64 sergio_matlab_chip.sbatch 0 -1
  ## loop backwards 72 -- 1
  sbatch  -p cpu2021  --array=01-64 sergio_matlab_chip.sbatch 1 -1
elif [[ "$2" -eq "2" ]]; then
  echo "second arg = 2 : submitting 64+ jobs, 64- jobs, to high_mem"
  ## loop forwards 1 -- 72
  sbatch  -p cpu2018--array=01-64 sergio_matlab_chip.sbatch 0 -1
  ## loop backwards 72 -- 1
  sbatch  -p cpu2018 --array=01-64 sergio_matlab_chip.sbatch 1 -1
else
  echo "rmer_rerun_all.sc $1 $2 : need second argument to be 1 (for cpu2021) or 2 (for cpu2018) "
fi
