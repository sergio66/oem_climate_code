#!/bin/bash

# condensed into master_sbatch.sc which does
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 1       for driver_computeERA5_16day_trends_desc_or_asc, iDorA = +1
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 2       for driver_computeERA5_16day_trends_desc_or_asc, iDorA = -1
#
# run this with   sbatch -p high_mem --array=1-20 sergio_matlab_jobB.sbatch 3    for driver_computeERA5_monthly_trends_desc_or_asc, iDorA = +1
# run this with   sbatch -p high_mem --array=1-20 sergio_matlab_jobB.sbatch 4    for driver_computeERA5_monthly_trends_desc_or_asc, iDorA = -1
#
# run this with   sbatch -p high_mem --array=1-20 sergio_matlab_jobB.sbatch 5    for driver_computeMERRA2_monthly_trends_desc_or_asc, iDorA = +1 (no separation)
#
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 6    for driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m which is 5 days/tons of memory, iDorA = +1, FAST
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 7    for driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m which is 5 days/tons of memory, iDorA = -1, FAST
#
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 8    for driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m which is 5 days/tons of memory, iDorA = +1, SLOW
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 9    for driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m which is 5 days/tons of memory, iDorA = -1, SLOW
#
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 10    for driver_compute_AIRSL3_trends_desc_or_asc.m which is 5 days/tons of memory, iDorA = +1, FAST
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 11    for driver_compute_AIRSL3_trends_desc_or_asc.m which is 5 days/tons of memory, iDorA = -1, FAST
#
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 12    for driver_computeMERRA2_monthly_OLR.m 1 .. 240
# run this with   sbatch -p high_mem --array=1 sergio_matlab_jobB.sbatch 13    for driver_computeERA5_monthly_OLR.m 1 .. 240
#

#  Name of the job:
#SBATCH --job-name=RUN_ERA_TRENDS

#  N specifies that 1 job step is to be allocated per instance of
#matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
###SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
##SBATCH --partition=strow
##SBATCH --partition=batch
##SBATCH --partition=high_mem
#SBATCH --partition=cpu2021

#  Specify the qos and run time (format:  dd-hh:mm:ss)
# trying 20 trapezoid layer jacs, less aggressive mem in MB
##SBATCH --qos=medium+
##SBATCH --time=15:59:00 
##SBATCH --mem-per-cpu=16000    
##SBATCH --cpus-per-task 1

### this is for sbatch 5 6 (ie the long laborious jobs)
##SBATCH --qos=long+
##SBATCH --time=10-00:00:00
##SBATCH --mem-per-cpu=360000    
##SBATCH --cpus-per-task 1

### this is for sbatch 12 13 (ie the short fast jobs)
#SBATCH --qos=short+
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=12000    
#SBATCH --cpus-per-task 1

########################################################################
### check the number of input arguments
if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi
########################################################################

#  Specify the job array (format:  start-stop:step)
if [[ "$1" -eq "" ]]; then
  echo "cmd line arg = DNE, generic run"
  matlab -nodisplay -r "iDorA=+1; driver_computeERA_16day_trends_desc_or_asc; exit"
elif [[ "$1" -eq "1" ]]; then
  echo "cmd line arg = 1; iDorA=+1; driver_computeERA_16day_trends_desc_or_asc"
  matlab -nodisplay -r "iCmdLineOption = $1; iDorA=+1; driver_computeERA_16day_trends_desc_or_asc; exit"
elif [[ "$1" -eq "2" ]]; then
  echo "cmd line arg = 2; iDorA=-1; driver_computeERA_16day_trends_desc_or_asc"
  matlab -nodisplay -r "iCmdLineOption = $1; iDorA=-1; driver_computeERA_16day_trends_desc_or_asc; exit"
elif [[ "$1" -eq "3" ]]; then
  echo "cmd line arg = 3; iDorA=+1; driver_computeERA5_monthly_trends_desc_or_asc"
  matlab -nodisplay -r "iCmdLineOption = $1; iDorA=+1; driver_computeERA5_monthly_trends_desc_or_asc; exit"
elif [[ "$1" -eq "4" ]]; then
  echo "cmd line arg = 4; iDorA=-1; driver_computeERA5_monthly_trends_desc_or_asc"
  matlab -nodisplay -r "iCmdLineOption = $1; iDorA=-1; driver_computeERA5_monthly_trends_desc_or_asc; exit"
elif [[ "$1" -eq "5" ]]; then
  echo "cmd line arg = 5; iDorA=+1; driver_computeMERRA2_monthly_trends_desc_or_asc"
  matlab -nodisplay -r "iCmdLineOption = $1; iDorA=1; driver_computeMERRA2_monthly_trends_desc_or_asc; exit"
elif [[ "$1" -eq "6" ]]; then
  echo "cmd line arg = 6; iDorA=1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN FAST"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = -1; iDorA=1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN; exit"
elif [[ "$1" -eq "7" ]]; then
  echo "cmd line arg = 7; iDorA=-1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN FAST"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = -1; iDorA=-1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN; exit"
elif [[ "$1" -eq "8" ]]; then
  echo "cmd line arg = 8;  iDorA=1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN SLOW"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = +1;  iDorA=1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN; exit"
elif [[ "$1" -eq "9" ]]; then
  echo "cmd line arg = 9;  iDorA=-1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN SLOW"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = +1;  iDorA=-1; driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN; exit"
elif [[ "$1" -eq "10" ]]; then
  echo "cmd line arg = 10;  iDorA=1; driver_compute_AIRSL3_trends_desc_or_asc.m FAST"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = -1;  iDorA=1; driver_compute_AIRSL3_trends_desc_or_asc; exit"
elif [[ "$1" -eq "11" ]]; then
  echo "cmd line arg = 11;  iDorA=-1; driver_compute_AIRSL3_trends_desc_or_asc.m FAST"
  matlab -nodisplay -r "iCmdLineOption = $1; iSlowORFast = -1;  iDorA=-1; driver_compute_AIRSL3_trends_desc_or_asc; exit"
##
elif [[ "$1" -eq "12" ]]; then
  echo "cmd line arg = 1, driver_computeMERRA2_monthly_OLR.m"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -r "driver_computeMERRA2_monthly_OLR; exit"
elif [[ "$1" -eq "13" ]]; then
  echo "cmd line arg = 1, driver_computeERA5_monthly_OLR.m"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -r "driver_computeERA5_monthly_OLR; exit"
fi
