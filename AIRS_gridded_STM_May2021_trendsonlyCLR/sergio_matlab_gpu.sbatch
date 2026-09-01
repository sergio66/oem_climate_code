#!/bin/bash
#
# clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly.m clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
#   can take two arguments, first one is switch for 1,2,3,4,5,6,7 and second is +/-1 for D/N or forward/backward
#
# run this with   kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 or 1  +/- 1   == TRENDS, looping forward(0) or backward(1)  +/- 1 is for N/D
#   or            kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 or 1  +/- 1   == TRENDS, looping forward(0) or backward(1)  +/- 1 is for N/D
#                 kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 or 1  +/- 1   == TRENDS, looping forward(0) or backward(1)  +/- 1 is for N/D
#                 kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 or 1  +/- 1   == TRENDS, looping forward(0) or backward(1)  +/- 1 is for N/D
#                 kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 or 1  +/- 1   == TRENDS, looping forward(0) or backward(1)  +/- 1 is for N/D
#
#       kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch F/B(0/1) N/D(+1,-1)  for forwards/backwards  night/day
# so eg kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch       for night
# so eg kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 +1  for night
# so eg kleenslurm; sbatch   --array=1-64 sergio_matlab_gpu.sbatch 0 -1  for day

########################################################################
#
# run this with   kleenslurm; sbatch   --array=1-40 sergio_matlab_gpu.sbatch      for clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2
#   or            kleenslurm; sbatch   --array=1-40 sergio_matlab_gpu.sbatch 2 
#
# run this with   kleenslurm; sbatch   --array=1-40 sergio_matlab_gpu.sbatch 3    for clust_convert_strowrates2oemrates_anomaly_nonucal
#                  though see driver_put_together_QuantileChoose_trends.m 
#                  though see driver_put_together_QuantileChoose_trends.m 
# run this with   kleenslurm; sbatch  --array=1-40 sergio_matlab_gpu.sbatch 4    for clust_convert_strowrates2oemrates_anomaly_nonucal_2012
#
#
#########################
#
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 5     for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 51    for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 6     for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 7     for driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
#
# run this with   kleenslurm; sbatch  --array=1-64 sergio_matlab_gpu.sbatch 5    for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --array=1-64 sergio_matlab_gpu.sbatch 51   for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
# run this with   kleenslurm; sbatch  --array=1-64 sergio_matlab_gpu.sbatch 6    for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --array=1-64 sergio_matlab_gpu.sbatch 7    for driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
#
# check stats after running using eg
# sacct --name=KCARTA_BREAKOUTS --state=CD -S 2016-05-07 --format=elapsed,maxrss
#
########################################################################
#
# to put together KCARTA T/WV/O3 97 layer jacs with CO2/N2O/CH4/CFC11CFC12 col jacs
# run this with   kleenslurm; sbatch  --array=1-64 sergio_matlab_gpu.sbatch 8    for clust_put_together_jacs_clrERA5.m
#
########################################################################
########################################################################
########################################################################

#SBATCH --job-name=ANOMALYorTREND     ## name
#SBATCH -N1                           ## number of job step is to be allocated per instance of matlab
#SBATCH --cpus-per-task 1             ## tasks per node/number of cores per matlab session will be

#SBATCH --gres=gpu:1
#SBATCH --partition=gpu              ## desired partition
##SBATCH --partition=match             ## desired partition
#SBATCH --cluster=chip-gpu            ## desired cluster
#SBATCH --account=pi_strow
#SBATCH --qos=shared                  ## qos to get as many cpu2024 as possible, else put pi_strow

########################################################################

#### this is ideal
##SBATCH --qos=short+
#SBATCH --time=0:59:00 

### this is hope but does not work
##SBATCH --qos=medium+
##SBATCH --time=1:59:00 

#  Specify the desired partition develop : use high mem/medium for 40, time = 4 hr
##SBATCH --qos=medium+
##SBATCH --time=3:59:00 

##  This is in MB, less aggressive
##SBATCH --mem-per-cpu=12000
##  This is in MB, very lean, don;t go below this as java starts crying
##SBATCH --mem-per-cpu=4000
##  This is in MB, very aggressive but I have been running outta memory, CRIS convolve
##SBATCH --mem-per-cpu=24000

##  This is in MB, easy peasy, 12000 works best for clust_run_retrieval_setlatbin_AIRS_loop_lonbin
##SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=12000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000 
##SBATCH --mem-per-cpu=8000
#SBATCH --mem=16000 

##  Specify the job array (format:  start-stop:step) 
## AIRS
##srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly; exit"
##srun matlab -nojvm -nodisplay -nosplash -r "clust_run_retrieval_latbins_AIRS_loop_anomaly; exit"
##
## IASI
##srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly_11year_iasitimespan; exit"
##

if [ $# -gt 1 ]; then
  echo "Your command line contains $# arguments : "
  echo $1
  echo $2
elif [ $# -eq 1 ]; then
  echo "Your command line contains $# argument : "
  echo $1
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

## for $1 = 0,1  iNorD = -1/+1
if [[ "$1" -eq ""  && "$2" -eq "" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = DNE, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Forward, N"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = +1; iForward = +1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "0" && "$2" -eq "" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 0, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Forward, N"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = +1; iForward = +1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
  
elif [[ "$1" -eq "0" && "$2" -eq "1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 0, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Forward, N"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = +1; iForward = +1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "0" && "$2" -eq "-1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 0, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Forward, D"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = -1; iForward = +1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"

elif [[ "$1" -eq "1" && "$2" -eq "" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 1, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Backward, N"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = +1; iForward = -1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "1" && "$2" -eq "1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 1, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Backward, N"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = +1; iForward = -1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "1" && "$2" -eq "-1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 1, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin, Backward, D"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "NorD = -1; iForward = -1; clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"

elif [[ "$1" -eq "2" ]]; then
  # this fixed latbin, loop over 72 lonbins
  echo "cmd line arg = 2, clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2 => JOB = latbin"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2; exit"
elif [[ "$1" -eq "3" ]]; then
  # makes spectral anomalies
  echo "cmd line arg = 3, clust_convert_strowrates2oemrates_anomaly"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_convert_strowrates2oemrates_anomaly_nonucal; exit"
elif [[ "$1" -eq "4" ]]; then
  # makes spectral anomalies from 2102
  echo "cmd line arg = 4, clust_convert_strowrates2oemrates_anomaly"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_convert_strowrates2oemrates_anomaly_nonucal_2012; exit"

#########################
  
elif [[ "$1" -eq "5" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
  echo "cmd line arg = 5, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "cd SyntheticTimeSeries_ERA5_AIRSL3_CMIP6; driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "51" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
  echo "cmd line arg = 5, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "cd SyntheticTimeSeries_ERA5_AIRSL3_CMIP6; driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH; exit"
elif [[ "$1" -eq "6" ]]; then
  # makes driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
  echo "cmd line arg = 6, driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "cd SyntheticTimeSeries_ERA5_AIRSL3_CMIP6; driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "7" ]]; then
  # makes driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
  echo "cmd line arg = 7, driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "cd SyntheticTimeSeries_ERA5_AIRSL3_CMIP6; driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2; exit"

#########################
  
elif [[ "$1" -eq "8" ]]; then
  # makes clust_put_together_jacs_clrERA5.m
  echo "cmd line arg = 8, clust_put_together_jacs_clrERA5.m"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_put_together_jacs_clrERA5; exit"
else
  echo "No valid choice made; quitting"
fi
