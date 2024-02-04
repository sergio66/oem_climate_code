#!/bin/bash
#
# clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly.m clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-64 sergio_matlab_jobB.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbin == TRENDS
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p high_mem                      --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p high_mem                      --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p cpu2021                       --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-40 sergio_matlab_jobB.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 2
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-40 sergio_matlab_jobB.sbatch 3     for clust_convert_strowrates2oemrates_anomaly_nonucal
# run this with   kleenslurm; sbatch  --exclude= --array=1-40 sergio_matlab_jobB.sbatch 4     for clust_convert_strowrates2oemrates_anomaly_nonucal_2012
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 5     for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 6     for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 7     for driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 8     for driver_check_WV_T_RH_MERRA_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 9     for driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 10    for driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 11    for driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2.m
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 12    for driver_spectral_trends_latbin_1_64_sarta; make sure you set iModel correctly
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 13    for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 14    for driver_spectral_trends_latbin_1_64_sarta_const_tracegas.m

# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_jobB.sbatch 15    for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2_constracegas.m
#
#########################################################################

# check stats after running using eg
# sacct --name=KCARTA_BREAKOUTS --state=CD -S 2016-05-07 --format=elapsed,maxrss
#
#  Name of the job:
#SBATCH --job-name=ANOMALY

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task=1

#  Specify the desired partition develop : use batch/short for 64, time = 1 hr
##SBATCH --partition=batch
##SBATCH --partition=high_mem
#SBATCH --qos=short+
#SBATCH --partition=cpu2021
#SBATCH --time=0:59:00 

#  Specify the desired partition develop : use high mem/medium for 40, time = 4 hr
##SBATCH --partition=high_mem
##SBATCH --qos=medium+
##SBATCH --time=3:59:00 

##  This is in MB, less aggressive
##SBATCH --mem-per-cpu=12000
##  This is in MB, very lean, don;t go below this as java starts crying
##SBATCH --mem-per-cpu=4000
##  This is in MB, very aggressive but I have been running outta memory, CRIS convolve
##SBATCH --mem-per-cpu=24000

##  This is in MB, easy peasy, 18000 works
#SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000 

#SBATCH --account=pi_strow

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments, argument 1 is "
  echo $1
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

if [[ "$1" -eq "" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = DNE, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 1, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "2" ]]; then
  # this fixed latbin, loop over 72 lonbins
  echo "cmd line arg = 2, clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2 => JOB = latbin"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2; exit"
elif [[ "$1" -eq "3" ]]; then
  # makes spectral anomalies
  echo "cmd line arg = 3, clust_convert_strowrates2oemrates_anomaly"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "clust_convert_strowrates2oemrates_anomaly_nonucal; exit"
elif [[ "$1" -eq "4" ]]; then
  # makes spectral anomalies from 2102
  echo "cmd line arg = 4, clust_convert_strowrates2oemrates_anomaly"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "clust_convert_strowrates2oemrates_anomaly_nonucal_2012; exit"

elif [[ "$1" -eq "5" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
  echo "cmd line arg = 5, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "6" ]]; then
  # makes driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
  echo "cmd line arg = 6, driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "7" ]]; then
  # makes driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
  echo "cmd line arg = 7, driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "8" ]]; then
  # makes driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2.m
  echo "cmd line arg = 8, driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "9" ]]; then
  # makes driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2.m
  echo "cmd line arg = 9, driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "10" ]]; then
  # makes driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m
  echo "cmd line arg = 10, driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "11" ]]; then
  # makes driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2.m
  echo "cmd line arg = 11, driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2; exit"

elif [[ "$1" -eq "12" ]]; then
  ## makes driver_spectral_trends_latbin_1_64_sarta.m, make sure you set iModel correctly
  echo "cmd line arg = 12, driver_spectral_trends_latbin_1_64_sarta"
  ## srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_spectral_trends_latbin_1_64_sarta; exit"
  %srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "iModel = +3; iNumYears = 20; driver_spectral_trends_latbin_1_64_sarta; exit"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "iModel = +51; iNumYears = 4; driver_spectral_trends_latbin_1_64_sarta; exit"

elif [[ "$1" -eq "13" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m
  echo "cmd line arg = 13, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas; exit"
elif [[ "$1" -eq "14" ]]; then
  # makes driver_spectral_trends_latbin_1_64_sarta_const_tracegas.m
  echo "cmd line arg = 14, driver_spectral_trends_latbin_1_64_sarta_const_tracegas"
  srun /usr/ebuild/software/MATLAB/2021a/bin/matlab -nodisplay -nosplash -r "driver_spectral_trends_latbin_1_64_sarta_const_tracegas; exit"

fi
