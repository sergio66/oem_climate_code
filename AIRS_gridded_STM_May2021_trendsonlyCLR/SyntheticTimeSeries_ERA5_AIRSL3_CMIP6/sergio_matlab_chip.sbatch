#!/bin/bash
#
# clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly.m clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-64 sergio_matlab_chip.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbin == TRENDS
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_chip.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p high_mem                      --array=1-64 sergio_matlab_chip.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p high_mem                      --array=1-64 sergio_matlab_chip.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p cpu2021                       --array=1-64 sergio_matlab_chip.sbatch 1     == TRENDS
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-40 sergio_matlab_chip.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_chip.sbatch 2
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-40 sergio_matlab_chip.sbatch 3     for clust_convert_strowrates2oemrates_anomaly_nonucal
# run this with   kleenslurm; sbatch  --exclude= --array=1-40 sergio_matlab_chip.sbatch 4     for clust_convert_strowrates2oemrates_anomaly_nonucal_2012
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 5     for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 51    for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m (can have WV const, CO2 const, T const, TH const
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 6     for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 7     for driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 8     for driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 81    for driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH.m (can have WV const, CO2 const, T const, TH const
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 9     for driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 10    for driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 11    for driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2.m
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 12    for driver_spectral_trends_latbin_1_64_sarta; make sure you set iModel correctly
#
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 13    for driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m
# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 14    for driver_spectral_trends_latbin_1_64_sarta_const_tracegas.m

# run this with   kleenslurm; sbatch  --exclude= --array=1-64 sergio_matlab_chip.sbatch 15    for driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2_constracegas.m
#
#########################################################################

# check stats after running using eg
# sacct --name=KCARTA_BREAKOUTS --state=CD -S 2016-05-07 --format=elapsed,maxrss
#

#SBATCH --job-name=ANOMALY            # Name of the job:
#SBATCH -N1                           # N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH --cpus-per-task=1             # number of cores per matlab session will be available for parallel jobs

##SBATCH --partition=2024             # desired partition
#SBATCH --partition=match             # desired partition
#SBATCH --cluster=chip-cpu            # desired cluster
#SBATCH --account=pi_strow
#SBATCH --qos=shared                  # qos to get as many cpu2024 as possible, else put pi_strow

##  This is in MB, easy peasy, 12000 works best for clust_run_retrieval_setlatbin_AIRS_loop_lonbin
##SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=12000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000
##SBATCH --mem-per-cpu=8000
#SBATCH --mem=16000

#SBATCH --time=3:59:00
#SBATCH --time=0:59:00 

########################################################################

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments, argument 1 is "
  echo $1
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

if [[ "$1" -eq "" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = DNE, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is fixed latbin, loop over 72 lonbins 
  echo "cmd line arg = 1, clust_run_retrieval_setlatbin_AIRS_loop_lonbin => JOB = latbin"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "clust_run_retrieval_setlatbin_AIRS_loop_lonbin; exit"
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

elif [[ "$1" -eq "5" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
  echo "cmd line arg = 5, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "51" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
  echo "cmd line arg = 51, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "idRH = 5; driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH; exit"

elif [[ "$1" -eq "6" ]]; then
  # makes driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
  echo "cmd line arg = 6, driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "7" ]]; then
  # makes driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
  echo "cmd line arg = 7, driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2; exit"

elif [[ "$1" -eq "8" ]]; then
  # makes driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2.m
  echo "cmd line arg = 8, driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "81" ]]; then
  # makes driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH.m
  echo "cmd line arg = 81, driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "idRH = 5; driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH; exit"

elif [[ "$1" -eq "9" ]]; then
  # makes driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2.m
  echo "cmd line arg = 9, driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AMIP6_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "10" ]]; then
  # makes driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m
  echo "cmd line arg = 10, driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2; exit"
elif [[ "$1" -eq "11" ]]; then
  # makes driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2.m
  echo "cmd line arg = 11, driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2; exit"

elif [[ "$1" -eq "12" ]]; then
  ## makes driver_spectral_trends_latbin_1_64_sarta.m, make sure you set iModel correctly
  echo "cmd line arg = 12, driver_spectral_trends_latbin_1_64_sarta"
  ## srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_spectral_trends_latbin_1_64_sarta; exit"
  %srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "iModel = +3; iNumYears = 20; driver_spectral_trends_latbin_1_64_sarta; exit"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "iModel = +51; iNumYears = 4; driver_spectral_trends_latbin_1_64_sarta; exit"

elif [[ "$1" -eq "13" ]]; then
  # makes driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m
  echo "cmd line arg = 13, driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas; exit"
elif [[ "$1" -eq "14" ]]; then
  # makes driver_spectral_trends_latbin_1_64_sarta_const_tracegas.m
  echo "cmd line arg = 14, driver_spectral_trends_latbin_1_64_sarta_const_tracegas"
  srun /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -nosplash -r "driver_spectral_trends_latbin_1_64_sarta_const_tracegas; exit"

fi
