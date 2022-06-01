#!/bin/bash
#
# /bin/rm OutputAnomaly_OBS/*/anomtest*.mat as needed
#
# run this with   sbatch --constraint=hpcf2013       --array=1-157 sergio_matlab_jobB.sbatch 1         for clust_run_retrieval_latbins_loop_anomaly
#   or            sbatch --exclude=cnode[204,225,260,267] --array=1-157 sergio_matlab_jobB.sbatch 1            JOB=TIMESTEP,loop over 40 latbins
#
# run this with   sbatch --constraint=hpcf2013       --array=1-40 sergio_matlab_jobB.sbatch  2         for clust_run_retrieval_latbins_loop_anomalyV2
#   or            sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 2             JOB=LATBIN,loop over 157 timesteps
#
# run this with   sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 3 for clust_convert_cris_strowrates2oemrates_anomaly
# run this with   sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 4 for clust_convert_cris_strowrates2oemrates_anomaly_testsergio
# run this with   sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 5 for clust_make_sergio_anomaly_nonucal_2012.m
#
# check stats after running using eg
# sacct --name=KCARTA_BREAKOUTS --state=CD -S 2016-05-07 --format=elapsed,maxrss
#
#  Name of the job:
#SBATCH --job-name=ANOMALY

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop : use batch/short for 157, time = 1 hr
#SBATCH --partition=batch
#SBATCH --qos=short+
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
##SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000 
#SBATCH --mem-per-cpu=6000 

##  Specify the job array (format:  start-stop:step) 
# AIRS
#srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly; exit"
#srun matlab -nojvm -nodisplay -nosplash -r "clust_run_retrieval_latbins_AIRS_loop_anomaly; exit"
#
# IASI
#srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly_11year_iasitimespan; exit"
#

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

if [[ "$1" -eq "" ]]; then
  # this is fixed timestep, loop over 40 latbins 
  echo "cmd line arg = DNE, clust_run_retrieval_latbins_loop_anomaly => JOB = timestep ==> need 157 processors"
  srun matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_loop_anomaly; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is fixed timestep, loop over 40 latbins 
  echo "cmd line arg = 1, clust_run_retrieval_latbins_loop_anomaly => JOB = timestep ==> need 157 processors"
  srun matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_loop_anomaly; exit"
elif [[ "$1" -eq "2" ]]; then
  # this fixed latbin, loop over 157 timesteps
  echo "cmd line arg = 2, clust_run_retrieval_latbins_loop_anomalyV2 => JOB = latbin ==> need 40 processors"
  srun matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_loop_anomalyV2; exit"
elif [[ "$1" -eq "3" ]]; then
  # makes spectral anomalies based on Larrabee's statlatbins
  echo "cmd line arg = 3, clust_convert_strowrates2oemrates_anomaly"
  #srun matlab -nodisplay -nosplash -r "clust_convert_cris_strowrates2oemrates_anomaly; exit"
  srun matlab -nodisplay -nosplash -r "clust_convert_cris_strowrates2oemrates_anomalyDEBUG; exit"
elif [[ "$1" -eq "4" ]]; then
  # makes spectral anomalies based on sergio latbins
  echo "cmd line arg = 4, clust_convert_strowrates2oemrates_anomaly_testsergio"
  srun matlab -nodisplay -nosplash -r "clust_convert_cris_strowrates2oemrates_anomaly_testsergio; exit"
elif [[ "$1" -eq "5" ]]; then
  # makes spectral anomalies based on sergio codes, should be good since I have no factor of 1000 in the radiances
  echo "cmd line arg = 5, clust_make_sergio_anomaly_nonucal_2012.m"
  srun matlab -nodisplay -nosplash -r "clust_make_sergio_anomaly_nonucal_2012; exit"
fi
