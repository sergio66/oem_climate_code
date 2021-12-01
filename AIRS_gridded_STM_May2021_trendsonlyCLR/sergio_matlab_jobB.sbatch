#!/bin/bash
#
# clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly.m clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-64 sergio_matlab_jobB.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbin == TRENDS
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#                 kleenslurm; sbatch  -p high_mem                      --array=1-64 sergio_matlab_jobB.sbatch 1     == TRENDS
#
# run this with   kleenslurm; sbatch  --constraint=hpcf2013       --array=1-40 sergio_matlab_jobB.sbatch            for clust_run_retrieval_setlatbin_AIRS_loop_lonbinV2
#   or            kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 2
#
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 3     for clust_convert_strowrates2oemrates_anomaly_nonucal
# run this with   kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 4     for clust_convert_strowrates2oemrates_anomaly_nonucal_2012
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
##SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000 
#SBATCH --mem-per-cpu=6000 

#SBATCH --account=pi_strow

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
fi
