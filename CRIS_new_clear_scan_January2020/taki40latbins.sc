#!/bin/bash

# loop : latbins within +/- 50
# for i in {6..35..1}

# loop : all latbins
# for i in {1..40..1}

# loop : one latbin
# for i in 20

# loop i = select the latbin(s), clust_run_retrieval_latbins_loop_anomalyV2 loop over all 157 timesteps for these latbins
for i in {20..20}
do 
  echo "exporting slurm JOB ID $i"
  export SLURM_ARRAY_TASK_ID=$i
  matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_loop_anomalyV2; exit" > outputrun${i} &
  #cat taki40latbins.sc > outputrun${i} &
done
