#!/bin/bash

# loop : latbins within +/- 50
# for i in {6..35..1}

# loop : all latbins
for i in {1..40..1}

do 
  echo "exporting $i times"
  export SLURM_ARRAY_TASK_ID=$i
  matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomalyV2; exit" > outputrun${i} &
  #cat taki40latbins.sc > outputrun${i} &
done
