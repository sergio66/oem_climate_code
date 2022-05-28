#!/bin/bash

# loop : latbins within +/- 50
# for i in {6..35..1}

# loop i = select the timestep(s), clust_run_retrieval_latbins_loop_anomaly loop over all 40 latbins for these timesteps, note no paranetheses
#for i in  1     2     3     4     5    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27   354   355   356   357   358   359   360   361   362   363   364   365
# loop i = select the timestep(s), clust_run_retrieval_latbins_loop_anomaly loop over all 40 latbins for these timesteps    2016.00 = okdates(85) = worst retrieval
for i in {85..85}
do 
  echo "exporting slurm JOB ID $i"
  export SLURM_ARRAY_TASK_ID=$i
  matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_loop_anomaly; exit" > outputrun${i} &
  #cat taki40latbins.sc > outputrun${i} &
done
