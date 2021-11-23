#!/bin/bash

# loop : latbins within +/- 50
# for i in {6..35..1}

# loop : all latbins
for i in { 1     2     3     4     5   275   276   277   278   279   280   281   282   354   355   356   357   358   359   360   361   362   363   364   365 }
do 
  echo "exporting $i times"
  export SLURM_ARRAY_TASK_ID=$i
  matlab -nodisplay -nosplash -r "clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly; exit" > outputrun${i} &
  #cat taki40latbins.sc > outputrun${i} &
done
