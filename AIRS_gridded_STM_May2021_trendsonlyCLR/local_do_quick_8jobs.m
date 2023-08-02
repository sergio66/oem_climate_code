%% Howard said it is fine to use 36 threads .... though if I launch 8 jpbs should that be 36/8 = 4.5 threads???
%% so lets use 8 threads, see do_quick_8jobs.sc

iUseStrowInteractCluster = -1;  %% on the superfast, for everyone to use, cluster
iUseStrowInteractCluster = +1;  %% for mere mortal groups at R1, see do_quick_8jobs.sc   CAN BE VERY DANGEROUS if your SCRIPT FORGETS TO QUIT!!!!!!!!!!!!!

%% we have 64 latbins, and planning on using 8 processors
if iUseStrowInteractCluster > 0
  iaLatbin_List2Process = 1:8;
  JOB = JOBM;
  iaLatbin_List2Process = (JOBM-1)*8 + iaLatbin_List2Process;
else
  iaLatbin_List2Process = JOB;
end

for iJ = 1 : length(iaLatbin_List2Process)
  JOBM = iaLatbin_List2Process(iJ);
  clust_run_retrieval_setlatbin_AIRS_loop_lonbin
end
