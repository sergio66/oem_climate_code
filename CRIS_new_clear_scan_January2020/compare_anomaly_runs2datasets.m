function [A,C] = compare_anomaly_runs2datasets(file1A,file2A,file1C,file2C);

%{
(1) old and new sarta with finite diff strow jacs
f1a = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_results_oldchans.mat';
f1b = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_cal_results_oldchans.mat';
compare_anomaly_runs(f1a,f1b)
  Method 2 smoothed stemp rates = 0.015246 0.008036 0.007211
  Method 2 smoothed co2 rates = 1.692519 -0.225477 1.917996
  Method 2 smoothed n2o rates = 0.605222 -0.184678 0.789901
  Method 2 smoothed ch4 rates = 3.903583 -1.076112 4.979695

f2a = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_results_newchans.mat';
f2b = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_cal_results_newchans.mat';
compare_anomaly_runs(f2a,f2b)
  Method 2 smoothed stemp rates = 0.015419 0.009667 0.005751
  Method 2 smoothed co2 rates = 1.703375 -0.249410 1.952784
  Method 2 smoothed n2o rates = 0.591044 -0.208141 0.799185
  Method 2 smoothed ch4 rates = 3.999539 -1.203400 5.202939
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 8; figure(ii); clf; end

disp(' >>>>>>> AIRS 365 timesteps')
A = compare_anomaly_runs(file1A,file2A);

disp('>>>>>>>> CRIS 157 timesteps')
C = compare_anomaly_runs(file1C,file2C);

