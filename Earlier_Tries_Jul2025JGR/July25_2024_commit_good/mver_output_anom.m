function mver_output_anom(str)

if exist('anomaly_0dayavg_results.mat')
  mver = ['!/bin/mv -i anomaly_0dayavg_results.mat anomaly_0dayavg_results_' str '.mat'];
  eval(mver)
else
  disp('anomaly_0dayavg_results.mat DNE');
end

if exist('anomaly_0dayavg_results_spectra.mat')
  mver = ['!/bin/mv -i anomaly_0dayavg_results_spectra.mat  anomaly_0dayavg_results_spectra_' str '.mat'];
  eval(mver)
else
  disp('anomaly_0dayavg_results_spectra.mat DNE')
end

if exist('anomaly_0dayavg_cal_results.mat')
  mver = ['!/bin/mv -i anomaly_0dayavg_cal_results.mat anomaly_0dayavg_cal_results_' str '.mat'];
  eval(mver)
else
  disp('anomaly_0dayavg_cal_results.mat DNE');
end

if exist('anomaly_0dayavg_results_spectra_cal.mat')
  mver = ['!/bin/mv -i anomaly_0dayavg_results_spectra_cal.mat  anomaly_0dayavg_results_spectra_cal_' str '.mat'];
  eval(mver)
else
  disp('anomaly_0dayavg_results_spectra_cal.mat DNE')
end
