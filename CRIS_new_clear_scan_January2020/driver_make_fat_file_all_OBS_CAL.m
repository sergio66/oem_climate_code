addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

rtime_obs  = [];
rtime_cal  = [];
btanom_obs = [];
btanom_cal = [];

if ~exist('fat_file_all_OBS_CAL.mat')
  for ii = 1 : 39
    
    f1 = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(ii) '.mat'];
    f2 = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(ii) '_cal.mat'];
    a1 = load(f1);
    a2 = load(f2);
    rtime_obs(ii,:) = a1.avg16_rtime;
    btanom_obs(ii,:,:) = a1.avg16_btanom;
    rtime_cal(ii,:) = a2.avg16_rtime;
    btanom_cal(ii,:,:) = a2.avg16_btanom;
  end
  
  saver = ['save fat_file_all_OBS_CAL.mat rtime_cal rtime_obs btanom_cal btanom_obs'];
  %eval(saver)
else
  load fat_file_all_OBS_CAL.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now do the trends

load f1305.mat

%% remember these are 16 day intervals
[yy,mm,dd,hh] = tai2utcSergio(rtime_cal(20,:));
daysSince2012 = change2days(yy,mm,dd,2012);

%% look at eg /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
warning off
for jj = 1 : 39
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  [yy,mm,dd,hh] = tai2utcSergio(rtime_cal(jj,:));
  daysSince2012 = change2days(yy,mm,dd,2012);
  tobs = squeeze(btanom_obs(jj,:,:));
  tobs = tobs';
  for ii = 1 : 1305
    [a,b] = Math_tsfit_lin_robust(daysSince2012,tobs(ii,:),4);
    b_obs(jj,ii,2) = a(2);
    b_err_obs(jj,ii,2) = b.se(2);
  end

  tcal = squeeze(btanom_cal(jj,:,:));
  tcal = tcal';
  for ii = 1 : 1305
    [a,b] = Math_tsfit_lin_robust(daysSince2012,tcal(ii,:),4);
    b_cal(jj,ii,2) = a(2);
    b_err_cal(jj,ii,2) = b.se(2);
  end

end
warning on

fprintf(1,'\n')

plot(f1305,nanmean(squeeze(b_obs(:,:,2))),f1305,nanmean(squeeze(b_cal(:,:,2))))

comment = 'see driver_make_fat_file_all_OBS_CAL.m';
save cris_npp_2012_o5_to_2019_04_trends_fron_anomaly.mat comment f1305 b_obs b_err_obs b_cal b_err_cal
