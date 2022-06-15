rtime_obs  = [];
rtime_cal  = [];
btanom_obs = [];
btanom_cal = [];

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
