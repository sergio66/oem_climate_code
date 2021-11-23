yyStart = 2012; yyEnd = 2017;
yyStart = 2002; yyEnd = 2017;

iNewOrOld = +1;  %% steves new stats
iNewOrOld = -1;  %% steve/strow old stats

robs = [];
rclr = [];
rtime = [];
if iNewOrOld > 0
  for yy0 = yyStart : yyEnd
    yy0
    fin = ['/asl/data/stats/airs/clear/rtp_airicrad_era_rad_kl_' num2str(yy0) '_clear_desc_ocean.mat'];
    a = load(fin);
    robs = cat(1,robs,a.robs);
    rclr = cat(1,rclr,a.rclr);
    rtime = cat(1,rtime,a.rtime_mean);
  end
else
  for ii = 1 : 40
    ii
    finX = ['/home/strow/Work/Airs/Stability/Data/Desc/statlat' num2str(ii) '.mat'];
    finY = ['/home/strow/Work/Airs/Stability/Data/Desc/statlat' num2str(ii) '.mat'];
    aX = load(finX);
    aY = load(finY);
    aX.robs = aY.robs;
    robs(ii,:,:) = aX.robs;
    rclr(ii,:,:) = aX.rclr;
    rtime(ii,:,:) = aX.rtime_mean;
  end
  rtime = rtime';
  robs = permute(robs,[2 1 3]);
  rclr = permute(rclr,[2 1 3]);
end

inst = 'airs';
inst = 'al1c';
saveopt = 'savebig';
    fit_type = 'rclr';
    fit_type = 'robs';

%% see ../../oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/clust_make_spectralrates_sergio*
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
whos robs rclr rtime

load_fairs;
for ii = 3 : 40
  [yy,mm,dd,hh] = tai2utcSergio(rtime(:,ii));
  if yyStart == 2012
    ooS = find(yy == yyStart & mm == 05 & dd == 01);
    ooE = find(yy == yyEnd   & mm == 04 & dd == 30);
  elseif yyStart == 2002
    ooS = find(yy == yyStart & mm == 09 & dd == 01);
    ooE = find(yy == yyEnd   & mm == 08 & dd == 30);
  end

  oo = ooS : ooE;
  xclr = squeeze(rclr(oo,ii,:));
  xobs = squeeze(robs(oo,ii,:));
  xtime = rtime(oo,ii);

  %% now run off stats
  if iNewOrOld == +1
    fout   = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin'];
    foutsm = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin'];
  elseif iNewOrOld == -1
    fout   = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/strow_latbin'];
    foutsm = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/strow_latbin'];
  end
  latid = ii;
  count = ones(size(xtime)) * 1000;
  fit_robust_one_lat_genericradiances(fout,foutsm,latid,fit_type,xtime,xobs,xclr,count,inst,saveopt);

  i1231 = find(fairs >= 1231,1)
  i790 = find(fairs >= 790.5,1);
  i792 = find(fairs >= 791.5,1)
  bonk = load(['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii)]);
  plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
  pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
check_sergio_2012_vs_2002_vs_strow_2012

