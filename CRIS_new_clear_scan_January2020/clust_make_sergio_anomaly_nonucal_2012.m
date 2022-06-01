%% see ../../oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/clust_make_spectralrates_sergio*
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
JOB = 20

yyStart = 2012; yyEnd = 2019;

iNewOrOld = -1;  %% strow orig stats

robs = [];
rclr = [];
rtime = [];
if iNewOrOld == -1
  for ii = JOB
    ii
    finX = ['/home/strow/Work/Cris/Stability/Data/Desc/statlat' num2str(ii) '.mat'];
    aX = load(finX);
    robs(ii,:,:)  = squeeze(aX.robs(:,5,:));
    rclr(ii,:,:)  = squeeze(aX.rclr(:,5,:));
    rtime(ii,:,:) = aX.rtime_mean(:,5);
  end
  rtime = rtime';
  robs = permute(robs,[2 1 3]);
  rclr = permute(rclr,[2 1 3]);
end

inst = 'airs';
inst = 'al1c';
inst = 'cris';  %% this is automatically CIRS NSR
saveopt = 'savebig';
    fit_type = 'robs';
    fit_type = 'rclr';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos robs rclr rtime

%% load_fairs;  %% load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/airs_f.mat'); basically give f2378 and f2645
usethese0 = 1 : 1317;

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/cris_ichan_vchan_nsr_fsr.mat
usethese = find(nsr.ichan <= 1305);
fairs = nsr.vchan(usethese);

for ii = JOB
  [yy,mm,dd,hh] = tai2utcSergio(rtime(:,ii));
  if yyStart == 2012
    ooS = find(yy == yyStart & mm == 05 & dd == 01);

    ooE = find(yy == yyEnd   & mm == 04 & dd == 30);
    ooE = find(yy == yyEnd   & mm == 03 & dd == 30);
    ooE = find(yy == yyEnd   & mm == 03); ooE = ooE(end);
  end

  oo = ooS : ooE;
  xclr = squeeze(rclr(oo,ii,usethese0));
  xobs = squeeze(robs(oo,ii,usethese0));
  xtime = rtime(oo,ii);

  %% now run off stats
  if iNewOrOld == -1 & fit_type == 'robs'
    fout   = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_' num2str(yyStart) '/sergio_robs_latbin'];
    foutsm = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_' num2str(yyStart) '/sergio_robs_latbin'];
  elseif iNewOrOld == -1 & fit_type == 'rclr'
    fout   = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_' num2str(yyStart) '/sergio_rclr_latbin'];
    foutsm = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_' num2str(yyStart) '/sergio_rclr_latbin'];
  end
  latid = ii;
  count = ones(size(xtime)) * 1000;
  fit_robust_one_lat_genericradiances(fout,foutsm,latid,fit_type,xtime,xobs,xclr,count,inst,saveopt);

  iC1231 = find(fairs >= 1231,1);
  iC791 = find(fairs >= 791,1);
  iC792 = find(fairs >= 792,1);
  iC791 = find(fairs >= 791-0.5,1);
  iC792 = find(fairs >= 792-0.5,1);
  bonk = load(['ANOM_16dayavg_2012_2019/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii)]);
  plot(smooth(1000*(bonk.all_bt_anom(:,iC791)-bonk.all_bt_anom(:,iC792)),180))
  pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%check_sergio_2012_vs_2002_vs_strow_2012

