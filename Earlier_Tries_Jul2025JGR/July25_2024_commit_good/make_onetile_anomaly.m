%% see driver_put_together_QuantileChoose_anomalies.m
%% see make_globalavg_and_onetile_anomaly.m

addpath /asl/matlib/science/
addpath /home/sergio/MATLABCODE/TIME

if ~exist('yy')
  junk = load('/asl/s1/sergio/JUNK/anomaly_ALL_Q03.mat','yy','mm','dd','hh','rtime');
  yy = junk.yy;
  mm = junk.mm;
  dd = junk.dd;
  hh = junk.hh;
  rtime = junk.rtime;
  daysSince2002 = change2days(yy,mm,dd,2002);
  load h2645structure.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear btavg* usethese

dirCode = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/';
dirData = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';

junk = input('Choose tile [LonBin LatBin] : ');
LonBin = junk(1);
LatBin = junk(2);

ii = LonBin;
jj = LatBin;

iQuant = input('Enter Quantile (3=default) : ');
if length(iQuant) == 0
  iQuant = 3;
end

fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_3_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
fname = [dirData fname '_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'];
a = load(fname,'bt_anom_desc','bt_desc');
for iCnt = 1 : 1
  clear junk
  iTile = (LatBin-1)*72 + LonBin;
  junk = squeeze(a.bt_anom_desc(iQuant,:,:));
  str = ['btavg' num2str(iCnt) ' = junk;'];
  eval(str);
end

btavgAnomFinal = btavg1;

disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
for iCnt = 1 : 1
  str = ['junk = btavg' num2str(iCnt) ';'];
  eval(str);
  figure(1); clf; pcolor(2002+daysSince2002/365,h.vchan,junk); colorbar; caxis([-1 +1]*2); shading flat; ylim([640 1640]); title(num2str(iCnt)); colormap(usa2);
  figure(2); clf; plot(h.vchan,nanmean(junk,2)); xlim([640 1640]); ylim([-1 +1]); title(num2str(iCnt)); plotaxis2;
  for ii = 1 : 2645
    P = polyfit(daysSince2002,junk(ii,:),1);
    quicktrend(ii) = P(1)*365;
  end
  figure(3); clf; plot(h.vchan,quicktrend); xlim([640 1640]); ylim([-1 +1]/10); ylabel('BT trend'); title(num2str(iCnt)); plotaxis2;
  trendavg(iCnt,:) = quicktrend;
  anomavg(iCnt,:) = nanmean(junk,2);

  figure(4); clf; plot(onidd,oni,'k',2002+daysSince2002/365,10*smooth(junk(1520,:),20),'linewidth',2); plotaxis2; xlim([2002 2023]); ylim([-2 +2])
  pause(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSave = input('save -1/+1  : ');
if iSave > 0
  comment = 'see make_global_avg_and_onetile_averages.m';
  fout = ['anomaly_tile_' num2str(iTile) '_timeseries_Q' num2str(iQuant,'%02d') '.mat'];
  saver = ['save -v7.3 ' fout ' btavgAnomFinal yy mm dd hh rtime comment anomavg iTile LonBin LatBin'];
  eval(saver);
end
