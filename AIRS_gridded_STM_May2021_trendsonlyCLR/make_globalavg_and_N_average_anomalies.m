addpath /asl/matlib/science/
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/TIME

%if ~exist('btanomD')
%  load /asl/s1/sergio/JUNK/anomaly_ALL_Q03.mat
%  daysSince2002 = change2days(yy,mm,dd,2002);
%  load h2645structure.mat
%end
if iDorA > 0
  if ~exist('daysSince2002A')
    daysSince2002A = change2days(yyD,mmD,ddD,2002);
  end
  btanomX = btanomD;
  fluxanomaly = fluxanomD;
  disp('make_globalavg_and_N_average_anomalies.m using btanomD')
else
  if ~exist('daysSince2002A')
    daysSince2002A = change2days(yyA,mmA,ddA,2002);
  end
  btanomX = btanomA;
  fluxanomaly = fluxanomA;
  disp('make_globalavg_and_N_average_anomalies.m using btanomA')
end

if ~exist('iNumAnomTimeSteps')
  [~,ind,iNumAnomTimeSteps] = size(btanomX);
  ind = 1 : ind;
  [h,ha,p,pa] = rtpread('/asl/rtp/airs/airs_l1c_v674/clear/2023/ecmwf_airicrad_day022_clear.rtp');
  iNumYears = iNumAnomTimeSteps/23;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlat65')
  do_XX_YY_from_X_Y
  [salti, landfrac] = usgs_deg10_dem(Y,X);  
  lf = landfrac; lf = lf(:);
  [mmjunk,nnjunk] = size(XX);
  if mmjunk == 1
    XX = XX';
    YY = YY';
  end
end

coslat = cos(YY*pi/180)*ones(1,iNumAnomTimeSteps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alone = load('anomaly_chID_1520_Q03.mat');
clear btavg* usethese anomavg trendavg

iCnt = 1;
usethese{iCnt} = 1 : 4608;
for ii = 1 : length(ind)
  moo = squeeze(btanomX(:,ii,:));
  btavg_cos(ii,:)  = nansum(moo.*coslat)./nansum(coslat);
  %xbtavg_cos(ii,:) = nanmean(moo.*coslat);
end
pcolor(2002+daysSince2002A/365,h.vchan(ind),btavg_cos); colorbar; shading flat;
%pcolor(2002+daysSince2002A/365,h.vchan(ind),btavg_cos-xbtavg_cos); colorbar; shading flat;
if iNumAnomTimeSteps < 460
  plot(2002+daysSince2002A/365,xbtavg_cos(1520,:),'bx-',2002+daysSince2002A/365,nanmean(alone.btanomX.*coslat,1),'linewidth',2); plotaxis2; xlim([2002 2023])
  plot(2002+daysSince2002A/365,btavg_cos(1520,:),'b',2002+daysSince2002A/365,nanmean(alone.btanomX.*coslat,1),'linewidth',2); plotaxis2; xlim([2002 2023])
else
  plot(2002+daysSince2002A/365,btavg_cos(1520,:),'b'); plotaxis2; title('BT 1231 cos wgt anomaly'); xlim([min(2002+daysSince2002A/365) max(2002+daysSince2002A/365)]); 
end

wow = btavg_cos(1520,:);  PW = polyfit(daysSince2002A,wow,1); YPW = polyval(PW,daysSince2002A);
plot(onidd,oni,'k',2002+daysSince2002A/365,10*smooth(btavg_cos(1520,:)-YPW,16*2),'r','linewidth',2); xlim([2002 2023])
plotaxis2; hl = legend('ONI','BT1231 anomaly vers2 CORRECT COSWGT','location','best');
%xow = xbtavg_cos(1520,:); PX = polyfit(daysSince2002A,xow,1); YPX = polyval(PX,daysSince2002A);
%plot(onidd,oni,'k',2002+daysSince2002A/365,10*smooth(xbtavg_cos(1520,:)-YPX,16*2),'b',2002+daysSince2002A/365,10*smooth(btavg_cos(1520,:)-YPW,16*2),'r','linewidth',2); xlim([2002 2023])
%plotaxis2; hl = legend('ONI','BT1231 anomaly vers1','BT1231 anomaly vers2 CORRECT COSWGT','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newLatGrid = [-90 -70 -50 -30 -10 +10 +30 +50 +70 +90];
newLatGrid = [-90:+10:+90];
newLatGrid = [-90 -75 -60 [-55:5:+55] +60 +75 +90];

for iCnt = 1 : length(newLatGrid)-1
  clear junk
  fprintf(1,'making avg %2i of %2i \n',iCnt, length(newLatGrid)-1);
  usethese{iCnt+1} = find(YY >= newLatGrid(iCnt) & YY < newLatGrid(iCnt+1));
  xusethese = usethese{iCnt+1};
  junkcoslat = coslat(xusethese,:);
  for ii = 1 : length(ind)
    moo = squeeze(btanomX(:,ii,:));
    moo = moo(xusethese,:);
    junk(ii,:) = nansum(moo.*junkcoslat)./nansum(junkcoslat);
  end
  str = ['btavg_newlat' num2str(iCnt) ' = junk;'];
  eval(str);
end

%% for iCnt = 1 : length(newLatGrid)-1
%%   str = ['btavg_newlat' num2str(iCnt) ' = btavg_newlat' num2str(iCnt+1) ';'];
%%   eval(str);
%% end

%btavgAnomFinal = [btavg0 btavg1 btavg2 btavg3 btavg4 btavg5 btavg6 btavg7 btavg8 btavg9];
btavgAnomFinal = [];
for ii = 1 : length(newLatGrid)-1
  str = ['btavgAnomFinal = [btavgAnomFinal btavg_newlat' num2str(ii) '];'];
  eval(str);
end

disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
for iCnt = 1 : length(newLatGrid)-1
  str = ['junk = btavg_newlat' num2str(iCnt) ';'];
  eval(str);
  figure(1); clf; pcolor(2002+daysSince2002A/365,h.vchan(ind),junk); colorbar; caxis([-1 +1]*2); shading flat; ylim([640 1640]); title(num2str(iCnt)); colormap(usa2);
  figure(2); clf; plot(h.vchan(ind),nanmean(junk,2)); xlim([640 1640]); ylabel('Mean Anomaly'); ylim([-1 +1]); 
    title(['Mean Anomaly for WideLatbin ' num2str(iCnt)]); plotaxis2;
  for ii = 1 : length(ind)
    P = polyfit(daysSince2002A,junk(ii,:),1);
    quicktrend(ii) = P(1)*365;
  end
  figure(3); clf; plot(h.vchan(ind),quicktrend,'r'); xlim([640 1640]); ylim([-1 +1]/10); ylabel('BT trend'); title(num2str(iCnt)); plotaxis2;
    title(['Mean BT Trend for WideLatbin ' num2str(iCnt)]); plotaxis2;
  trendavg(iCnt,:) = quicktrend;
  anomavg(iCnt,:) = nanmean(junk,2);
  pause(1);
end

figure(4); clf; pcolor(h.vchan(ind),meanvaluebin(newLatGrid),anomavg); colormap(usa2); colorbar; caxis([-1 +1]); xlim([640 1640]); shading flat; 
  title('mean spectral anomaly')
figure(5); clf; pcolor(h.vchan(ind),meanvaluebin(newLatGrid),anomavg/(iNumYears/2)); colormap(usa2); colorbar; caxis([-1 +1]/(iNumYears/2)); xlim([640 1640]); shading flat; 
  title(['mean of mean spectral anomaly'])
figure(6); clf; pcolor(h.vchan(ind),meanvaluebin(newLatGrid),trendavg); colormap(usa2); colorbar; caxis([-1 +1]*0.15); xlim([640 1640]); shading flat; 
  title('spectral trends')

iCnt = 1;
  %str = ['junk = btavg' num2str(iCnt) ';'];
  %eval(str);
  junk = btavg_cos;
  figure(1); clf; pcolor(2002+daysSince2002A/365,h.vchan(ind),junk); colorbar; caxis([-1 +1]*2); shading flat; ylim([640 1640]); title('Cosine Avg Anom'); colormap(usa2);
  figure(2); clf; pcolor(2002+daysSince2002A/365,h.vchan(ind),junk); colorbar; caxis([-1 +1]*1); shading flat; ylim([780  980]); title('Cosine Avg Anom'); colormap(usa2);
  figure(3); clf; plot(2002+daysSince2002A/365,junk);                title('Cosine Avg Anom');   xlim([min(2002+daysSince2002A/365) max(2002+daysSince2002A/365)]); 

  i900 = find(h.vchan >= 900,1); 
  figure(4); plot(2002+daysSince2002A/365,smooth(junk(i900,:),3)); plotaxis2;title('Cosine Avg Anom 900 cm-1');   xlim([min(2002+daysSince2002A/365) max(2002+daysSince2002A/365)]); 

  %% see list from Mitchell_Goldberg-Dissertation.pdf at beginning of driver_put_together_QuantileChoose_anomalies.m');
  iboo = 0; 
  iboo = iboo + 1; ijunk = find(h.vchan >= 0667.019,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0681.457,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0704.436,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0723.029,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0801.099,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0900.000,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 0960.000,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 1040.083,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 1231.300,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 1419.100,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 1519.070,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(h.vchan >= 1598.490,1); iaChan(iboo) = ijunk;
  figure(5); plot(2002+daysSince2002A/365,junk(iaChan,:)); plotaxis2;title('Cosine Avg Anom at some channels');   xlim([min(2002+daysSince2002A/365) max(2002+daysSince2002A/365)]); 
  figure(5); clf; for ii = 1 : length(iaChan); plot(2002+daysSince2002A/365,smooth(junk(iaChan(ii),:),3),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan)),'location','best','fontsize',8);

  figure(5); clf; baboo = find(h.vchan(iaChan) <= 802); for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('Temperature Sounding');
  figure(6); clf; baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('Window Channels');
  figure(7); clf; baboo = find(h.vchan(iaChan) > 1400); for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),5),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('WV Sounding');
  figure(8); clf; baboo = find(h.vchan(iaChan) <= 802); for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('Temperature Sounding'); xlim([2020 2025])
  figure(9); clf; baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('Window Channels'); xlim([2020 2025])
  figure(10); clf; baboo = find(h.vchan(iaChan) > 1400); for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,smooth(junk(iaChan(baboo(ii)),:),5),'linewidth',2); hold on; end; hold off
    plotaxis2; hl = legend(num2str(h.vchan(iaChan(baboo))),'location','best','fontsize',8); title('WV Sounding'); xlim([2020 2025])

  for ii = 1 : length(ind)
    P = polyfit(daysSince2002A,junk(ii,:),1);
    quicktrend(ii) = P(1)*365;
  end
  figure(11); clf; plot(h.vchan(ind),quicktrend); xlim([640 1640]); ylim([-1 +1]/10); ylabel('BT trend from Cosine Avg Anom'); title(num2str(iCnt)); plotaxis2;
  figure(12); clf; plot(h.vchan(ind),2*nanmean(junk,2)/iNumYears,'b',h.vchan(ind),quicktrend,'r'); 
    xlim([640 1640]); ylim([-1 +1]/10); plotaxis2; legend('Cosine Avg Anom/iNumYears/2','BT trend from Cosine Avg Anom','location','best','fontsize',10); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modis_ceres_era5_anoms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%iSave = input('save -1/+1  : ');
iSave = -1; 
disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  
if iSave > 0
  comment = 'see make_global_avg_and_N_averages_anomalies.m';
  fout = ['anomaly_globalavg_and_' num2str(length(newLatGrid)-1) '_averages_timeseries_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat'];
  saver = ['save -v7.3 ' fout ' btavgAnomFinal yy mm dd hh rtime comment usethese newLatGrid anomavg'];
  eval(saver);
end
%}
