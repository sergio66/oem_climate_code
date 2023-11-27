dirCode = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/';
dirData = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/driver_bt1231_timeseries_YY0MM0DD0_YYEMMEDDE_0001_4608.m')
disp('see Mitchell_Goldberg-Dissertation.pdf for list of chans : wget https://aosc.umd.edu/sites/default/files/dissertations-theses/Mitchell%20Goldberg-Dissertation.pdf')

  disp('Fig 5.8        667.766 cm-1       1 mb')
  disp('Fig 5.9        667.715 cm-1       2 mb')
  disp('Fig 5.9        667.270 cm-1      15 mb')
  disp('Fig 5.8        667.018 cm-1      25 mb')
  disp('Fig 4.4        681.457 cm-1      90 mb')
  disp('Fig 4.5        704.436 cm-1     350 mb')
  disp('Fig 4.6        723.029 cm-1     700 mb')
  disp('Fig 4.7        801.099 cm-1     850 mb')
  disp(' ')
  disp('Fig 4.8        1519.07 cm-1     315 mb')
  disp('Fig 4.9        1598.49 cm-1     490 mb  MidTropWV')
  disp('Fig 5.10       1519.07 cm-1     315 mb  UTWV')
  disp('Fig 5.1        1520.87 cm-1     315 mb???')
  disp(' ')
  disp('Fig 5.4        1040.03 cm-1     80 mb???')
  disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iHuge = -1;
iHuge = +1;
iHuge = input('Enter (-1/default) save one channel, all tiles, all timesteps or (+1) save all channels, all tiles, all timesteps (huge memory) : ');
if length(iHuge) == 0
  iHuge = -1;
end

iNumYears = 20; iNumAnomTimeSteps = 454; % 385/16 = 22.8125; floor(iNumYears*22.8125 - 2) = 454!!!! (about 23 timesteps/year, 2 are missing because of missing AIRS data)
if iHuge == +1
  disp('wow : lotsa memory!')
  btanom = nan(4608,2645,454); %% too much memory!
  monitor_memory_whos
else
  btanom = nan(4608,454);
end
bt1231 = nan(1,4608);
btChID = nan(1,4608);

disp('reading in anomalies')

load h2645structure.mat

iCnt = 0;
iQuant = 3;
iChID  = 0213;  %% 0704 cm-1     find(h.vchan >= 0704,1)
iChID  = 1145;  %% 1040 cm-1     find(h.vchan >= 1040,1)
iChID  = 2025;  %% 1519 cm-1     find(h.vchan >= 1519,1)
iChID  = 1861;  %% 1419 cm-1     find(h.vchan >= 1419,1)
iChID  = 1520;  %% 1231 cm-1     find(h.vchan >= 1231,1)
iChID  = 1511;  %% 1226.82 cm-1  find(h.vchan >= 1226.65,1)

%{
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
[h,ha,p,pa] = rtpread('/asl/rtp/airs/airs_l1c_v674/clear/2023/ecmwf_airicrad_day022_clear.rtp');

% ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/tile_fits_quantiles.m
% for qi = qi1 : qi2
%    % r1231 = squeeze(d.rad_quantile_desc(:,1520,qi));
%    % r1228 = squeeze(d.rad_quantile_desc(:,1513,qi));
%    r1231 = squeeze(d.rad_quantile_desc(k_desc,1520,qi));
%    r1228 = squeeze(d.rad_quantile_desc(k_desc,1513,qi));
%    bt1231 = rad2bt(fairs(1520),r1231);
%    bt1228 = rad2bt(fairs(1513),r1228);
% end

boo = [1511 1520];

addpath /home/sergio/MATLABCODE
load /home/sergio/KCARTA/L2SComparisons/l2s_kc122_H16_605_2830.mat
[fc,qc] = convolve_airs(double(w),double(d),double(h.ichan));
semilogy(w,d(:,1),fc,qc(:,1),fc(boo),qc(boo,1),'ro'); xlim([1220 1240])

plot(p.rlon,p.rlat,'.')
lala = find(p.rlat >= 0,1);

robs = nanmean(p.robs1,2);
robs = p.robs1(:,lala);

plot(h.vchan,rad2bt(h.vchan,robs),h.vchan(boo),rad2bt(h.vchan(boo),robs(boo)),'ro'); xlim([1220 1235])

%}

if iHuge < 0
  junk = -1;
  while junk < 0
    junk = input('Enter wavenumber to look for : ');
    iChID = find(h.vchan >= junk,1);
    fprintf(1,'iChID = %4i corresponds to SARTA ID %4i %6.2f cm-1 \n',iChID,h.ichan(iChID),h.vchan(iChID))
    junk = input('Proceed (-1/+1[default]) : ');
    if length(junk) == 0
      junk = +1;
    end
    if junk < 0
      return
    end
  end
end

for jj = 1 : 64
  fprintf(1,'latbin %2i of 64 \n',jj);
  figure(1); clf; pcolor(reshape(bt1231,72,64)'); title('BT 1231'); shading interp; colorbar; pause(0.1);
  figure(2); clf; pcolor(reshape(btChID,72,64)'); title('BT ChID'); shading interp; colorbar; pause(0.1);
  for ii = 1 : 72
    iCnt = iCnt + 1;
    fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_3_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
    fname = [dirData fname '_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'];
    iaFound(iCnt) = 0;   
    if exist(fname)
      iaFound(iCnt) = +1;
      a = load(fname,'bt_anom_desc','bt_desc');
      junk = a.bt_anom_desc;
      if iHuge == +1
        btanom(iCnt,:,:) = squeeze(junk(iQuant,:,:));
      else
        btanom(iCnt,:) = squeeze(junk(iQuant,iChID,:));
      end
      bt1231(iCnt) = a.bt_desc(1520,iQuant);
      btChID(iCnt) = a.bt_desc(iChID,iQuant);
    end
  end
end
figure(1); clf; pcolor(reshape(bt1231,72,64)'); title('BT 1231'); colormap jet; shading interp; colorbar; pause(0.1);
figure(2); clf; pcolor(reshape(btChID,72,64)'); title('BT ChID'); colormap jet; shading interp; colorbar; pause(0.1);

a = load(fname,'yy_desc','mm_desc','dd_desc','hh_desc','rtime_desc');
yy = a.yy_desc;
mm = a.mm_desc;
dd = a.dd_desc;
hh = a.hh_desc;
rtime = a.rtime_desc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/COLORMAP

do_XX_YY_from_X_Y
[salti, landfrac] = usgs_deg10_dem(Y,X);
lf = landfrac; lf = lf(:);
[mmjunk,nnjunk] = size(XX);
if mmjunk == 1
  XX = XX';
  YY = YY';
end

figure(1); clf; scatter_coast(XX(:),YY(:),100,bt1231(:)); title('BT 1231'); colormap jet; shading interp; colorbar; pause(0.1);
figure(2); clf; scatter_coast(XX(:),YY(:),100,btChID(:)); title('BT ChID'); colormap jet; shading interp; colorbar; pause(0.1);

coslat = cos(YY*pi/180)*ones(1,454);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iHuge < 0
  figure(3);
  daysSince2002 = change2days(yy,mm,dd,2002);
  P = polyfit(daysSince2002,nanmean(btanom,1),1);         fprintf(1,'quick BTtrend = %8.6f K/year \n',P(1)*365)
  P = polyfit(daysSince2002,nanmean(btanom.*coslat,1),1); fprintf(1,'wgted BTtrend = %8.6f K/year \n',P(1)*365)
  
  pcolor(2002+daysSince2002/365,rlat,squeeze(nanmean(reshape(btanom,72,64,454),1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
    xlabel('Time'); ylabel('Latitude'); title(['BT' num2str(h.vchan(iChID),'%6.2f') ' cm-1 anomaly']); 
  hold on; plot(2002+daysSince2002/365,100*nanmean(btanom,1),'k','linewidth',2); hold off; plotaxis2;
  
  figure(4);
  oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
  plot(onidd,10*oni,'k',2002+daysSince2002/365,smooth(100*nanmean(btanom,1),10),'b','linewidth',2); plotaxis2; xlim([2002 2023]); ylim([-1 +1]*50)
  
  P = polyfit(daysSince2002,nanmean(btanom,1),1); YP = polyval(P,daysSince2002);
  PX = polyfit(daysSince2002,nanmean(btanom.*coslat,1),1); YPX = polyval(PX,daysSince2002);
  
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(btanom,1)),10),'b',...
       2002+daysSince2002/365,smooth(100*(nanmean(btanom,1)-YP),10),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(btanom,1)-YP),24),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(btanom.*coslat,1)-YPX),24),'r','linewidth',2); xlim([2002 2023]);
  plotaxis2;
  hl = legend('ONI',['BT' num2str(h.vchan(iChID),'%6.2f') ' cm-1 anomaly'],'location','best');
  hl = legend('ONI',['BT' num2str(h.vchan(iChID),'%6.2f')],'location','best');

  [U,S,V] = svd(btanom,'econ');

elseif iHuge > 0

  figure(3);
  daysSince2002 = change2days(yy,mm,dd,2002);
  moo = squeeze(nanmean(btanom,1));
  P = polyfit(daysSince2002,moo(1520,:),1);         fprintf(1,'quick BTtrend = %8.6f K/year \n',P(1)*365)
  moox = squeeze(btanom(:,1520,:));
  P = polyfit(daysSince2002,nanmean(moox.*coslat,1),1); fprintf(1,'wgted BTtrend = %8.6f K/year \n',P(1)*365)
  
  pcolor(2002+daysSince2002/365,rlat,squeeze(nanmean(reshape(moox,72,64,454),1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
    xlabel('Time'); ylabel('Latitude'); title(['BT' num2str(h.vchan(iChID),'%6.2f') ' cm-1 anomaly']); 
  hold on; plot(2002+daysSince2002/365,100*nanmean(moox,1),'k','linewidth',2); hold off; plotaxis2;
  
  figure(4);
  oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
  plot(onidd,10*oni,'k',2002+daysSince2002/365,smooth(100*nanmean(moox,1),10),'b','linewidth',2); xlim([2002 2023]);
  
  P = polyfit(daysSince2002,nanmean(moox,1),1); YP = polyval(P,daysSince2002);
  PX = polyfit(daysSince2002,nanmean(moox.*coslat,1),1); YPX = polyval(PX,daysSince2002);
  
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(moox,1)),10),'b',...
       2002+daysSince2002/365,smooth(100*(nanmean(moox,1)-YP),10),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(moox,1)-YP),24),'r','linewidth',2); xlim([2002 2023]);
  plot(onidd,10*oni,'k',...
       2002+daysSince2002/365,smooth(100*(nanmean(moox.*coslat,1)-YPX),24),'r','linewidth',2); xlim([2002 2023]);
  plotaxis2;
  hl = legend('ONI',['BT' num2str(h.vchan(iChID),'%6.2f') ' cm-1 anomaly'],'location','best');
  hl = legend('ONI',['BT' num2str(h.vchan(iChID),'%6.2f')],'location','best');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSave = input('save (-1/+1) : ');
if iSave > 0
  comment = 'look at driver_put_together_QuantileChoose_anomalies.m';
  if iHuge < 0
    fout = ['anomaly_chID_' num2str(iChID,'%04d') '_Q' num2str(iQuant,'%02d') '.mat'];
    saver = ['save ' fout ' btanom yy mm dd hh rtime comment bt1231 btChID'];
  else
    fout = ['/asl/s1/sergio/JUNK/anomaly_ALL_Q' num2str(iQuant,'%02d') '.mat'];
    saver = ['save -v7.3 ' fout ' btanom yy mm dd hh rtime comment bt1231'];
  end
  eval(saver)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iHuge > 0
  iX = input('Do you wish to do global and (+1/default) N zonal averages or (-1) one tile : ');
  iX = 1;
  if lenth(iX) == 0
    iX = 1;
  end
  if iX > 0
    make_globalavg_and_N_average_anomalies
  else
    make_globalavg_and_onetile_anomaly
    %make_onetile_anomaly
  end
end
