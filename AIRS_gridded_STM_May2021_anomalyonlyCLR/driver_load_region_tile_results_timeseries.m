addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/matlib/science/
addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools

%clear all

%for ii = 1 : 10
%  figure(ii); colormap jet
%end

iNumLay = 33;
iNumLay = 20;

if ~exist('iQuantile')
  iQuantile = 16;
  iQuantile = 08;
  iQuantile = input('Enter (08) 50% or (16) hottest : ');
  if length(iQuantile) == 0
    iQuantile = 16;
  end
end

if ~exist('p')
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
  ctmp = load('/home/sergio/MATLABCODE/PLOTTER/coast.mat');
end

if ~exist('region')
  %% see driver_select_tiles_for_clust_quickBTanom.m
  disp('Region (0) All');
  disp('       (1) India');
  disp('       (2) ENSO 3.4');
  disp('       (3) Artic Oscillation lat -90 to +90, around Greenland');
  disp('       (4) Southern Oscillation lon -180 to 180, lat around -60 S : ');
  disp('       (5) WV taprecorder lon -180 to 180, lat around +/- 8 EQUATOR : ');
  disp('       (6) ENSO 3.0');
  iaTile = input('Enter Region : ');
  if iaTile == 0
    region = 1 : 4608;
  elseif iaTile == 1
    %% india
    region = find(p.rlat >= +5 & p.rlat <= +35 & p.rlon >= +70 & p.rlon <= +90);
  elseif iaTile == 2
    %% Nino SST Region 3.4
    %% https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni
    region = find(p.rlat >= -5 & p.rlat <= +5 & p.rlon >= -170 & p.rlon <= -120);
  elseif iaTile == 3
    %% artic oscillation
    region = find(p.rlat >= -90 & p.rlat <= +90 & p.rlon >= -60 & p.rlon <= -20);
  elseif iaTile == 4
    %% southern ocean
    region = find(p.rlat >= -75 & p.rlat <= -55 & p.rlon >= -180 & p.rlon <= +180);
  elseif iaTile == 5
    %% wv taperecorder around equator
    region = find(p.rlat >= -8 & p.rlat <= +8 & p.rlon >= -180 & p.rlon <= +180);
  elseif iaTile == 6
    %% Nino SST Region 3.0
    %% https://journals.ametsoc.org/view/journals/clim/15/18/1520-0442_2002_015_2616_tttvcb_2.0.co_2.xml
    region = find(p.rlat >= -5 & p.rlat <= +5 & p.rlon >= -150 & p.rlon <= -90);
  end
end

if ~exist('results')
  results = zeros(length(region),412,6);
  resultsT  = zeros(length(region),412,20);
  resultsWV = zeros(length(region),412,20);
  resultsO3 = zeros(length(region),412,20);
end

if ~exist('iaFound')
  iaFound = zeros(length(region),412);
end

figure(1); clf;
whos region
plot(p.rlon,p.rlat,'.',p.rlon(region),p.rlat(region),'rx')
hold on; plot(ctmp.long, ctmp.lat, 'k','linewidth',2); hold off

iJ0 = input('Start reading tile number X .... (the defult is 1) : ');
if length(iJ0) == 0
  iJ0 = 1;
end
for iJ = iJ0 : length(region)
  iTile = region(iJ);
  YY = floor((iTile-1)/72) + 1;
  XX = iTile-(YY-1)*72;
  [iTile XX YY]

  data_trends = load(['ANOM/LatBin' num2str(YY,'%02d') '/LonBin' num2str(XX,'%02d') '/v1_Quantile16/desc_anaomaly_LatBin' num2str(YY,'%02d') '_LonBin' num2str(XX,'%02d') '_Quantile16_timesetps_001_412_V1.mat']);
  obs_anom(iJ,:,:) = data_trends.bt_anom_all;

  iWarning = 0;
  for ii = 1 : 412
    fname = ['OutputAnomaly_OBS/' num2str(iTile,'%04d') '/anomtest_timestep' num2str(ii) '.mat'];
    fname = ['OutputAnomaly_OBS/' num2str(iTile) '/anomtest_timestep' num2str(ii) '.mat'];
    fname = ['OutputAnomaly_OBS/Quantile'  num2str(iQuantile,'%02d') '/' num2str(iTile) '/anomtest_timestep' num2str(ii) '.mat'];
    fname = ['OutputAnomaly_OBS_May12_2021_AIRS_STM//Quantile'  num2str(iQuantile,'%02d') '/' num2str(iTile) '/anomtest_timestep' num2str(ii) '.mat'];
    if exist(fname) & iaFound(iJ,ii) == 1
      fprintf(1,'%s already loaded in \n',fname);
    elseif exist(fname) & iaFound(iJ,ii) == 0
      iExist = +1;
      loader = ['load ' fname];
      eval(loader);
      iaFound(iJ,ii) = 1;
      results(iJ,ii,1:6) = oem.finalrates(1:6);
      [mmn,nn] = size(oem.ak_water);
  
      nlays(ii) = nn;
      nn0 = min(nn,iNumLay);
      if nn0 == nn
        resultsWV(iJ,ii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
        resultsT(iJ,ii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
        resultsO3(iJ,ii,1:nn) = oem.finalrates((1:nn)+6+nn*2);
      else
        iWarning = iWarning + 1
        iaWarning(iWarning) = ii;
        wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(iJ,ii,1:nn0) = wah(1:nn0);
        wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(iJ,ii,1:nn0) = wah(1:nn0);
        wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(iJ,ii,1:nn0) = wah(1:nn0);
      end
 
      rates(iJ,ii,:) = rateset.rates;
      fits(iJ,ii,:)  = oem.fit';
    else
      iExist = -1;
      iaFound(iJ,ii) = 0;
      results(iJ,ii,1:6) = NaN;
      resultsWV(iJ,ii,:) = NaN;
      resultsT(iJ,ii,:)  = NaN;
      resultsO3(iJ,ii,:) = NaN;
      rates(iJ,ii,:) = NaN;
      fits(iJ,ii,:) = NaN;
    end
    if mod(ii,72) == 0
      fprintf(1,'+ %2i\n',ii/72);
    elseif iExist == 1
      fprintf(1,'.');
    elseif iExist == -1
      fprintf(1,' ');
    end
  end
  fprintf(1,'\n');
  fprintf(1,'tile %3i of %3i tiles for region %2i : found %4i of %4i timesteps \n',iJ,length(region),iaTile,sum(iaFound(iJ,:)),412)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);
plevs = flipud(plevs);

clear pavg
for ii = 1 : iNumLay
  %iavg = jacobian.wvjaclays_used{iNumLay}-6;
  iavg = jacobian.wvjaclays_used{ii}-jacobian.wvjaclays_offset;
  pavg(ii) = mean(plays(iavg));
  pavg20(ii) = mean(plays(iavg));
  plevs20_top(ii) = min(plevs(iavg));
  plevs20_bot(ii) = max(plevs(iavg));
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
figure(1); scatter_coast(X(:),Y(:),50,landfrac); colorbar; title('landfrac');  caxis([0 1])

datime = (1:412);
datime = datime*16/365 + 2002.75;

figure(1); plot(datime,squeeze(nanmean(results(:,:,1))));  title('CO2');     xlim([min(datime) max(datime)])
figure(2); plot(datime,squeeze(nanmean(results(:,:,2))));  title('N2O');     xlim([min(datime) max(datime)])
figure(3); plot(datime,squeeze(nanmean(results(:,:,3))));  title('CH4');     xlim([min(datime) max(datime)])
figure(4); plot(datime,squeeze(nanmean(results(:,:,4))));  title('CFC11');   xlim([min(datime) max(datime)])
figure(5); plot(datime,squeeze(nanmean(results(:,:,5))));  title('CFC12');   xlim([min(datime) max(datime)])
figure(6); plot(datime,squeeze(nanmean(results(:,:,6))),datime,squeeze(nanmean(obs_anom(:,1520,:))));  
  plotaxis2; title('ST'); xlim([min(datime) max(datime)]); hl = legend('SST retrieve','BT1231 anom','location','best','fontsize',10);
iaUTlay = find(pavg > 100 & pavg < 200);
  wah = nanmean(nanmean(resultsT(:,:,iaUTlay),1),3); whos wah
  figure(6); plot(datime,squeeze(nanmean(results(:,:,6))),'b',datime,squeeze(nanmean(obs_anom(:,1520,:))),'r',datime,wah,'k','linewidth',2);  
  plotaxis2; title('ST and UT temp'); xlim([min(datime) max(datime)]); hl = legend('SST retrieve','BT1231 anom','UT Temp retrieve','location','best','fontsize',10);
iaLTlay = find(pavg > 200);
  wah = nanmean(nanmean(resultsT(:,:,iaLTlay),1),3); whos wah
  figure(6); plot(datime,squeeze(nanmean(results(:,:,6))),'b',datime,squeeze(nanmean(obs_anom(:,1520,:))),'r',datime,wah,'k','linewidth',2);  
    plotaxis2; title('ST and LT temp'); xlim([min(datime) max(datime)]); hl = legend('SST retrieve','BT1231 anom','LT Temp retrieve','location','best','fontsize',10);
  %% 23 timesteps/year ==> 0.5 yr ~ 11 steps
  figure(6); plot(datime,smooth(squeeze(nanmean(results(:,:,6))),11),'b',...
                  datime,smooth(squeeze(nanmean(obs_anom(:,1520,:))),11),'r',datime,smooth(wah,11),'k','linewidth',2);  
    plotaxis2; title('ST and LT temp'); xlim([min(datime) max(datime)]); hl = legend('SST retrieve','BT1231 anom','LT Temp retrieve','location','best','fontsize',10);
iaBDRlay = find(pavg > 800);
  wahWV = nanmean(nanmean(resultsWV(:,:,iaBDRlay),1),3); whos wah
  figure(6); plot(datime,smooth(squeeze(nanmean(results(:,:,6))),11),'b',...
                  datime,smooth(squeeze(nanmean(obs_anom(:,1520,:))),11),'r',datime,smooth(wah,11),'k',datime,smooth(wahWV,11)*10,'g','linewidth',2);  
    plotaxis2; title('ST and LT temp'); xlim([min(datime) max(datime)]); 
    hl = legend('SST retrieve','BT1231 anom','LT Temp retrieve','10*WV BDRY','location','best','fontsize',10);
 
figure(7); pcolor(datime,pavg,squeeze(nanmean(resultsT,1))');  title('T'); caxis([-5 +5]);
figure(8); pcolor(datime,pavg,squeeze(nanmean(resultsWV,1))'); title('WV'); xlim([min(datime) max(datime)]); colormap(usa2); caxis([-1 +1]);
figure(9); pcolor(datime,pavg,squeeze(nanmean(resultsO3,1))'); title('O3'); xlim([min(datime) max(datime)]); colormap(usa2); caxis([-1 +1])
for ii = 7 : 9; figure(ii); xlim([min(datime) max(datime)]);  colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; colorbar; end
figure(7); ylim([10 1000]);  caxis([-2 +2])
figure(8); ylim([100 1000]); caxis([-0.5 +0.5])
figure(9); ylim([1 100]);    caxis([-0.5 +0.5])
Xpt = X(:); Xpt = Xpt(region); Ypt = Y(:); Ypt = Ypt(region);
figure(10); scatter_coast(X(:),Y(:),50,landfrac); colorbar off; title('landfrac');  caxis([-10 +10]); colormap(usa2); hold on; plot(Xpt,Ypt,'kx','Linewidth',4,'Markersize',10); hold off

addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
figure(11);
for ii = 1 : length(region)
 data = squeeze(results(ii,:,6));
 [B, stats] = Math_tsfit_lin_robust_finitewrapper((1:412)*16,data,4);
 st_trend(ii)     = B(2);
 st_trend_err(ii) = stats.se(2);
end
figure(11); 
  scatter_coast(p.rlon(region),p.rlat(region),200,st_trend); colorbar; title('ST Trend K/yr')
  simplemap(p.rlat(region),p.rlon(region),st_trend,5); colorbar; title('ST Trend K/yr')
  axis([min(p.rlon(region))-5 max(p.rlon(region))+5 min(p.rlat(region))-5 max(p.rlat(region))+5])
  caxis([-0.05 +0.05]); colormap(usa2)
  caxis([-0.1 +0.1]); colormap(usa2)

if iaTile == 0 | iaTile == 3
  raLats = unique(Y(:));
  for ii = 1 : length(raLats)
    miaow = find(abs(p.rlat(region) - raLats(ii)) <= 1);
    num_found(ii) = length(miaow);
    mean_st(ii,:) = nanmean(squeeze(results(miaow,:,6)),1);
    mean_T(ii,:,:) = nanmean(squeeze(resultsT(miaow,:,1:20)),1);
    mean_WV(ii,:,:) = nanmean(squeeze(resultsWV(miaow,:,1:20)),1);
    mean_O3(ii,:,:) = nanmean(squeeze(resultsO3(miaow,:,1:20)),1);
  end

  figure(12);
  %% see https://acp.copernicus.org/articles/13/4563/2013/acp-13-4563-2013.pdf Fig 2
  iiTrop = find(abs(raLats) <= 08);
  iiTrop = find(abs(raLats) <= 30);
  pcolor(datime,1:20,(squeeze(mean_WV(iiTrop(1),:,:)))'); shading interp; colorbar; caxis([-0.5 +0.5])
  boo = squeeze(nanmean(mean_WV(iiTrop,:,:),1));
  for ii = 1 : 20
    smooth_mean_WVTrop(:,ii) = smooth(boo(:,ii),8);
  end
  pcolor(datime,pavg,smooth_mean_WVTrop'); shading interp; colorbar; caxis([-0.25 +0.25]); set(gca,'ydir','reverse'); set(gca,'yscale','log');
  ylim([10 300]); colormap(usa2)

  figure(13)
  %% see https://acp.copernicus.org/articles/13/4563/2013/acp-13-4563-2013.pdf Fig 3
  ii100 = find(abs(pavg-100) == min(abs(pavg-100)));
  pcolor(datime,raLats,squeeze(mean_WV(:,:,ii100))); shading interp; colorbar; caxis([-0.5 +0.5])
  boo = squeeze(mean_WV(:,:,ii100));
  for ii = 1 : 64
    smooth_mean_WV100mb(ii,:) = smooth(boo(ii,:),8);
    smooth_mean_WV100mb(ii,:) = smooth(boo(ii,:),16);
    smooth_mean_WV100mb(ii,:) = smooth(boo(ii,:),32);
    smooth_mean_WV100mb(ii,:) = smooth(boo(ii,:),4);
  end
  pcolor(datime,raLats,smooth_mean_WV100mb); shading interp; colorbar; caxis([-0.25 +0.25]); colormap(usa2)
  for ii = 2003:2020
    line([ii ii],[-90 +90],'color','k')
  end
end
error('lkjgd')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
addpath /asl/matlib/plotutils
figure(10); aslprint(['/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/region_map_location_' num2str(iaTile) '.pdf']);
figure(11); aslprint(['/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/region_map_st_trend_' num2str(iaTile) '.pdf']);
figure(06); aslprint(['/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/region_map_st_anomaly_' num2str(iaTile) '.pdf']);
figure(07); aslprint(['/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/region_map_T_z_anomaly_' num2str(iaTile) '.pdf']);
%}
