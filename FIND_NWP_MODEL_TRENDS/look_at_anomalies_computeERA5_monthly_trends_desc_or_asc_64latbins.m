addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE_Git/COLORMAP/
addpath /home/sergio/MATLABCODE_Git/COLORMAP/COLORBREWER/cbrewer/cbrewer

cmap = usa2;
cmap = cbrewer('div', 'RdBu', 128);
cmap = cbrewer('div', 'RdBu', 16);
cmap(cmap < 0) = 0; cmap(cmap > 1) = 1;
cmap = flipud(cmap);
jett = jet(128); jett(1,:) = 1;

load ERA5_atm_data_2002_09_to_2024_08_trends_desc_64latbins.mat

junkcnt = 0;
for yyjunk = StartY : StopY
  mmS = 1; mmE = 12;
  if yyjunk == StartY
    mmS = 9;
  elseif yyjunk == StopY
    mmE = 6;
    mmE = 8;
  end
  for mmjunk = mmS : mmE
    junkcnt = junkcnt + 1;
    yy(junkcnt) = yyjunk;
    mm(junkcnt) = mmjunk;
  end
end
yymm = yy + (mm-1)/12;

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = flipud(meanvaluebin(plevs));
hgt   = p2h(plays)/1000;
do_XX_YY_from_X_Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hunga Tonga
%% see https://www.nature.com/articles/s43247-024-01620-3#data-availability
%% Strong persistent cooling of the stratosphere after the Hunga eruption
%% Matthias Stocker, Andrea K. Steiner, Florian Ladstädter, Ulrich Foelsche & William J. Randel 
%% Communications Earth & Environment volume 5, Article number: 450 (2024) Cite this article

%% see ../AIRS_gridded_STM_May2021_trendsonlyCLR]/driver_gather_gridded_anomaly_retrieval_results.m
%% see ../AIRS_gridded_STM_May2021_trendsonlyCLR]/plot_anomalies_All.m
%% see ../AIRS_gridded_STM_May2021_trendsonlyCLR]/plot_Nature2024.m

figure(1); 
boo = find(rlat >= -45 & rlat <= 20);
wah = squeeze(nanmean(anom64_ptemp(:,boo,:),2));
pcolor(yymm,hgt,wah); shading interp; title('ERA5 2002/09-2024/08 T anomalies between -45 S and 20 N'); colorbar; colormap(cmap); caxis([-1 +1]*2)
xlim([2022-1/12 2022+2/12]); ylim([18 34]);

figure(2)
boo = find(rlat >= -30 & rlat <= 10);
wah = squeeze(nanmean(anom64_ptemp(:,boo,:),2));
pcolor(yymm,hgt,wah); shading interp; title('ERA5 2002/09-2024/08 T anomalies between -30S S and 10 N'); colorbar; colormap(cmap); caxis([-1 +1]*4)
xlim([2022-1/12 2024-1/12]); ylim([16 38]);

figure(3);
boo = find(hgt <= 19,1);
wah = squeeze(nanmean(anom64_ptemp(boo,:,:),1));
pcolor(yymm,rlat,wah); shading interp; title('ERA5 2002/09-2024/08 T anomalies at 19 km'); colorbar; colormap(cmap); caxis([-1 +1]*4)
xlim([2022-1/12 2024-1/12]); ylim([-1 +1]*45)

figure(4);
boo = find(hgt <= 27,1);
wah = squeeze(nanmean(anom64_ptemp(boo,:,:),1));
pcolor(yymm,rlat,wah); shading interp; title('ERA5 2002/09-2024/08 T anomalies at 27 km'); colorbar; colormap(cmap); caxis([-1 +1]*4)
xlim([2022-1/12 2024-1/12]); ylim([-1 +1]*45)

figure(5);
boo = find(hgt <= 32,1);
wah = squeeze(nanmean(anom64_ptemp(boo,:,:),1));
pcolor(yymm,rlat,wah); shading interp; title('ERA5 2002/09-2024/08 T anomalies at 32 km'); colorbar; colormap(cmap); caxis([-1 +1]*4)
xlim([2022-1/12 2024-1/12]); ylim([-1 +1]*45)

figure(6)
boo = find(rlat >= -30 & rlat <= 13);
wah = squeeze(nanmean(anom64_gas_1(:,boo,:),2));
pcolor(yymm,hgt,wah); shading interp; title('ERA5 2002/09-2024/08 WV anomalies between -30S S and +30 N'); colorbar; colormap(cmap); caxis([-1 +1]*0.25)
xlim([2022-1/12 2024-1/12]); ylim([16 38]);
xlim([2022-1/12 2025-1/12]); ylim([10 25]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE_Git/COLORMAP/LLS
load llsmap5

figure(7)
wah = squeeze(nanmean(anom64_ptemp,3));
wah = trend64_ptemp;
pcolor(rlat,plays,wah); shading interp; title('ERA5 2002/09-2024/08 T trends');colormap(llsmap5); caxis([-1 +1]*0.15); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')
set(gca,'yscale','log')

figure(8)
wah = squeeze(nanmean(anom64_gas_1,3));
wah = trend64_gas_1;
pcolor(rlat,plays,wah); shading interp; title('ERA5 2002/09-2024/08 WV trends');colormap(llsmap5); caxis([-1 +1]*0.015); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')

figure(9)
wah = squeeze(nanmean(anom64_gas_3,3));
wah = trend64_gas_3;
pcolor(rlat,plays,wah); shading interp; title('ERA5 2002/09-2024/08 O3 trends');colormap(llsmap5); caxis([-1 +1]*0.015); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')
set(gca,'yscale','log')

figure(10)
wah = squeeze(nanmean(anom64_RH,3));
wah = trend64_RH;
pcolor(rlat,plays,wah); shading interp; title('ERA5 2002/09-2024/08 RH trends');colormap(llsmap5); caxis([-1 +1]*0.5); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11)
wah = squeeze(nanmean(anom64_ptemp,2));
pcolor(yymm,plays,wah); shading interp; title('ERA5 2002/09-2024/08 T anoms');colormap(llsmap5); caxis([-1 +1]*2); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')
set(gca,'yscale','log')

figure(12)
wah = squeeze(nanmean(anom64_gas_1,2));
pcolor(yymm,plays,wah); shading interp; title('ERA5 2002/09-2024/08 WV anoms');colormap(llsmap5); caxis([-1 +1]*0.15); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')

figure(11)
wah = squeeze(nanmean(anom64_ptemp,2));
pcolor(yymm,plays,smoothn(wah,1)); shading interp; title('ERA5 2002/09-2024/08 T anoms');colormap(llsmap5); caxis([-1 +1]*2); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')
set(gca,'yscale','log')

figure(12)
wah = squeeze(nanmean(anom64_gas_1,2));
pcolor(yymm,plays,smoothn(wah,1)); shading interp; title('ERA5 2002/09-2024/08 WV anoms');colormap(llsmap5); caxis([-1 +1]*0.15); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')

tropics = find(abs(rlat) < 30);
figure(13)
wah = squeeze(nanmean(anom64_ptemp(:,tropics,:),2));
pcolor(yymm,plays,smoothn(wah,1)); shading interp; title('ERA5 TROPICS 2002/09-2024/08 T anoms');colormap(llsmap5); caxis([-1 +1]*2); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')
set(gca,'yscale','log')

figure(14)
wah = squeeze(nanmean(anom64_gas_1(:,tropics,:),2));
pcolor(yymm,plays,smoothn(wah,1)); shading interp; title('ERA5 TROPICS 2002/09-2024/08 WV anoms');colormap(llsmap5); caxis([-1 +1]*0.15); colorbar
ylim([10 1000]); set(gca,'ydir','reverse')

for ii = 1 : 14
  figure(ii); set(gca,'fontsize',12);
end
figure(12); ylim([50 1000]); set(gca,'yscale','log')
figure(14); ylim([50 1000]); set(gca,'yscale','log')
