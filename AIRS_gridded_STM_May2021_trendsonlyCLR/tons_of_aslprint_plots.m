
addpath /asl/matlib/plotutils
figure(5); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/cfc12_trends_global.pdf');
figure(27); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/clrOLR_trends_global.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%
load /home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat

cmax = 0.1;
cmax = 0.15
cmap = usa2;
cmap = llsmap5;

load umbc_trends.mat
figure(6); aslmap(6,rlat65,rlon73,smoothn(umbc_st_trend4608',1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt UMBC'); caxis([-cmax +cmax])
figure(6); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/st_trends_global.pdf');

figure(7); aslmap(7,rlat65,rlon73,smoothn(airs_quantile16_bt1231_trend4608',1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt BT1231'); caxis([-cmax +cmax])
figure(7); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/bt1231_trends_global.pdf');

%%%%%%%%%%
giss = load('giss_trends.mat');
figure(8); aslmap(8,rlat65,rlon73,smoothn(giss.giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt GISS'); caxis([-cmax +cmax])
figure(8); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/giss_trends_global.pdf');

plot(giss.giss_trend4608(:),results(:,6),'.'); xlabel('GISS SKT trend K/yr'); ylabel('UMBC SKT trend K/yr');
addpath /home/sergio/MATLABCODE/NANROUTINES
[r,chisqr,P,sigP,numpts] = nanlinearcorrelation(giss.giss_trend4608(:),results(:,6));

%%%%%%%%%%
airsL3 = load('airsL3_trends.mat');
figure(9); aslmap(9,rlat65,rlon73,smoothn(airsL3.airsL3_trend4608',1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt AIRSL3'); caxis([-cmax +cmax])
figure(9); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_trends_global.pdf');

plot(airsL3.airsL3_trend4608(:),results(:,6),'.'); xlabel('AIRSL3 SKT trend K/yr'); ylabel('UMBC SKT trend K/yr');
addpath /home/sergio/MATLABCODE/NANROUTINES
[r,chisqr,P,sigP,numpts] = nanlinearcorrelation(airsL3.airsL3_trend4608(:),results(:,6));

%%%%%%%%%%
figure(7); aslmap(7,rlat65,rlon73,smoothn(airs_quantile16_bt1231_trend4608'-umbc_st_trend4608',1), [-90 +90],[-180 +180]);  
  colormap(cmap);  title('d/dt (BT1231-UMBC)'); caxis([-cmax +cmax]/5)
  figure(7); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/bt1231-umbc_st_global.pdf');

figure(8); aslmap(8,rlat65,rlon73,smoothn(giss.giss_trend4608'-(reshape(results(:,6),72,64)'),1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt GISS-UMBC'); 
caxis([-cmax +cmax])
figure(8); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/giss_VS_umbc_trends_global.pdf');

figure(9); aslmap(9,rlat65,rlon73,smoothn(airsL3.airsL3_trend4608'-(reshape(results(:,6),72,64)'),1), [-90 +90],[-180 +180]);  colormap(cmap);  title('d/dt AIRSL3-UMBC'); 
caxis([-cmax +cmax])
figure(9); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_VS_umbc_trends_global.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%

%% then go to eg find_RH_trends 
figure(25); clf
pcolor(rlat,pavgLAY(1:97,1000),deltaTlat(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaTlat)); junk = cos(rlat) * ones(1,101);
area_wgtT = nansum(deltaTlat.*junk,1)./nansum(junk,1);
hold on; plot(area_wgtT(1:97)*1000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([1 1000]); caxis([-0.1 +0.1]); colorbar; plotaxis2;

figure(26); clf
pcolor(rlat,pavgLAY(1:97,1000),fracWVlat(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracWVlat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracWV = nansum(fracWVlat.*junk,1)./nansum(junk,1);
hold on; plot(area_wgt_fracWV(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-2 +2]*1e-3); colorbar; plotaxis2;

figure(27); clf
pcolor(rlat,pavgLAY(1:97,1000),deltaRHlat(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaRHlat)); junk = cos(rlat) * ones(1,100);
area_wgtRH = nansum(deltaRHlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtRH(1:97)*100,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-0.5 +0.5]); colorbar; plotaxis2;

%% and to find_O3_trends
figure(28); clf
pcolor(rlat,pavgLAY(1:97,1000),fracO3lat(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracO3lat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracO3 = nansum(fracO3lat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgt_fracO3(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([1 100]); caxis([-2 +2]*1e-3); colorbar; plotaxis2; colormap(usa2)

figure(29); clf
pcolor(rlat,pavgLAY(1:97,1000),deltaO3lat(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaO3lat)); junk = cos(rlat) * ones(1,97);
area_wgtO3 = nansum(deltaO3lat.*junk,1)./nansum(junk,1);
hold on; plot(area_wgtO3(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([1 100]); caxis([-0.5 +0.5]*1e-2); colorbar; plotaxis2; colormap(usa2)

for ii=25:29; figure(ii); colormap(usa2); xlabel('Latitude'); ylabel('p(mb)'); end
figure(25); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/T_z_trends_zonal.pdf');
figure(26); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/fracWV_z_trends_zonal.pdf');
figure(27); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/RH_z_trends_zonal.pdf');
figure(28); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/fracO3_z_trends_zonal.pdf');
figure(29); aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/ppmO3_z_trends_zonal.pdf');
