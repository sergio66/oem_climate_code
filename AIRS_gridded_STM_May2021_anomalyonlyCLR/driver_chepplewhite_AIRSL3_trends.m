addpath /asl/matlib/maps
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load /home/chepplew/data/rates_anomalies/tiled/airs_l3/airs_l3_v7_lin_rates_interp.mat

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

load('llsmap5.mat');

plevs = flipud(plevs);

for ii = 1 : 64
  t_trends_zonal(ii,:) = nanmean(squeeze(infit.tz_d_coef(:,ii,:)),2);
end

figure(2); clf; pcolor(rlat,plevs,t_trends_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar('horizontal')
  title('Zonal d/dt T AIRS L3')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_T_zonal_trends.pdf');

i500 = find(plevs >= 500,1);
[plevs(i500-1) plevs(i500)]
%%%i500 = 6
T500 = squeeze(0.5*(infit.tz_d_coef(i500-1,:,:) + infit.tz_d_coef(i500,:,:)));
T500 = squeeze(infit.tz_d_coef(i500,:,:));
figure(1); aslmap(1,rlat65,rlon73,smoothn(T500,1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]/20)
title('d/dt AIRS L3 T(500 mb)');
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_rh_500mb_global_trends.pdf');

%{
clear data dataMap
data = smoothn(t_trends_zonal',1); 
dataMap = smoothn(T500,1);
save airsL3_T_zonal_trends.mat rlat plevs data dataMap rlat65 rlon73
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load /home/chepplew/data/rates_anomalies/tiled/airs_l3/airs_l3_v7_rh_lin_rates_interp.mat

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

load('llsmap5.mat');

for ii = 1 : 64
  rh_trends_zonal(ii,:) = nanmean(squeeze(infit.rh_d_coef(:,ii,:)),2);
end

hlevs = flipud(hlevs);

figure(2); clf; pcolor(rlat,hlevs,smoothn(rh_trends_zonal',1)); colormap(llsmap5); caxis([-0.5 +0.5])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar('horizontal')
  title('Zonal d/dt RH AIRS L3')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_rh_zonal_trends.pdf');

cmap = llsmap5; 

i600 = find(hlevs >= 600,1);
[hlevs(i600-1) hlevs(i600)]
rh600 = squeeze(0.5*(infit.rh_d_coef(i600-1,:,:) + infit.rh_d_coef(i600,:,:)));
figure(1); aslmap(1,rlat65,rlon73,smoothn(rh600,1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.5 +0.5])
title('d/dt AIRS L3 RH(600 mb)');
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_rh_600mb_global_trends.pdf');

i500 = find(hlevs >= 500,1);
[hlevs(i500-1) hlevs(i500)]
rh500 = squeeze(0.5*(infit.rh_d_coef(i500-1,:,:) + infit.rh_d_coef(i500,:,:)));
figure(1); aslmap(1,rlat65,rlon73,smoothn(rh500,1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.5 +0.5])
title('d/dt AIRS L3 RH(500 mb)');
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_rh_500mb_global_trends.pdf');

%{
clear data dataMap
data = smoothn(rh_trends_zonal',1); 
dataMap = smoothn(rh500,1);
save airsL3_rh_zonal_trends.mat rlat hlevs data dataMap rlat65 rlon73
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
load /home/chepplew/data/rates_anomalies/tiled/airs_l3/airs_l3_v7_rh_o3ppm_lin_rates_interp.mat

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

load('llsmap5.mat');

for ii = 1 : 64
  o3ppm_trends_zonal(ii,:) = nanmean(squeeze(infit.o3ppm_d_coef(:,ii,:)),2);
end

cmap = llsmap5; 

plevs = plevs;         %%% CHRIS WENT AND FIXED THE UPSIDE DOWN PROBLEM HERE UGH
i025 = find(plevs <= 25,1);

%plevs = flipud(plevs); %%% CHRIS DID not fix the problem
%i025 = find(plevs >= 25,1);

[plevs(i025-1) plevs(i025)]
o3ppm025 = squeeze(0.5*(infit.o3ppm_d_coef(i025-1,:,:) + infit.o3ppm_d_coef(i025,:,:)));
figure(1); aslmap(1,rlat65,rlon73,smoothn(o3ppm025,1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]/100)
title('d/dt AIRS L3 O3(025 mb)');
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_o3ppm_025mb_global_trends.pdf');

%figure(1); clf; pcolor(rlat,plevs,       (o3ppm_trends_zonal')); colormap(llsmap5); caxis([-0.005 +0.005]*2)
%  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([1 100]); colorbar('horizontal')
figure(2); clf; pcolor(rlat,plevs,smoothn(o3ppm_trends_zonal',1)); colormap(llsmap5); caxis([-0.005 +0.005]*2)
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([1 100]); colorbar('horizontal')
  title('Zonal d/dt O3 AIRS L3')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/airsL3_o3_zonal_trends.pdf');

%{
clear data dataMap
data = smoothn(o3ppm_trends_zonal',1); 
dataMap = smoothn(o3ppm025,1);
save airsL3_o3_zonal_trends.mat rlat plevs data dataMap rlat65 rlon73
%}
