addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/maps

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /home/chepplew/data/rates_anomalies/tiled/merra2/merra2_rh_trend.mat

plevs = load('merra2_rhlevs42.dat'); plevs = plevs(:,2);
i050 = find(plevs >= 050,1);
i500 = find(plevs >= 500,1);
cmap = llsmap5;

cmax = 0.5;
figure(1); clf
aslmapSergio(rlat65,rlon73,100*smoothn(squeeze(fit.rha_sub_coef(:,:,i500,2)'),1), [-90 +90],[-180 +180]);
colormap(cmap); colorbar('south'); caxis([-cmax +cmax]);
title('MERRA RH 500 mb from Chris')

figure(2); clf
wah = 100*squeeze(fit.rha_sub_coef(:,:,:,2));
bah = smoothn(squeeze(nanmean(wah,1)),1);
pcolor(rlat,plevs,bah');
cmap = llsmap5;
set(gca,'ydir','reverse');
set(gca,'yscale','log');
set(gca,'colormap',cmap,'CLim',[-cmax +cmax]);
set(gca,'Ylim',[100 1000]);
shading interp
colorbar
title('MERRA RH from Chris')

%{
clear data dataMap
data    = bah;
dataMap = 100*smoothn(squeeze(fit.rha_sub_coef(:,:,i500,2)'),1);
plevsM = plevs;
save merra_rh_zonal_trends.mat rlat plevsM data dataMap rlat65 rlon73 
%}

disp('ret to cont'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /home/chepplew/data/rates_anomalies/tiled/merra2/merra2_tz_trend.mat

plevs = load('merra2_tlevs72.dat'); plevs = plevs(:,2);
i050 = find(plevs >= 050,1);
i500 = find(plevs >= 500,1);
cmap = llsmap5;

cmax = 0.15;
figure(1); clf
aslmapSergio(rlat65,rlon73,smoothn(squeeze(fit.tza_sub_coef(:,:,i500,2)'),1), [-90 +90],[-180 +180]);
colormap(cmap); colorbar('south'); caxis([-cmax +cmax]);
title('MERRA TZ 500 mb from Chris')

figure(2); clf
wah = squeeze(fit.tza_sub_coef(:,:,:,2));
bah = smoothn(squeeze(nanmean(wah,1)),1);
pcolor(rlat,plevs,bah');
cmap = llsmap5;
set(gca,'ydir','reverse');
set(gca,'yscale','log');
set(gca,'colormap',cmap,'CLim',[-cmax +cmax]);
set(gca,'Ylim',[10 1000]);
shading interp
colorbar
title('MERRA TZ from Chris')

%{
clear data dataMap
data    = bah;
dataMap = smoothn(squeeze(fit.tza_sub_coef(:,:,i500,2)'),1);
plevsM = plevs;
save merra_T_zonal_trends.mat rlat plevsM data dataMap rlat65 rlon73 
%}

disp('ret to cont'); pause;
error('kjs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmax = 0.1;
figure(1); clf
aslmapSergio(rlat65,rlon73,smoothn(reshape(squeeze(fits.o3pp_coef(i050,:,2)),64,72),1), [-90 +90],[-180 +180]);
colormap(cmap); colorbar('south'); caxis([-cmax +cmax]);
title('ERA O3PPM 050 mb from Chris')

figure(2); clf
wah = reshape(squeeze(fits.o3pp_coef(:,:,2)),98,72,64);
bah = smoothn(squeeze(nanmean(wah,2)),1);
pcolor(rlat,play(1:97),bah(1:97,:));
cmap = llsmap5;
set(gca,'ydir','reverse');
set(gca,'yscale','log');
set(gca,'colormap',cmap,'CLim',[-cmax +cmax]);
set(gca,'Ylim',[1 100]);
shading interp
colorbar
title('ERA O3PPM from Chris')

disp('ret to cont'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
