addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/maps

load /home/chepplew/data/rates_anomalies/tiled/era_interim/erai_rh_tz_03ppm_rates_interp.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% plev = flipud(plev); apparently uses correct plev

playN = plev(1:100)-plev(2:101);
playD = log(plev(1:100)./plev(2:101));
play = playN./playD;

i050 = find(plev >= 050,1);
i500 = find(plev >= 500,1);
cmap = llsmap5;

i050 = find(play >= 050,1);
i500 = find(play >= 500,1);
cmap = llsmap5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmax = 0.5;
figure(1); clf
aslmapSergio(rlat65,rlon73,smoothn(reshape(squeeze(fits.rh_coef(i500,:,2)),64,72),1), [-90 +90],[-180 +180]);
colormap(cmap); colorbar('south'); caxis([-cmax +cmax]);
title('ERA RH 500 mb from Chris')

figure(2); clf
wah = reshape(squeeze(fits.rh_coef(:,:,2)),98,72,64);
bah = smoothn(squeeze(nanmean(wah,2)),1);
pcolor(rlat,play(1:97),bah(1:97,:));
cmap = llsmap5;
set(gca,'ydir','reverse');
set(gca,'yscale','log');
set(gca,'colormap',cmap,'CLim',[-cmax +cmax]);
set(gca,'Ylim',[100 1000]);
shading interp
colorbar
title('ERA RH from Chris')

%{
clear data dataMap
data    = bah(1:97,:);
dataMap = smoothn(reshape(squeeze(fits.rh_coef(i500,:,2)),64,72),1);
play97 = play(1:97);
save era_rh_zonal_trends_donebyCLH.mat rlat play97 data dataMap rlat65 rlon73 
%}

disp('ret to cont'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmax = 0.15;
figure(1); clf
aslmapSergio(rlat65,rlon73,smoothn(reshape(squeeze(fits.tz_coef(i500,:,2)),64,72),1), [-90 +90],[-180 +180]);
colormap(cmap); colorbar('south'); caxis([-cmax +cmax]);
title('ERA TZ 500 mb from Chris')

figure(2); clf
wah = reshape(squeeze(fits.tz_coef(:,:,2)),99,72,64);
bah = smoothn(squeeze(nanmean(wah,2)),1);
pcolor(rlat,play(1:97),bah(1:97,:));
cmap = llsmap5;
set(gca,'ydir','reverse');
set(gca,'yscale','log');
set(gca,'colormap',cmap,'CLim',[-cmax +cmax]);
set(gca,'Ylim',[10 1000]);
shading interp
colorbar
title('ERA TZ from Chris')

%{
clear data dataMap
data    = bah(1:97,:);
dataMap = smoothn(reshape(squeeze(fits.tz_coef(i500,:,2)),64,72),1);
play97 = play(1:97);
save era_T_zonal_trends_donebyCLH.mat rlat play97 data dataMap rlat65 rlon73 
%}

disp('ret to cont'); pause;

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
