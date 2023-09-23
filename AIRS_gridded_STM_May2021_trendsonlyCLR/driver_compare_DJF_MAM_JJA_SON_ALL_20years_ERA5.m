%% cp /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m compare_05_10_15_20_spectral_rates.m

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps

load llsmap5
clear plotoptions

iAbsOrPercent = +1;
if iAbsOrPercent == 1
  strdiff = ' absolute';
else
  strdiff = ' percent';
end

if ~exist('aALL')
  dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc'

  aDJF = load([dir0 '_DJF.mat']);
  aMAM = load([dir0 '_MAM.mat']);
  aJJA = load([dir0 '_JJA.mat']);
  aSON = load([dir0 '_SON.mat']);
  aALL = load([dir0 '.mat']);
  
  loader = ['load h2645structure.mat'];
  eval(loader)

  [h,ha,p,pa] = rtpread(['summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp']);
end
  
rlat = p.rlat;
rlat = cos(rlat*pi/180);

load latB64.mat
rlon73 = -180 : 5 : +180;  rlat65 = latB2;
rlon = 0.5*(rlon73(1:end-1)+rlon73(2:end));
rlat = 0.5*(rlat65(1:end-1)+rlat65(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

plays = p.plays(:,1000);
pavgLAY = p.plays;
iFig = 6;
rlat100x64 = ones(100,1)*rlat';
plays100x64 = plays(1:100)*ones(1,64);

plotoptions.str11 = 'DJF';
plotoptions.str12 = 'MAM';
plotoptions.str21 = 'JJA';
plotoptions.str22 = 'SON';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = 0;
iFig = iFig + 1; mooDJF = squeeze(aDJF.trend_stemp); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooDJF,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('DJF BT1231 trend')
iFig = iFig + 1; mooMAM = squeeze(aMAM.trend_stemp); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooMAM,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('MAM BT1231 trend')
iFig = iFig + 1; mooJJA = squeeze(aJJA.trend_stemp); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooJJA,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('JJA BT1231 trend')
iFig = iFig + 1; mooSON = squeeze(aSON.trend_stemp); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooSON,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('SON BT1231 trend')
iFig = iFig + 1; mooALL = squeeze(aALL.trend_stemp); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooALL,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('ALL BT1231 trend')
aXYZ.trend_stemp = (aDJF.trend_stemp + aMAM.trend_stemp + aJJA.trend_stemp + aSON.trend_stemp)/4;
iFig = iFig + 1; mooALL = squeeze(aXYZ.trend_stemp);      aslmap(iFig,rlat65,rlon73,smoothn((reshape(aXYZ.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('<SEASON> BT1231 trend')
iFig = iFig;  junk = aALL.trend_stemp - aXYZ.trend_stemp;
cx = [-1 +1]*0.15/10;
if iAbsOrPercent == -1
  junk = 100./mooALL;
  cx = [-1 +1]*100;
end
junk = junk';
aslmap(iFig,rlat65,rlon73,smoothn((reshape(junk,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis(cx); title(['(ALL-SUM) dSKT/dt : ' strdiff])

iFig = iFig + 1; 
plotoptions.yLinearOrLog = +1;
plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
aslmap_2x2tiledlayout(mooDJF,mooMAM,mooJJA,mooSON,iFig,plotoptions);

iFig = iFig+1;
aslmap_polar_2x2tiledlayout(mooDJF,mooMAM,mooJJA,mooSON,iFig,plotoptions);

disp('ret to continue'); pause

%iFig = iFig + 1; 
%  mooALL = squeeze(aXYZ.trend_stemp)./squeeze(aALL.trend_stemp); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL-1,5); colormap(usa2); caxis([-1 +1]*0.5); title('1- ALL<4>/ALL BT1231 trend')
%  mooALL = squeeze(aXYZ.trend_stemp) - squeeze(aALL.trend_stemp0); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL,5); colormap(usa2); caxis([-1 +1]*0.015); title('1- ALL<4>/ALL BT1231 trend')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.trend_ptemp; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('DJF dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.trend_ptemp; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('MAM dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.trend_ptemp; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('JJA dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.trend_ptemp; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('SON dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.trend_ptemp; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('ALL dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.trend_ptemp + aMAM.trend_ptemp + aJJA.trend_ptemp + aSON.trend_ptemp); moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('<SEASONAL> dT/dt');

moo = aALL.trend_ptemp - 1/4*(aDJF.trend_ptemp + aMAM.trend_ptemp + aJJA.trend_ptemp + aSON.trend_ptemp); 
cx = [-1 +1]*0.15/10;
if iAbsOrPercent == -1
  moo = 100./aALL.trend_ptemp;
  cx = [-1 +1]*100;
end
moo = squeeze(nanmean(reshape(moo,100,72,64),2));
pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis(cx); shading interp; ylim([10 1000]); colormap(usa2); colorbar; 
title(['ALL-<SEASONAL> dT/dt : ' strdiff]);

iFig = iFig + 1; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [10 1000];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dT/dt trends [K/yr]';
plotoptions.cx = [-1 +1]*0.151;
miaow11 = squeeze(nanmean(reshape(aDJF.trend_ptemp,100,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.trend_ptemp,100,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.trend_ptemp,100,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.trend_ptemp,100,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.trend_gas_1; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('DJF d(trend_gas_1)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.trend_gas_1; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('MAM d(trend_gas_1)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.trend_gas_1; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('JJA d(trend_gas_1)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.trend_gas_1; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('SON d(trend_gas_1)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.trend_gas_1; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('ALL d(trend_gas_1)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.trend_gas_1 + aMAM.trend_gas_1 + aJJA.trend_gas_1 + aSON.trend_gas_1); moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('<SEASONAL> d(trend_gas_1)/dt');

moo = aALL.trend_gas_1 - 1/4*(aDJF.trend_gas_1 + aMAM.trend_gas_1 + aJJA.trend_gas_1 + aSON.trend_gas_1); 
cx = [-1 +1]*0.015/10;
if iAbsOrPercent == -1
  moo = 100./aALL.trend_gas_1;
  cx = [-1 +1]*100;
end
moo = squeeze(nanmean(reshape(moo,100,72,64),2));
pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis(cx); shading interp; ylim([10 1000]); colormap(usa2); colorbar; 
title(['ALL-<SEASONAL> dtrend_gas_1/dt : ' strdiff]);

iFig = iFig + 1; 
plotoptions.yLinearOrLog = +1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [100 1000];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dfracWV/dt trends [1/yr]';
plotoptions.cx = [-1 +1]*0.151/10;
miaow11 = squeeze(nanmean(reshape(aDJF.trend_gas_1,100,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.trend_gas_1,100,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.trend_gas_1,100,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.trend_gas_1,100,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.trend_gas_3; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('DJF d(trend_gas_3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.trend_gas_3; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('MAM d(trend_gas_3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.trend_gas_3; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('JJA d(trend_gas_3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.trend_gas_3; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('SON d(trend_gas_3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.trend_gas_3; moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('ALL d(trend_gas_3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.trend_gas_3 + aMAM.trend_gas_3 + aJJA.trend_gas_3 + aSON.trend_gas_3); moo = squeeze(nanmean(reshape(moo,100,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('<SEASONAL> d(trend_gas_3)/dt');

moo = aALL.trend_gas_3 - 1/4*(aDJF.trend_gas_3 + aMAM.trend_gas_3 + aJJA.trend_gas_3 + aSON.trend_gas_3); 
cx = [-1 +1]*0.015/10;
if iAbsOrPercent == -1
  moo = 100./aALL.trend_gas_3;
  cx = [-1 +1]*100;
end
moo = squeeze(nanmean(reshape(moo,100,72,64),2));
pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis(cx); shading interp; ylim([10 1000]); colormap(usa2); colorbar; 
title(['ALL-<SEASONAL> dtrend_gas_3/dt : ' strdiff]);

iFig = iFig + 1; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [0.010 100];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dfracO3/dt trends [1/yr]';
plotoptions.cx = [-1 +1]*0.151/10;
miaow11 = squeeze(nanmean(reshape(aDJF.trend_gas_3,100,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.trend_gas_3,100,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.trend_gas_3,100,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.trend_gas_3,100,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);


disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = -1;
if iPrint > 0
  figure(07); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_skt_trends_lat_lon_4panels20_years_seasonal.pdf');
  figure(08); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_skt_trends_lat_lon_4panels20_years_seasonal_polar.pdf');
  figure(15); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_tz_trends_zonal_p_4panels20_years_seasonal.pdf')
  figure(22); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_gas_1_trends_zonal_p_4panels20_years_seasonal.pdf');
  figure(29); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_gas_3_trends_zonal_p_4panels20_years_seasonal.pdf');
end
