addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /umbc/xfs2/strow/asl/matlib/maps
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS

load llsmap5
cmap = llsmap5;

load umbc_trends_iQuantile50.mat  %% from ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/find_T_RH_trends.m 

aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xumbc.trend_ST,72,64)') ,1), [-90 +90],[-180 +180]); title('UMBC dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_T',1)); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T UMBC'])
colormap(cmap); ylim([1 1000])

figure(3); 
pcolor(xumbc.rlat,xumbc.pavgLAY,xumbc.trend_RH'); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_RH',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH UMBC'])
colormap(cmap); ylim([100 1000])

figure(4); 
pcolor(xumbc.rlat,xumbc.pavgLAY,xumbc.trend_WV'); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_WV',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([100 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac UMBC'])
colormap(cmap); ylim([100 1000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load giss_airsV7_climcps_xmip_era5_merra2.mat %% from see/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/compare_12years_vs_19years_hiatus_warming.m

figure(5)
aslmap(5,xumbc.rlat65,xumbc.rlon73,smoothn(xgiss.giss_trend4608' ,1), [-90 +90],[-180 +180]); title('GISS dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xamip6_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('AMIP6 dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xamip6_12.trend_ptemp(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T AMIP6'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xamip6_12.trend_RH(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH AMIP6'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xamip6_12.trend_gas_1(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac AMIP6'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xcmip6_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('CMIP6 dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xcmip6_12.trend_ptemp(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T CMIP6'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xcmip6_12.trend_RH(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH CMIP6'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xcmip6_12.trend_gas_1(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac CMIP6'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xmerra2_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('MERRA2 dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xmerra2_12.trend_ptemp(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T MERRA2'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xmerra2_12.trend_RH(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH MERRA2'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xmerra2_12.trend_gas_1(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac MERRA2'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xera5_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xera5_12.trend_ptemp(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T ERA5'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xera5_12.trend_RH(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH ERA5'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xera5_12.trend_gas_1(1:97,:); junk = reshape(junk,97,72,64); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac ERA5'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xairsV7_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('AIRSV7 dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xairsV7_12.trend_ptemp; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xairsV7_12.Tlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T AIRSV7'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xairsV7_12.trend_RH; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xairsV7_12.Qlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH AIRSV7'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xairsV7_12.trend_gas_1; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xairsV7_12.Qlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac AIRSV7'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue'); pause(0.1)
aslmap(1,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xclimcaps_12.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('CLIMCAPS dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(2); 
junk = xclimcaps_12.trend_ptemp; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xclimcaps_12.Tlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([1 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T CLIMCAPS'])
colormap(cmap); ylim([1 1000])

figure(3); 
junk = xclimcaps_12.trend_RH; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xclimcaps_12.Qlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH CLIMCAPS'])
colormap(cmap); ylim([100 1000])

figure(4); 
junk = xclimcaps_12.trend_gas_1; junk = permute(junk,[3 1 2]); junk = squeeze(nanmean(junk,2));
pcolor(xumbc.rlat,xclimcaps_12.Qlevs,smoothn(junk,1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac CLIMCAPS'])
colormap(cmap); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z11 = xumbc.trend_ST;
z12 = xairsV7_12.trend_stemp;
z13 = xclimcaps_12.trend_stemp;
z21 = xera5_12.trend_stemp;
z22 = xmerra2_12.trend_stemp;
z23 = xgiss.giss_trend4608;
xx = xamip6_12.trend_stemp;
plotoptions.cx = [-1 +1]*0.15;
plotoptions.yReverseDir = -1;
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'UMBC';
plotoptions.str12 = 'AIRSv7';
plotoptions.str13 = 'CLIMCAPS';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.str23 = 'GISS';
plotoptions.xstr = ' ';
plotoptions.ystr = ' ';
plotoptions.maintitle = 'dST/dt K/yr';
figure(6); clf; aslmap_2x3tiledlayout(z11,z12,z13,z21,z22,z23,6,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%

z11 = xumbc.trend_T;
z12 = xairsV7_12.trend_ptemp;    z12 = squeeze(nanmean(z12,1)); 
  clear junk; for ii=1:64; junk(ii,:) = interp1(log10(xairsV7_12.Tlevs),z12(ii,:),log10(xumbc.pavgLAY),[],'extrap'); end; z12 = junk;
z13 = xclimcaps_12.trend_ptemp;  z13 = squeeze(nanmean(z13,1)); z13 = z13(:,1:97);
z21 = xera5_12.trend_ptemp(1:97,:); z21 = reshape(z21,97,72,64); z21 = squeeze(nanmean(z21,2))';
z22 = xmerra2_12.trend_ptemp(1:97,:); z22 = reshape(z22,97,72,64); z22 = squeeze(nanmean(z22,2))';
z23 = xamip6_12.trend_ptemp; z23 = z23(1:97,:); z23 = reshape(z23,97,72,64); z23 = squeeze(nanmean(z23,2))';

plotoptions.cx = [-1 +1]*0.15;
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'UMBC';
plotoptions.str12 = 'AIRSv7';
plotoptions.str13 = 'CLIMCAPS';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.str23 = 'AMIP6';
plotoptions.ylimits = [+1 1000];
plotoptions.yReverseDir = +1;
plotoptions.yLinearOrLog = -1;
plotoptions.xstr = 'Latitude';
plotoptions.ystr = 'Pressure (mb)';
plotoptions.maintitle = 'dT/dt K/yr';
x = xumbc.rlat;
y = xumbc.pavgLAY;
figure(7); clf; tiled_2x3layout(z11,z12,z13,z21,z22,z23,7,plotoptions,x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%

z11 = xumbc.trend_RH;
z12 = xairsV7_12.trend_RH;    z12 = squeeze(nanmean(z12,1)); 
  clear junk; for ii=1:64; junk(ii,:) = interp1(log10(xairsV7_12.Qlevs),z12(ii,:),log10(xumbc.pavgLAY),[],'extrap'); end; z12 = junk; clear junk
z13 = xclimcaps_12.trend_RH;  z13 = squeeze(nanmean(z13,1)); clear junk; junk = z13; z13 = zeros(size(z12)); z13(:,97-66+1:97) = junk;
z21 = xera5_12.trend_RH(1:97,:); z21 = reshape(z21,97,72,64); z21 = squeeze(nanmean(z21,2))';
z22 = xmerra2_12.trend_RH(1:97,:); z22 = reshape(z22,97,72,64); z22 = squeeze(nanmean(z22,2))';
z23 = xamip6_12.trend_RH; z23 = z23(1:97,:); z23 = reshape(z23,97,72,64); z23 = squeeze(nanmean(z23,2))';

plotoptions.cx = [-1 +1]*0.5;
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'UMBC';
plotoptions.str12 = 'AIRSv7';
plotoptions.str13 = 'CLIMCAPS';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.str23 = 'AMIP6';
plotoptions.ylimits = [+100 1000];
plotoptions.ylimits = [+1 1000];
plotoptions.yReverseDir = +1;
plotoptions.yLinearOrLog = +1;
plotoptions.xstr = 'Latitude';
plotoptions.ystr = 'Pressure (mb)';
plotoptions.maintitle = 'dRH/dt pc/yr';
x = xumbc.rlat;
y = xumbc.pavgLAY;
figure(8); clf; tiled_2x3layout(z11,z12,z13,z21,z22,z23,8,plotoptions,x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%

z11 = xumbc.trend_WV;
z12 = xairsV7_12.trend_gas_1;    z12 = squeeze(nanmean(z12,1)); 
  clear junk; for ii=1:64; junk(ii,:) = interp1(log10(xairsV7_12.Qlevs),z12(ii,:),log10(xumbc.pavgLAY),[],'extrap'); end; z12 = junk; clear junk
z13 = xclimcaps_12.trend_gas_1;  z13 = squeeze(nanmean(z13,1)); clear junk; junk = z13; z13 = zeros(size(z12)); z13(:,97-66+1:97) = junk;
z21 = xera5_12.trend_gas_1(1:97,:); z21 = reshape(z21,97,72,64); z21 = squeeze(nanmean(z21,2))';
z22 = xmerra2_12.trend_gas_1(1:97,:); z22 = reshape(z22,97,72,64); z22 = squeeze(nanmean(z22,2))';
z23 = xamip6_12.trend_gas_1; z23 = z23(1:97,:); z23 = reshape(z23,97,72,64); z23 = squeeze(nanmean(z23,2))';

plotoptions.cx = [-1 +1]*0.02;
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'UMBC';
plotoptions.str12 = 'AIRSv7';
plotoptions.str13 = 'CLIMCAPS';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.str23 = 'AMIP6';
plotoptions.ylimits = [+100 1000];
plotoptions.ylimits = [+1 1000];
plotoptions.yReverseDir = +1;
plotoptions.yLinearOrLog = +1;
plotoptions.xstr = 'Latitude';
plotoptions.ystr = 'Pressure (mb)';
plotoptions.maintitle = 'dWVfrac/dt /yr';
x = xumbc.rlat;
y = xumbc.pavgLAY;
figure(9); clf; tiled_2x3layout(z11,z12,z13,z21,z22,z23,9,plotoptions,x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%

z11 = xumbc.trend_ST;
z12 = xairsV7_12.trend_stemp;
z13 = xclimcaps_12.trend_stemp;
z21 = xera5_12.trend_stemp;
z22 = xmerra2_12.trend_stemp;
%z23 = xgiss.giss_trend4608;
z23 = xamip6_12.trend_stemp;
plotoptions.cx = [-1 +1]*0.15;
plotoptions.yReverseDir = -1;
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'UMBC';
plotoptions.str12 = 'AIRSv7';
plotoptions.str13 = 'CLIMCAPS';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
%plotoptions.str23 = 'GISS';
plotoptions.str23 = 'AMIP6';
plotoptions.xstr = ' ';
plotoptions.ystr = ' ';
plotoptions.maintitle = 'dST/dt K/yr';
figure(10); clf; aslmap_2x3tiledlayout(z11,z12,z13,z21,z22,z23,10,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
addpath /asl/matlib/plotutils
figure(6); aslprint('stemp_trends_2002_2014.pdf');
figure(7); aslprint('tz_trends_2002_2014.pdf');
figure(8); aslprint('RHz_trends_2002_2014.pdf');
figure(9); aslprint('wvfrac_trends_2002_2014.pdf');
figure(10); aslprint('stemp_trends_2002_2014_amip6.pdf');
%}
