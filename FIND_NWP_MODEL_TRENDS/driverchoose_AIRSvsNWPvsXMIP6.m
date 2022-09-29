function [airsChoice,nwpChoice,xmip6Choice] = driverchoose_AIRSvsNWPvsXMIP6(iA,iNWP,iXMIP6,iNumYears)

% this serves to supercede
%   plot_ERA_ERA5_AIRSL3_AMIP6_trends.m
%   plot_ERA_ERA5_AIRSL3_CMIP6_trends.m
% function [airsChoice,nwpChoice,xmip6Choice] = driverchoose_AIRSvsNWPvsXMIP6(iA,iNWP,iXMIP6)
%   iA  = airs L3 choice = +1 Joel Susskind AIRS L3 (default)
%                          -1 Chris Barnet CLIMCAPS
%   iNWP = 1   ERA-I
%          2   MERRA2
%          5   ERA5 (default)
%   iXMIP6 = CMIP -1 default
%            AMIP +1

if nargin == 0
  iA     = 1;
  iNWP   = 5;
  iXMIP6 = -1;
  iNumYears = 19;
elseif nargin == 1
  iNWP   = 5;
  iXMIP6 = -1;
  iNumYears = 19;
elseif nargin == 2
  iXMIP6 = -1;
  iNumYears = 19;
elseif nargin == 3
  iNumYears = 19;
end

iNorD = 1;
iAorOrL = 0;

%%iNumYears = 19;

fprintf(1,'driverchoose_AIRSvsNWPvsXMIP6.m : iNumYears = %2i \n',iNumYears)

if iNumYears ~= 12 | iNumYears ~= 19
  airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL,19);
else
  airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL,iNumYears);
end

nwpChoice   = getdata_NWP(iNWP,iNorD,iAorOrL,iNumYears);

xmip6Choice = getdata_XMIP6(iXMIP6,iNorD,iAorOrL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
see eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m

xyz = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','nwp_spectral_trends_cmip6_era5_airsL3_umbc');
xyz = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5_uncX3.mat','nwp_spectral_trends_cmip6_era5_airsL3_umbc');

xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends
nwpChoice

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5
figure(1); clf; colormap(llsmap5)
figure(2); clf; colormap(llsmap5)
figure(3); clf; colormap(llsmap5)
figure(4); clf; colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%

sum(sum(nwpChoice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp))

figure(1); clf; pcolor(nwpChoice.trend_ptemp); shading interp; colorbar;
  caxis([-1 +1]/10)
figure(2); clf; pcolor(nwpChoice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp); shading interp; colorbar; 
  caxis([-1 +1]/100)

figure(1); clf; pcolor(squeeze(nanmean(reshape(nwpChoice.trend_ptemp,100,72,64),2))); shading interp; colorbar;
  caxis([-1 +1]/10); set(gca,'ydir','reverse')
figure(2); clf; pcolor(squeeze(nanmean(reshape(nwpChoice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp,100,72,64),2))); shading interp; colorbar; 
  caxis([-1 +1]/100); set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%

sum(sum(xmip6Choice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp))

figure(3); clf; pcolor(xmip6Choice.trend_ptemp); shading interp; colorbar;
  caxis([-1 +1]/10)
figure(4); clf; pcolor(xmip6Choice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp); shading interp; colorbar; 
  caxis([-1 +1]/100)

figure(3); clf; pcolor(squeeze(nanmean(reshape(xmip6Choice.trend_ptemp,100,72,64),2))); shading interp; colorbar;
  caxis([-1 +1]/10); set(gca,'ydir','reverse')
figure(4); clf; pcolor(squeeze(nanmean(reshape(xmip6Choice.trend_ptemp - xyz.nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp,100,72,64),2))); shading interp; colorbar; 
  caxis([-1 +1]/100); set(gca,'ydir','reverse')
%}
