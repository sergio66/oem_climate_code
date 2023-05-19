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

iNorD = 1;
iAorOrL = 0;

airsL3  = getdata_AIRSL3vsCLIMCAPSL3(+1,iNorD,iAorOrL);
climsL3 = getdata_AIRSL3vsCLIMCAPSL3(+1,iNorD,iAorOrL);

merra2  = getdata_NWP(2,iNorD,iAorOrL);
era5    = getdata_NWP(5,iNorD,iAorOrL);

cmip = getdata_XMIP6(-1,iNorD,iAorOrL);
umbc = getdata_XMIP6(+1,iNorD,iAorOrL);

iFig = 0;
plotoptions.str11 = 'UMBC';    plotoptions.str12 = 'AMIP6';       
plotoptions.str21 = 'ERA5';    plotoptions.str22 = 'MERRA2';   
plotoptions.str31 = 'AIRS L3'; plotoptions.str32 = 'CLIMCAPS';

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end
plotoptions.plotcolors = llsmap5;
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir = +1;
plotoptions.xstr = 'Latitude (deg)';
plotoptions.ystr = 'Pressure (mb)';

x = merra2.trend_rlat64;
y = merra2.trend_plays;

iFig = iFig + 1;
z11 = squeeze(nanmean(reshape(merra2.trend_ptemp,100,72,64),2));
z12 = z11;
z21 = z11;
z22 = z11;
z31 = z11;
z32 = z11;
plotoptions.maintitle = 'dT(z,lat)/dt'; plotoptions.cx = [-1 +1]*0.15; plotoptions.yLimits = [10 1000]; plotoptions.xLimits = [-85 +85];
profile_plots_3x2tiledlayout(x,y,z11,z12,z21,z22,z31,z32,iFig,plotoptions);
