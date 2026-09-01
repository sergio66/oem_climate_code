%% LAYER 1
%% topts.set_era5_cmip6_airsL3 = 6;           %% use AMIP6    a priori
%% topts.set_era5_cmip6_airsL3 = 3;           %% use AIRSL3   a priori
%% topts.set_era5_cmip6_airsL3 = -3;          %% use CLIMCAPS a priori
%% topts.set_era5_cmip6_airsL3 = 2;           %% use MERRA2   a priori
%% topts.set_era5_cmip6_airsL3 = 5;           %% use ERA5     a priori
%% topts.set_era5_cmip6_airsL3 = 8;           %% use MLS      a priori >>>>
%% topts.set_era5_cmip6_airsL3 = 0;           %% use 0        a priori, DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% then check these ind settings
%%   check_settings.m:7:settings.set_era5_cmip6_airsL3_WV_T_O3 
%%     is using ERA5 or AIRS L3 or MERRA to set rates, you can choose to see -1:WV/T/ST/O3 or +1/+2/+3/+4/+5/+40  for WV/T+ST/O3/T/ST/LowerT only
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == -1   : reset alllayers/all T/WV/O3
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +1   : set WV
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +2   : set T/ST
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +3   : set O3
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +4   : set T
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +5   : set ST
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +10  : set lower WV
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +100 : set lower WV/upper WV with MLS
%%   settings.set_era5_cmip6_airsL3_WV_T_O3 == +40  : set lower T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbINITIAL = xb;

if settings.set_era5_cmip6_airsL3 == 0
  disp('apriori will be using ZERO (no NWP trends being used, nothing in lower at or upper atm)')

elseif settings.set_era5_cmip6_airsL3 == 5
  disp(' apriori will be using ERA5 trends')
  set_apriori_usingERA5

elseif settings.set_era5_cmip6_airsL3 == 6
  disp(' apriori will be using CMIP6 trends')
  set_apriori_usingCMPI6
  
elseif abs(settings.set_era5_cmip6_airsL3) == 3
  set_apriori_usingAIRSL3
  
elseif settings.set_era5_cmip6_airsL3 == 2
  disp(' apriori will be using MERRA2 trends')
  set_apriori_usingMERRA2  

elseif settings.set_era5_cmip6_airsL3 == 8
  disp(' apriori will be using MLS trends')
  set_apriori_usingMLS

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbFINAL = xb;

printarray([(1:length(xbFINAL))' xbINITIAL xbFINAL xbINITIAL-xbFINAL],'xbINITIAL and xbFINAL in set_apriori_ERA5_MERRA2_or_AIRSL3_MLS_geophysical.m')

