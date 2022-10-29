addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/re_plot_zonal_trends_umbc_airsL3_models.m
disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');

ixAorOorL = 0;  %% +1 = ocean, -1 = Land, 0 = both

ixAorOorL = input('Enter region : ');
if length(ixAorOorL) == 0 
  ixAorOorL = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iL3orCLIMCAPS')
  error('iL3orCLIMCAPS DNE')
end

if iL3orCLIMCAPS == -1
  Airs_Lat0 = Airs_Lat;
  Airs_Lon0 = Airs_Lon;
  [Airs_Lon,Airs_Lat] = meshgrid(Airs_Lon0,Airs_Lat0);
end

iDorA = input('Enter Asc(-1) or Desc (+1) : ');
if length(iDorA) == 0
  iDorA = +1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = input('do NATIVE 180 bins? (-1/+1) (-1 = default) : ');
if length(iDo) == 0
  iDo = -1;
end

if iDo > 0
  woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan
  
  iXStart = 1; iXEnd = 180;

  disp('NATIVE 180 ZONES')
  clear latbins* save_lat*
  latbins = -90:+90;
  save_lat = 0.5*(latbins(1:end-1)+latbins(2:end));
  
  disp(' ')
  disp('turning AIRS L3 data from gridded Lat/Lon bins to Zonal Lat bins');
  plot_L3_data_native               %% copied from /home/sergio/MATLABCODE/AIRS_L3/plot_L3_data.m; need to modify this so I save it and can reuse it for eg shorter time spans!!!!
  
  disp(' ')
  disp('doing the geophysical rates (do_profilerate_fit includes lag1) and anoms')
  do_the_fits_airsL3_ratesv7_native %% copied from /home/sergio/MATLABCODE/AIRS_L3/do_the_fits_airsL3_ratesv7
  
  %disp(' ')
  %disp('doing the radiances, including lag1')
  %do_the_rads_airsL3
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = input('do EQ AREA 40 latbins? (-1/+1) (-1 = default) : ');
if length(iDo) == 0
  iDo = -1;
end

if iDo > 0
  woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan
  
  iXStart = 1; iXEnd = 40;

  disp('40 EQUAL AREA ZONAL RESULTS')
  clear latbins* save_lat*
  latbins = equal_area_spherical_bands(20);    
  save_lat = 0.5*(latbins(1:end-1)+latbins(2:end));
  
  disp(' ')
  disp('turning AIRS L3 data from gridded Lat/Lon bins to Zonal Lat bins');
  plot_L3_data_zonal               %% copied from /home/sergio/MATLABCODE/AIRS_L3/plot_L3_data.m

  disp(' ')
  disp('doing the geophysical rates (do_profilerate_fit includes lag1) and anoms')
  do_the_fits_airsL3_ratesv7_zonal %% copied from /home/sergio/MATLABCODE/AIRS_L3/do_the_fits_airsL3_ratesv7
  
  %disp(' ')
  %disp('doing the radiances, including lag1')
  %do_the_rads_airsL3
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = input('do 64x72 TILES? (-1/+1) (-1 = default) : ');
if length(iDo) == 0
  iDo = -1;
end

if iDo > 0
  woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan
  
  disp('TILE RESULTS')
  clear latbins* save_lat*
  clear save64*
  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  
  drlon = 5; 
  rlon = -180 : drlon : +180;      %% 73 edges, so 72 bins
  rlat = latB2;                    %% 65 edges, so 64 bins
  save_lat64x72 = 0.5*(rlat(1:end-1)+rlat(2:end));
  save_lon64x72 = 0.5*(rlon(1:end-1)+rlon(2:end));
  
  disp(' ')
  disp('turning AIRS L3 data from gridded Lat/Lon bins to Tiles bins VERY VERY SLOOOOOOOOWWWWWWWWW');
  plot_L3_data_tiles
  
  disp(' ')
  disp('doing the geophysical rates (do_profilerate_fit includes lag1) and anoms')
  iSlowORFast = +1;
  do_the_fits_airsL3_ratesv7_tiles
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = input('do 64x72 TILES fastgribd? (-1/+1) (-1 = default) : ');
if length(iDo) == 0
  iDo = -1;
end

if iDo > 0
  woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan
  
  disp('TILE RESULTS')
  clear latbins* save_lat*
  clear save64*
  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  
  drlon = 5; 
  rlon = -180 : drlon : +180;      %% 73 edges, so 72 bins
  rlat = latB2;                    %% 65 edges, so 64 bins
  save_lat64x72 = 0.5*(rlat(1:end-1)+rlat(2:end));
  save_lon64x72 = 0.5*(rlon(1:end-1)+rlon(2:end));
  
  disp(' ')
  disp('turning AIRS L3 data from gridded Lat/Lon bins to Tiles bins VERY VERY SLOOOOOOOOWWWWWWWWW');
  plot_L3_data_tiles_fast_grib
  
  disp(' ')
  disp('doing the geophysical rates (do_profilerate_fit includes lag1) and anoms')
  iSlowORFast = -1;
  do_the_fits_airsL3_ratesv7_tiles
end
