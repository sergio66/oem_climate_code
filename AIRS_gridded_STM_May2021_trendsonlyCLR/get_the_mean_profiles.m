%% now need to get in mean profiles

%% see eg SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_gather_RH_rates_AIRSL3_NWP_XMIP.m
%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/plot_driver_gather_gridded_retrieval_results.m
%% see eg FIND_NWP_MODEL_TRENDS/driver_computeERA_16day_trends_desc_or_asc.m

%[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');
[hMean17years,ha,pMean17years,pa]     = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
iLoad = 1;
  iDorA = 1;
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  end
  era_prof = load(fin);
  hTimeStep1 = era_prof.hnew_op;
  pTimeStep1 = era_prof.pnew_op;

%% now see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/find_T_RH_trends.m
h = hTimeStep1; p = pTimeStep1;      %% I been using this in eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_gather_gridded_retrieval_results
h = hMean17years; p = pMean17years;  %% I think I should use this
