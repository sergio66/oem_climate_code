%% check_settings.m:7:settings.set_era5_cmip6_airsL3_WV_T_O3 = -1; %% is using ERA5 or AIRS L3 or MERRA to set rates, you can choose to see -1:WV/T/ST/O3 or +1/+2/+3/+4/+5/+40  for WV/T+ST/O3/T/ST/LowerT only
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == -1  : set all
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +1  : set WV
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +2  : set T/ST
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +3  : set O3
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +4  : set T
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +5  : set ST
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +10 : set lower WV
%% settings.set_era5_cmip6_airsL3_WV_T_O3 == +40 : set lower T

if settings.set_era5_cmip6_airsL3 == 5
  disp(' apriori will be using ERA5 trends')
  vars_cmip6_era5_airsL3_umbc = whos('-file','nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');
  xrates = load('nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat','era5_100_layertrends');
  xrates = xrates.era5_100_layertrends;

  bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  bad = find(isnan(xrates.gas_1)); xrates.gas_1(bad) = 0;
  bad = find(isnan(xrates.gas_3)); xrates.gas_3(bad) = 0;

  ix = driver.iLon;
  iy = driver.iLat;
  iz = (iy-1)*72 + ix;
  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 5
    %% stemp
    boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  else
    boo = 6;
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 10
    %% WV(z)
    wvscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; wvscale = ones(size(wvscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 10
      wvscale(end-3:end) = [1 2 3 4]/4; wvscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* wvscale; 
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 4 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
    %% T(z)
    tscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; tscale = ones(size(tscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
      tscale(end-3:end) = [1 2 3 4]/4; tscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* tscale;
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 3
    %% O3(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  xb = reshape(xb,length(xb),1);

elseif settings.set_era5_cmip6_airsL3 == 6
  disp(' apriori will be using CMIP6 trends')
  vars_cmip6_era5_airsL3_umbc = whos('-file','nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');
  xrates = load('nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat','cmip6_100_layertrends');
  xrates = xrates.cmip6_100_layertrends;

  bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  bad = find(isnan(xrates.gas_1)); xrates.gas_1(bad) = 0;
  bad = find(isnan(xrates.gas_3)); xrates.gas_3(bad) = 0;

  ix = driver.iLon;
  iy = driver.iLat;
  iz = (iy-1)*72 + ix;
  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 5
    %% ST
    boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  else
    boo = 6;
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 1
    %% WV(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 4 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
    %% T(z)
    tscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; tscale = ones(size(tscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
      tscale(end-3:end) = [1 2 3 4]/4; tscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* tscale;
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 3
    %% O3(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end
  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  xb = reshape(xb,length(xb),1);

elseif settings.set_era5_cmip6_airsL3 == 3
  disp(' apriori will be using AIRS L3 trends')
  vars_cmip6_era5_airsL3_umbc = whos('-file','nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');
  xrates = load('nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat','airsL3_100_layertrends');
  xrates = xrates.airsL3_100_layertrends;

  bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  bad = find(isnan(xrates.gas_1)); xrates.gas_1(bad) = 0;
  bad = find(isnan(xrates.gas_3)); xrates.gas_3(bad) = 0;

  ix = driver.iLon;
  iy = driver.iLat;
  iz = (iy-1)*72 + ix;
  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 5
    %% ST
    boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  else
    boo = 6;
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 1
    %% WV(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 4 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
    %% T(z)
    tscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; tscale = ones(size(tscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
      tscale(end-3:end) = [1 2 3 4]/4; tscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* tscale;
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 3
    %% O3(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end
  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  xb = reshape(xb,length(xb),1);

elseif settings.set_era5_cmip6_airsL3 == 2
  disp(' apriori will be using MERRA2 trends')
  zrates = load('../FIND_NWP_MODEL_TRENDS/MERRA2_atm_data_2002_09_to_2021_08_trends.mat');
  xrates.stemp = zrates.trend_stemp;
  xrates.ptemp = zrates.trend_ptemp;
  xrates.gas_1 = zrates.trend_gas_1;
  xrates.gas_3 = zrates.trend_gas_3;

  bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  bad = find(isnan(xrates.gas_1)); xrates.gas_1(bad) = 0;
  bad = find(isnan(xrates.gas_3)); xrates.gas_3(bad) = 0;

  ix = driver.iLon;
  iy = driver.iLat;
  iz = (iy-1)*72 + ix;
  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 5
    %% ST
    boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  else
    boo = 6;
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 1
    %% WV(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 4 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
    %% T(z)
    tscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; tscale = ones(size(tscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
      tscale(end-3:end) = [1 2 3 4]/4; tscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* tscale;
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 3
    %% O3(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  xb = reshape(xb,length(xb),1);

elseif settings.set_era5_cmip6_airsL3 == 8
  disp(' apriori will be using MLS trends')
  zrates = load('../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
  xrates.stemp = zrates.trend_stemp * 0;  
  xrates.ptemp = zrates.trend_ptemp * 0; bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  xrates.gas_1 = zrates.trend_gas_1; xrates.gas_1(65:100,:) = 0; %% xrates.gas_1(60:100,:) = 0; %% about 250 mb to GND 
  xrates.gas_3 = zrates.trend_gas_3 * 0;

  bad = find(isnan(xrates.ptemp)); xrates.ptemp(bad) = 0;
  bad = find(isnan(xrates.gas_1)); xrates.gas_1(bad) = 0;
  bad = find(isnan(xrates.gas_3)); xrates.gas_3(bad) = 0;

  ix = driver.iLon;
  iy = driver.iLat;
  iz = (iy-1)*72 + ix;
  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 5
    %% ST
    boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  else
    boo = 6;
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 1
    %% WV(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 2 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 4 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
    %% T(z)
    tscale = boo(end)+1 : boo(end)+1 + iNlays_retrieve-1; tscale = ones(size(tscale)); 
    if settings.set_era5_cmip6_airsL3_WV_T_O3 == 40
      tscale(end-3:end) = [1 2 3 4]/4; tscale(1:end-4) = 0;
    end
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve) .* tscale;
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  if settings.set_era5_cmip6_airsL3_WV_T_O3 == -1 | settings.set_era5_cmip6_airsL3_WV_T_O3 == 3
    %% O3(z)
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  else
    boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); 
  end

  % boo = 6;                                             xb(boo)     = xrates.stemp(iz);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_1(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve); 
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.ptemp(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  % boo = (boo(end)+1 : boo(end)+1 + iNlays_retrieve-1); xb(boo) = average_over_5(xrates.gas_3(:,iz),floor(100/iNlays_retrieve),iNlays_retrieve);
  xb = reshape(xb,length(xb),1);
end
