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

%%%%%%%%%%%%%%%%%%%%%%%%%
