%% cluster_do_the_fits_airsL3_ratesv7_tiles_radiances.m
iZonalorAll = +1;
iSave = +1;

iDo = -1;
if iDo > 0  
  disp('reading in AIRS L3 fluxes')
  for ibah = 1 : 64;
    if mod(ibah,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end
    filein = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_latbin_' num2str(ibah,'%02i') '.mat'];
    a = load(filein);
    if iZonalorAll == -1
      %% huge file, eg 72x64x2645x264
      flux64x72ta.anomflux(:,ibah,:,:) = a.thestatsradtrend64x72.anomflux;   %% 72 x 14 x 262 --> 72 x 64 x 14 x 262
      flux64x72ta.trendflux(:,ibah,:)  = a.thestatsradtrend64x72.trendflux;  %% 72 x 14       --> 72 x 64 x 14
      flux64x72ta.BTtrend(:,ibah,:)    = a.thestatsradtrend64x72.BTtrend;
      flux64x72ta.BTtrenderr(:,ibah,:) = a.thestatsradtrend64x72.BTtrenderr;
      flux64x72ta.radanom(:,ibah,:,:)  = a.thestatsradtrend64x72.radanom;    %% 72 x 2645 x 262 --> 72 x 64 x 2645 x 262
      flux64x72ta.BTanom(:,ibah,:,:)   = a.thestatsradtrend64x72.BTanom;     %% 72 x 2645 x 262 --> 72 x 64 x 2645 x 262 
    else
      %% much smaller, zonal avg
      flux64x72ta.anomflux(ibah,:,:) = nanmean(a.thestatsradtrend64x72.anomflux,1);   %% 72 x 14 x 262 --> 64 x 14 x 262
      flux64x72ta.trendflux(ibah,:)  = nanmean(a.thestatsradtrend64x72.trendflux,1);  %% 72 x 14       --> 64 x 14
      flux64x72ta.BTtrend(ibah,:)    = nanmean(a.thestatsradtrend64x72.BTtrend,1);
      flux64x72ta.BTtrenderr(ibah,:) = nanmean(a.thestatsradtrend64x72.BTtrenderr,1);
      flux64x72ta.radanom(ibah,:,:)  = nanmean(a.thestatsradtrend64x72.radanom,1);    %% 72 x 2645 x 262 --> 64 x 2645 x 262
      flux64x72ta.BTanom(ibah,:,:)   = nanmean(a.thestatsradtrend64x72.BTanom,1);     %% 72 x 2645 x 262 --> 64 x 2645 x 262 
    end
  end
  fprintf(1,'\n');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  if iSave > 0
    
    comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_flux.m';
    anomflux = flux64x72ta.anomflux; 
    trendflux = flux64x72ta.trendflux; 
    RRTM_bands = a.thestatsradtrend64x72.RRTM_bands;
      if iZonalorAll == +1
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat comment anomflux trendflux
      else
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all72x64_anomflux_14RRTMbands.mat comment anomflux trendflux
      end
    
    BTtrend = flux64x72ta.BTtrend;
    BTtrenderr = flux64x72ta.BTtrenderr;
      if iZonalorAll == +1
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all64_BTtrend_2645chans.mat comment BTtrend BTtrenderr
      else
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all72x64_BTtrend_2645chans.mat comment BTtrend BTtrenderr
      end
    
    radanom = flux64x72ta.radanom; 
      if iZonalorAll == +1
        save -v7.3 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_radanom_2645chans_all64.mat comment radanom
      else
        save -v7.3 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_radanom_2645chans_all72x64.mat comment radanom
      end
    
    BTanom = flux64x72ta.BTanom; 
      if iZonalorAll == +1
        save -v7.3 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_BTanom_2645chans_all64.mat comment BTanom
      else
        save -v7.3 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_BTanom_2645chans_all72x64.mat comment BTanom
      end
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% driver_computeERA5_monthly_trends_desc_or_asc_64latbins.m
era5 = load('ERA5_atm_data_2002_09_to_2024_08_trends_desc_64latbins.mat');

%% driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m calls
%% plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m

%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_spectral_trends_latbin_1_64_sarta.m

iDo = +1;
if iDo > 0

  disp('reading in ERA5 fluxes')
  for ibah = 1 : 64;
    if mod(ibah,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end

    %%% look at era5x.thesave
    era5x = load('~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat17_2002_09_2024_08.mat');

    fileERA5 = ['/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ibah,'%02i') '_2002_09_2024_06.mat'];
    fileERA5 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ibah,'%02i') '_2002_09_2024_08.mat'];
  
    era5x = load(fileERA5);
  
    %% plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
    fluxERA5.RRTM_bands = era5x.thesave.RRTM_bands;
    if iZonalorAll == -1
      fluxERA5.stanom(ibah,:,:)        = era5x.thesave.sxt_anom;      %% 72 x 264 --> 64 x 72 x 264
      fluxERA5.bt1231anom(ibah,:,:)    = era5x.thesave.xbt1231_anom;  %% 72 x 264 --> 64 x 72 x 264
  
      fluxERA5.anomflux(ibah,:,:,:)    = era5x.thesave.xanomflux;      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxERA5.trendflux(ibah,:,:)     = era5x.thesave.xtrendflux;     %% 72 x 14       --> 64 x 72 x 14
      fluxERA5.trendflux_unc(ibah,:,:) = era5x.thesave.xtrendflux_unc; %% 72 x 14       --> 64 x 72 x 14
  
      fluxERA5.BTanom(ibah,:,:,:)    = era5x.thesave.xanomSpectral;        %% 2645 x 72 x 264
      fluxERA5.BTtrend(ibah,:,:)     = era5x.thesave.xtrendSpectral;       %% 2645 x 72  
      fluxERA5.BTtrend_unc(ibah,:,:) = era5x.thesave.xtrendSpectral_unc;   %% 2645 x 72
  
      fluxERA5.BTanom_zonal(ibah,:,:)    = era5x.thesave.xanom;            %% 2645 x 264
      fluxERA5.BTtrend_zonal(ibah,:)     = era5x.thesave.xtrend;           %% 2645 x 1
      fluxERA5.BTtrend_zonal_unc(ibah,:) = era5x.thesave.xtrend_unc;       %% 2645 x 1
    else
      fluxERA5.stanom(ibah,:)        = nanmean(era5x.thesave.sxt_anom,1);      %% 72 x 264 --> 64 x 72 x 264
      fluxERA5.bt1231anom(ibah,:)    = nanmean(era5x.thesave.xbt1231_anom,1);  %% 72 x 264 --> 64 x 72 x 264
  
      fluxERA5.anomflux(ibah,:,:)    = nanmean(era5x.thesave.xanomflux,1);      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxERA5.trendflux(ibah,:)     = nanmean(era5x.thesave.xtrendflux,1);     %% 72 x 14       --> 64 x 72 x 14
      fluxERA5.trendflux_unc(ibah,:) = nanmean(era5x.thesave.xtrendflux_unc,1); %% 72 x 14       --> 64 x 72 x 14
  
      fluxERA5.BTanom(ibah,:,:)    = nanmean(era5x.thesave.xanomSpectral,2);        %% 2645 x 72 x 264
      fluxERA5.BTtrend(ibah,:)     = nanmean(era5x.thesave.xtrendSpectral,2);       %% 2645 x 72  
      fluxERA5.BTtrend_unc(ibah,:) = nanmean(era5x.thesave.xtrendSpectral_unc,2);   %% 2645 x 72
  
      %fluxERA5.BTanom_zonal(ibah,:)    = nanmean(era5x.thesave.xanom,1);            %% 2645 x 264
      %fluxERA5.BTtrend_zonal(ibah)     = nanmean(era5x.thesave.xtrend,1);           %% 2645 x 1
      %fluxERA5.BTtrend_zonal_unc(ibah) = nanmean(era5x.thesave.xtrend_unc,1);       %% 2645 x 1
      fluxERA5.BTanom_zonal(ibah,:,:)    = era5x.thesave.xanom;            %% 2645 x 264
      fluxERA5.BTtrend_zonal(ibah,:)     = era5x.thesave.xtrend;           %% 2645 x 1
      fluxERA5.BTtrend_zonal_unc(ibah,:) = era5x.thesave.xtrend_unc;       %% 2645 x 1
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%


end
