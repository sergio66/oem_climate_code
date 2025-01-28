iZonalorAll = +1;
iSave = +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cluster_do_the_fits_airsL3_ratesv7_tiles_radiances.m

iDo = +1;
iDo = -1;
if iDo > 0  
  disp('reading in AIRS L3 fluxes   + for 10, . for 1')
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
      flux64x72ta.anomflux(:,ibah,:,:) = a.thestatsradtrend64x72.anomflux;          %% 72 x 14 x 262 --> 72 x 64 x 14 x 262
      flux64x72ta.trendflux(:,ibah,:)  = a.thestatsradtrend64x72.trendflux;         %% 72 x 14       --> 72 x 64 x 14
      flux64x72ta.trendflux_unc(:,ibah,:) = a.thestatsradtrend64x72.trendflux_unc;  %% 72 x 14       --> 72 x 64 x 14
      flux64x72ta.BTtrend(:,ibah,:)    = a.thestatsradtrend64x72.BTtrend;
      flux64x72ta.BTtrenderr(:,ibah,:) = a.thestatsradtrend64x72.BTtrenderr;
      flux64x72ta.radanom(:,ibah,:,:)  = a.thestatsradtrend64x72.radanom;    %% 72 x 2645 x 262 --> 72 x 64 x 2645 x 262
      flux64x72ta.BTanom(:,ibah,:,:)   = a.thestatsradtrend64x72.BTanom;     %% 72 x 2645 x 262 --> 72 x 64 x 2645 x 262 
    else
      %% much smaller, zonal avg
      flux64x72ta.anomflux(ibah,:,:) = nanmean(a.thestatsradtrend64x72.anomflux,1);            %% 72 x 14 x 262 --> 64 x 14 x 262
      flux64x72ta.trendflux(ibah,:)  = nanmean(a.thestatsradtrend64x72.trendflux,1);           %% 72 x 14       --> 64 x 14
      flux64x72ta.trendflux_unc(:,ibah,:) = nanmean(a.thestatsradtrend64x72.trendflux_unc,1);  %% 72 x 14       --> 72 x 64 x 14
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
    trendflux_unc = flux64x72ta.trendflux_unc; 
    RRTM_bands = a.thestatsradtrend64x72.RRTM_bands;
      if iZonalorAll == +1
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat comment anomflux trendflux RRTMbands
      else
        save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all72x64_anomflux_14RRTMbands.mat comment anomflux trendflux RRTMbands
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

%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m calls
%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m

iDo = -1;
iDo = +1;
if iDo > 0

  disp('reading in MERRA2 fluxes  + for 10, . for 1')
  for ibah = 1 : 64;
    if mod(ibah,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end

    %%% look at era5x.thesave
    %era5x = load('~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat17_2002_09_2024_08.mat');
    merra2x = load('~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_rlat17_2002_09_2022_08.mat');

    %% made by driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
    fileERA5 = ['/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ibah,'%02i') '_2002_09_2024_06.mat'];
    fileERA5 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ibah,'%02i') '_2002_09_2024_08.mat'];

    %% made by driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH.m
    fileMERRA2 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_idRH_5_rlat' num2str(ibah,'%02i') '_2002_09_2022_08.mat'];
    fileMERRA2 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_idRH_5_rlat' num2str(ibah,'%02i') '_2002_09_2024_08.mat'];

    %era5x = load(fileERA5);
    merra2x = load(fileMERRA2);
  
    %% plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
    fluxMERRA2.RRTM_bands = merra2x.thesave.RRTM_bands;
    if iZonalorAll == -1
      fluxMERRA2.stanom(ibah,:,:)        = merra2x.thesave.sxt_anom;      %% 72 x 264 --> 64 x 72 x 264
      fluxMERRA2.bt1231anom(ibah,:,:)    = merra2x.thesave.xbt1231_anom;  %% 72 x 264 --> 64 x 72 x 264
  
      fluxMERRA2.anomflux(ibah,:,:,:)    = merra2x.thesave.xanomflux;      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxMERRA2.trendflux(ibah,:,:)     = merra2x.thesave.xtrendflux;     %% 72 x 14       --> 64 x 72 x 14
      fluxMERRA2.trendflux_unc(ibah,:,:) = merra2x.thesave.xtrendflux_unc; %% 72 x 14       --> 64 x 72 x 14
  
      fluxMERRA2.BTanom(ibah,:,:,:)    = merra2x.thesave.xanomSpectral;        %% 2645 x 72 x 264
      fluxMERRA2.BTtrend(ibah,:,:)     = merra2x.thesave.xtrendSpectral;       %% 2645 x 72  
      fluxMERRA2.BTtrend_unc(ibah,:,:) = merra2x.thesave.xtrendSpectral_unc;   %% 2645 x 72
  
      fluxMERRA2.BTanom_zonal(ibah,:,:)    = merra2x.thesave.xanom;            %% 2645 x 264
      fluxMERRA2.BTtrend_zonal(ibah,:)     = merra2x.thesave.xtrend;           %% 2645 x 1
      fluxMERRA2.BTtrend_zonal_unc(ibah,:) = merra2x.thesave.xtrend_unc;       %% 2645 x 1
    else
      fluxMERRA2.stanom(ibah,:)        = nanmean(merra2x.thesave.sxt_anom,1);      %% 72 x 264 --> 64 x 72 x 264
      fluxMERRA2.bt1231anom(ibah,:)    = nanmean(merra2x.thesave.xbt1231_anom,1);  %% 72 x 264 --> 64 x 72 x 264
  
      fluxMERRA2.anomflux(ibah,:,:)    = nanmean(merra2x.thesave.xanomflux,1);      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxMERRA2.trendflux(ibah,:)     = nanmean(merra2x.thesave.xtrendflux,1);     %% 72 x 14       --> 64 x 72 x 14
      fluxMERRA2.trendflux_unc(ibah,:) = nanmean(merra2x.thesave.xtrendflux_unc,1); %% 72 x 14       --> 64 x 72 x 14
  
      fluxMERRA2.BTanom(ibah,:,:)    = nanmean(merra2x.thesave.xanomSpectral,2);        %% 2645 x 72 x 264        data = squeeze(tcalc(iii,ilon,:));
      fluxMERRA2.BTtrend(ibah,:)     = nanmean(merra2x.thesave.xtrendSpectral,2);       %% 2645 x 72  
      fluxMERRA2.BTtrend_unc(ibah,:) = nanmean(merra2x.thesave.xtrendSpectral_unc,2);   %% 2645 x 72
  
      %fluxMERRA2.BTanom_zonal(ibah,:)    = nanmean(merra2x.thesave.xanom,1);            %% 2645 x 264
      %fluxMERRA2.BTtrend_zonal(ibah)     = nanmean(merra2x.thesave.xtrend,1);           %% 2645 x 1
      %fluxMERRA2.BTtrend_zonal_unc(ibah) = nanmean(merra2x.thesave.xtrend_unc,1);       %% 2645 x 1
%%    fluxMERRA2.BTanom_zonal(ibah,:,:)    = merra2x.thesave.xanom;            %% 2645 x 264        data = tcalcavg(iii,:); ????
      fluxMERRA2.BTtrend_zonal(ibah,:)     = merra2x.thesave.xtrend;           %% 2645 x 1
      fluxMERRA2.BTtrend_zonal_unc(ibah,:) = merra2x.thesave.xtrend_unc;       %% 2645 x 1
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%

  if iSave > 0
    
    comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_flux.m';
    anomflux = fluxMERRA2.anomflux; 
    trendflux = fluxMERRA2.trendflux; 
    trendflux_unc = fluxMERRA2.trendflux_unc; 
    RRTM_bands = merra2x.thesave.RRTM_bands;
      if iZonalorAll == +1
        save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_anomflux_14RRTMbands.mat comment anomflux trendflux*
      else
        save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all72x64_anomflux_14RRTMbands.mat comment anomflux trendflux*
      end
    
    BTtrend = fluxMERRA2.BTtrend;
    BTtrenderr = fluxMERRA2.BTtrend_unc;
      if iZonalorAll == +1
        save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_BTtrend_2645chans.mat comment BTtrend BTtrenderr
      else
        save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all72x64_BTtrend_2645chans.mat comment BTtrend BTtrenderr
      end
    
    %radanom = fluxMERRA2.radanom; 
    %  if iZonalorAll == +1
    %    save -v7.3../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all64.mat comment radanom
    %  else
    %    save -v7.3../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all72x64.mat comment radanom
    %  end
    
    BTanom = fluxMERRA2.BTanom; 
    BTtrend_zonal = fluxMERRA2.BTtrend_zonal;
    BTtrend_zonalerr = fluxMERRA2.BTtrend_zonal_unc;
      if iZonalorAll == +1
        save -v7.3  ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64.mat comment BTanom BTtrend_zonal BTtrend_zonalerr
      else
        save -v7.3  ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all72x64.mat comment BTanom BTtrend_zonal BTtrend_zonalerr
      end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% driver_computeERA5_monthly_trends_desc_or_asc_64latbins_deltaRH.m
era5 = load('ERA5_atm_data_2002_09_to_2024_08_trends_desc_64latbins.mat');

%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m calls
%% AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m

%% ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
  idRH = +1;   %% keep WV   constant
  idRH = +2;   %% keep CO2  constant
  idRH = +3;   %% keep T,ST constant
  idRH = +4;   %% keep RH   constant
  idRH = +5;   %% put in everything, including clouds
idRH = input('Enter idRH (1:5) : ');

iDo = +1;
if iDo > 0

  fprintf(1,'reading in ERA5 fluxes deltaRH = %2i    + for 10, . for 1 \n',idRH)
  for ibah = 1 : 64;
    if mod(ibah,10) == 0
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end

    %%% look at merra2x.thesave
    merra2x = load(['~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_idRH_'  num2str(idRH) '_rlat17_2002_09_2024_08.mat']);

    fileMERRA2 = ['/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ibah,'%02i') '_2002_09_2024_06.mat'];

    if idRH ~= 2
      fileMERRA2 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/reconstruct_merra2_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ibah,'%02i') '_2002_09_2024_08.mat'];
    else
      fileMERRA2 = ['../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2_ConstG/reconstruct_merra2_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ibah,'%02i') '_2002_09_2024_08.mat'];
    end 

    merra2x = load(fileMERRA2);
  
    %% plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
    fluxMERRA2.RRTM_bands = merra2x.thesave.RRTM_bands;
    if iZonalorAll == -1
      fluxMERRA2.stanom(ibah,:,:)        = merra2x.thesave.sxt_anom;      %% 72 x 264 --> 64 x 72 x 264
      fluxMERRA2.bt1231anom(ibah,:,:)    = merra2x.thesave.xbt1231_anom;  %% 72 x 264 --> 64 x 72 x 264
  
      fluxMERRA2.anomflux(ibah,:,:,:)    = merra2x.thesave.xanomflux;      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxMERRA2.trendflux(ibah,:,:)     = merra2x.thesave.xtrendflux;     %% 72 x 14       --> 64 x 72 x 14
      fluxMERRA2.trendflux_unc(ibah,:,:) = merra2x.thesave.xtrendflux_unc; %% 72 x 14       --> 64 x 72 x 14
  
      fluxMERRA2.BTanom(ibah,:,:,:)    = merra2x.thesave.xanomSpectral;        %% 2645 x 72 x 264
      fluxMERRA2.BTtrend(ibah,:,:)     = merra2x.thesave.xtrendSpectral;       %% 2645 x 72  
      fluxMERRA2.BTtrend_unc(ibah,:,:) = merra2x.thesave.xtrendSpectral_unc;   %% 2645 x 72
  
      fluxMERRA2.BTanom_zonal(ibah,:,:)    = merra2x.thesave.xanom;            %% 2645 x 264
      fluxMERRA2.BTtrend_zonal(ibah,:)     = merra2x.thesave.xtrend;           %% 2645 x 1
      fluxMERRA2.BTtrend_zonal_unc(ibah,:) = merra2x.thesave.xtrend_unc;       %% 2645 x 1
    else
      fluxMERRA2.stanom(ibah,:)        = nanmean(merra2x.thesave.sxt_anom,1);      %% 72 x 264 --> 64 x 72 x 264
      fluxMERRA2.bt1231anom(ibah,:)    = nanmean(merra2x.thesave.xbt1231_anom,1);  %% 72 x 264 --> 64 x 72 x 264
  
      fluxMERRA2.anomflux(ibah,:,:)    = nanmean(merra2x.thesave.xanomflux,1);      %% 72 x 14 x 264 --> 64 x 72 x 14 x 264
      fluxMERRA2.trendflux(ibah,:)     = nanmean(merra2x.thesave.xtrendflux,1);     %% 72 x 14       --> 64 x 72 x 14
      fluxMERRA2.trendflux_unc(ibah,:) = nanmean(merra2x.thesave.xtrendflux_unc,1); %% 72 x 14       --> 64 x 72 x 14
  
      fluxMERRA2.BTanom(ibah,:,:)    = nanmean(merra2x.thesave.xanomSpectral,2);        %% 2645 x 72 x 264        data = squeeze(tcalc(iii,ilon,:));
      fluxMERRA2.BTtrend(ibah,:)     = nanmean(merra2x.thesave.xtrendSpectral,2);       %% 2645 x 72  
      fluxMERRA2.BTtrend_unc(ibah,:) = nanmean(merra2x.thesave.xtrendSpectral_unc,2);   %% 2645 x 72
  
      %fluxMERRA2.BTanom_zonal(ibah,:)    = nanmean(merra2x.thesave.xanom,1);            %% 2645 x 264
      %fluxMERRA2.BTtrend_zonal(ibah)     = nanmean(merra2x.thesave.xtrend,1);           %% 2645 x 1
      %fluxMERRA2.BTtrend_zonal_unc(ibah) = nanmean(merra2x.thesave.xtrend_unc,1);       %% 2645 x 1
%%    fluxMERRA2.BTanom_zonal(ibah,:,:)    = merra2x.thesave.xanom;            %% 2645 x 264        data = tcalcavg(iii,:); ????
      fluxMERRA2.BTtrend_zonal(ibah,:)     = merra2x.thesave.xtrend;           %% 2645 x 1
      fluxMERRA2.BTtrend_zonal_unc(ibah,:) = merra2x.thesave.xtrend_unc;       %% 2645 x 1
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%

  if iSave > 0
    
    comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_flux.m';
    anomflux = fluxMERRA2.anomflux; 
    trendflux = fluxMERRA2.trendflux; 
    trendflux_unc = fluxMERRA2.trendflux_unc; 
    RRTM_bands = merra2x.thesave.RRTM_bands;
      if iZonalorAll == +1
        saver = ['save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_idRH_' num2str(idRH) '_anomflux_14RRTMbands.mat comment anomflux trendflux* '];
        eval(saver)
      else
        saver = ['save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all72x64_idRH_' num2str(idRH) '_anomflux_14RRTMbands.mat comment anomflux trendflux* '];
        eval(saver)
      end
    
    BTtrend = fluxMERRA2.BTtrend;
    BTtrenderr = fluxMERRA2.BTtrend_unc;
      if iZonalorAll == +1
        saver = ['save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_idRH_' num2str(idRH) '_BTtrend_2645chans.mat comment BTtrend BTtrenderr '];
        eval(saver)
      else
        saver = ['save ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all72x64_idRH_' num2str(idRH) '_BTtrend_2645chans.mat comment BTtrend BTtrenderr '];
        eval(saver)
      end
    
    %radanom = fluxMERRA2.radanom; 
    %  if iZonalorAll == +1
    %    save -v7.3../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all64.mat comment radanom
    %  else
    %    save -v7.3../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all72x64.mat comment radanom
    %  end
    
    BTanom = fluxMERRA2.BTanom; 
    BTtrend_zonal = fluxMERRA2.BTtrend_zonal;
    BTtrend_zonalerr = fluxMERRA2.BTtrend_zonal_unc;
      if iZonalorAll == +1
        saver = ['save -v7.3  ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64_idRH_' num2str(idRH) '.mat comment BTanom BTtrend_zonal BTtrend_zonalerr '];
        eval(saver)
      else
        saver = ['save -v7.3  ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all72x64_idRH_' num2str(idRH) '.mat comment BTanom BTtrend_zonal BTtrend_zonalerr '];
        eval(saver)
      end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
>> !ls -lth /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_btanom_all64_anomflux_14RRTMbands.mat
-rw-rw-r-- 1 sergio pi_strow 1.3M Sep 14 14:22 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_btanom_all64_anomflux_14RRTMbands.mat
>> !ls -lth //asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_btanom_all64_BTtrend_2645chans.mat
-rw-rw-r-- 1 sergio pi_strow 1.1M Sep 14 14:22 //asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_btanom_all64_BTtrend_2645chans.mat
>> !ls -lth /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all64.mat
-rw-rw-r-- 1 sergio pi_strow 264M Sep 14 14:23 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_radanom_2645chans_all64.mat
>> !ls -lth /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64.mat
-rw-rw-r-- 1 sergio pi_strow 280M Sep 14 14:23 /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64.mat

>> !ls -lth ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_anomflux_14RRTMbands.mat
-rw-rw-r-- 1 sergio pi_strow 1.3M Sep 14 14:15 ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_MERRA2_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_anomflux_14RRTMbands.mat
>> !ls -lth ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_BTtrend_2645chans.mat
-rw-rw-r-- 1 sergio pi_strow 1.2M Sep 14 14:16 ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_MERRA2_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_btanom_all64_BTtrend_2645chans.mat
>> !ls -lth ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64.mat
-rw-rw-r-- 1 sergio pi_strow 253M Sep 14 14:16 ../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_MERRA2_AIRSL3_CMIP6/STS/NIGHTorAVG/MERRA2/merra2_64x72_Sept2002_Aug2022_20yr_desc_BTanom_2645chans_all64.mat
%}
