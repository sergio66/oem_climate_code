function driver = get_rates(driver,settings,iNoiseType)

%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3, -3 for CLIMCAPS L3
%%         = +5 for ERA5, +2 for MERRA2
%%         = +6 for CMIP6, -6 for AMIP6
%%             1 for N, 2 for D

if ~isfield(driver,'iDebugRatesUseNWP')
  iDebugRatesUseNWP = -1;
else
  iDebugRatesUseNWP = driver.iDebugRatesUseNWP;
end

%ix = driver.iibin;
ix = driver.iLon;
iy = driver.iLat;

if driver.i16daytimestep < 0
  %% usual data
  load(driver.rateset.datafile)

  if settings.dataset == 30
    %% AMSU AMSU AMSU
    b_asc  = b_asc(:,:,1:13);  %% for some reason SARTA only has 13 chans
    b_desc = b_desc(:,:,1:13); %% for some reason SARTA only has 13 chans
    b_err_asc  = b_err_asc(:,:,1:13);  %% for some reason SARTA only has 13 chans
    b_err_desc = b_err_desc(:,:,1:13); %% for some reason SARTA only has 13 chans
  end

  switch driver.rateset.ocb_set

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'obs','tracegas'}
      % driver.rateset.rates = real(squeeze(b_obs(ix,:,2))');
      % driver.rateset.unc_rates = real(squeeze(b_err_obs(ix,:,2))');
      if driver.NorD > 0 
        driver.rateset.rates = real(squeeze(b_desc(ix,iy,:)));
      elseif driver.NorD < 0 
        driver.rateset.rates = real(squeeze(b_asc(ix,iy,:)));
      end

    iNnoise = +1;  %% use AIRS stability paper ideas
    iNnoise = -1;  %% use what successfully worked for May 2021 AIRS STM, and still works
    iNnoise = iNoiseType;
    if iNnoise < 0 
      if driver.NorD > 0
        driver.rateset.unc_rates = real(squeeze(b_err_desc(ix,iy,:)));  %% this is what was used for AIRS STM May 2021
      elseif driver.NorD < 0
        driver.rateset.unc_rates = real(squeeze(b_err_desc(ix,iy,:)));  %% this is what was used for AIRS STM May 2021, gives great results for daytime trends as well
        driver.rateset.unc_rates = real(squeeze(b_err_asc(ix,iy,:)));   %% equivalent of what was used for AIRS STM May 2021
      end
    else
      %% use the ideas from AIRS stability paper
      %{
        %% see driver_put_together_QuantileChoose_trends.m for conversion from mean_BT to NedT(T(\nu))
        addpath /home/sergio/KCARTA/MATLAB
        airs_noise = instr_chans2645('airs',2);
        load h2645structure.mat
        junk = squeeze(mean_BT(ix,iy,:));
        junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      %}
      junk = squeeze(airs_noiseTtrue(ix,iy,:));
      N = 0.01 * 12000; %% assume 12000 obs per 16 day interval per tile, and d(quantile) ~ 0.01 ... 
      driver.rateset.unc_rates = junk/sqrt(N);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'cal'
      % error('cal')   
      % driver.rateset.rates = real(squeeze(b_cal(ix,:,2))');
      % driver.rateset.unc_rates = real(squeeze(b_err_cal(ix,:,2))');

      % this is assuming I am reading in  dataset = 4, ocbset = +1
      %  driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
      %  driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
      if settings.dataset == 4 | settings.dataset == 5 | settings.dataset == 6 | settings.dataset == 7
        driver.rateset.rates = real(thesave.xtrend(:,ix));
        driver.rateset.unc_rates = real(thesave.xtrendErr(:,ix));

      % this is assuming I am reading in  dataset = 9, ocbset = +1
      %  driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '_2002_09_2022_08.mat']; %% co2/n2o/ch4 unchanging
      elseif settings.dataset == 9
        driver.rateset.rates = real(thesave.xtrend(:,ix));
        driver.rateset.unc_rates = real(thesave.xtrendErr(:,ix));

      else
        error('settings.dataset')
      end      

      % this is assuming I am reading in eg SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbinX/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat'];
      % dataset = 6,7,8,9, ocbset = +1
      % driver.rateset.rates = real(squeeze(b_cal_desc(ix,iy,:)));
      % driver.rateset.unc_rates = real(squeeze(b_cal_err_desc(ix,iy,:)));  %% this is what was used for AIRS STM May 2021

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'bias'
      error('bias')
      % driver.rateset.rates = real(squeeze(b_bias(ix,:,2))');
      % driver.rateset.unc_rates = real(squeeze(b_err_bias(ix,:,2))');

  end

elseif driver.i16daytimestep > 0
  if ~strfind(driver.anomalydatafile,'_tile_')
    anom = load(driver.rateset.datafile);
    iVersAnom = +1;
    if iVersAnom == 1
      iiAnomBins = (driver.anomalylatbin-1)*driver.anomalytimesteps + (1:driver.anomalytimesteps);
      fprintf(1,' get_rates.m : subsetting iiAnomBins between %6i and %6i \n',iiAnomBins(1),iiAnomBins(end));
      anom.avg16_btanom = (anom.btavgAnomFinal(:,iiAnomBins))';
    end
    get_anom_data_V0   %% orig, before Sept 2023 only had one anomaly time series
  else
    anom = load(driver.rateset.datafile);
    driver.rateset.rates = anom.btavgAnomFinal(:,driver.i16daytimestep);
    %anom_noise = load('noise_16day_avg_mission.mat');
    anom_noise = load('btn_avg.mat');
    anom_noise = squeeze(nanmean(squeeze(nanmean(anom_noise.btn_avg,1)),1));
    driver.rateset.unc_rates = anom_noise;
  end
end

%% disp('TESTING driver.rateset.unc_rates = 0.01')
%% driver.rateset.unc_rates = 0.01 * ones(size(driver.rateset.unc_rates));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now see if you need to overriode the spectral arates with reconstricted NWP rates
% see plot_profile_trends2.m

if iDebugRatesUseNWP > 0 
  %iDebugRatesUseNWP = -1;
  %iDebugRatesUseNWP = 52;
  if iDebugRatesUseNWP >= 30 & iDebugRatesUseNWP <= 39
    if driver.NorD == +1
      iDebugRatesUseNWP = 31;
    elseif driver.NorD == -1
      iDebugRatesUseNWP = 32;
    end
  elseif iDebugRatesUseNWP >= 50 & iDebugRatesUseNWP <= 59
    if driver.NorD == +1
      iDebugRatesUseNWP = 51;
    elseif driver.NorD == -1
      iDebugRatesUseNWP = 52;
    end
  elseif iDebugRatesUseNWP >= 60 & iDebugRatesUseNWP <= 69
    if driver.NorD == +1
      iDebugRatesUseNWP = 61;
    elseif driver.NorD == -1
      iDebugRatesUseNWP = 62;
    end
  end
  
  if iDebugRatesUseNWP > 0 & settings.ocb_set == 1
    disp('iDebugRatesUseNWP > 0 & settings.ocb_set == 1 so NOT loading in reconstrtcted sets, ,which hasve increasing CO2/CH4/N2O')
  end

  if iDebugRatesUseNWP > 0 & settings.ocb_set == 0
    iz = (iy-1)*72 + ix;
    fprintf(1,'  ---> get_rates.m settings.ocb_set == 0 AND iDebugRatesUseNWP=%2i ix/iy= %2i/%2i iz/iibin=%4i/%4i \n',[iDebugRatesUseNWP ix iy iz driver.iibin])
    if iz ~= driver.iibin
      error('driver.iibin ~= iz')
    end
  end
  if settings.ocb_set == 0 
    if iDebugRatesUseNWP == -1
      %% do nothing, stick to AIRS obs
    elseif iDebugRatesUseNWP == 31 | iDebugRatesUseNWP == 32 
      disp('  --> --> get_rates ... iDebugRatesUseNWP = 3 so use AIRS L3 reconstructed rates')
      if iDebugRatesUseNWP == 31
        load reconstructed_spectral_trends_nwp_night.mat
      elseif iDebugRatesUseNWP == 32
        load reconstructed_spectral_trends_nwp_day.mat
      end
      driver.rateset.rates = the_nwp_trends.airsL3(:,iz);
    elseif iDebugRatesUseNWP == 51 | iDebugRatesUseNWP == 52
      disp('  --> --> get_rates ... iDebugRatesUseNWP = 5 so use ERA5 reconstructed rates')
      if iDebugRatesUseNWP == 51
        load reconstructed_spectral_trends_nwp_night.mat
      elseif iDebugRatesUseNWP == 52
        load reconstructed_spectral_trends_nwp_day.mat
      end
      driver.rateset.rates = the_nwp_trends.era5(:,iz);
    elseif iDebugRatesUseNWP == 61 | iDebugRatesUseNWP == 62
      disp('  --> --> get_rates ... iDebugRatesUseNWP = 6 so use CMIP6 reconstructed rates')
      if iDebugRatesUseNWP == 61
        load reconstructed_spectral_trends_nwp_night.mat
      elseif iDebugRatesUseNWP == 62
        load reconstructed_spectral_trends_nwp_day.mat
      end
      driver.rateset.rates = the_nwp_trends.cmip6(:,iz);
    end
  end

end
