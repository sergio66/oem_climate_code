function driver = get_rates(driver,iNoiseType)

%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3
%%         = +5 for ERA5
%%         = +6 for CMIP6
%%             1 for N, 2 for D

if ~isfield(driver,'iDebugRatesUseNWP')
  iDebugRatesUseNWP = -1;
else
  iDebugRatesUseNWP = driver.iDebugRatesUseNWP;
end


%ix = driver.iibin;
ix = driver.iLon;
iy = driver.iLat;
%keyboard_nowindow

if driver.i16daytimestep < 0
  %% usual data
  load(driver.rateset.datafile)
  switch driver.rateset.ocb_set
    case 'bias'
      error('bias')
      % driver.rateset.rates = real(squeeze(b_bias(ix,:,2))');
      % driver.rateset.unc_rates = real(squeeze(b_err_bias(ix,:,2))');
    case 'cal'
      error('cal')
      % driver.rateset.rates = real(squeeze(b_cal(ix,:,2))');
      % driver.rateset.unc_rates = real(squeeze(b_err_cal(ix,:,2))');
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

  end

elseif driver.i16daytimestep > 0
  anom = load(driver.rateset.datafile);
  [mmjunk,nnjunk] = size(anom.avg16_btanom);
  if driver.i16daytimestep > mmjunk
    fprintf(1,'oh oh looking for timestep %3i but there are only %3i data points in anomdata file for latbin %2i \n',driver.i16daytimestep,mmjunk,ix)
    driver.rateset.rates = zeros(2645,1);
    driver.rateset.unc_rates = 0.01*ones(size(driver.rateset.rates));
    return
  else

    iSetType = 1;  % only use facx
    iSetType = 2;  % use noise_16day_avg_mission.mat plus facx
    iSetType = 3;  % use btn_avg.mat                 plus facx

    facx = 0.01; %%% till june 21, 2019

    %% easiest : use same noise for all channels
    facx = 0.1;    %%% test
    facx = 0.001;  %%% test --> has big problems around 2007 (some anomalies are not fit)
    facx = 0.005;  %%% test
    facx = 0.0025; %%% test --> maybe best YAY YAY YAY

    if iSetType == 2
      disp('  get_rates : unc iSetType = 2')
      %% this is next step : spectral, but for any time step scale by sqrt(average counts x 16)
	anom_noise = load('noise_16day_avg_mission.mat');
      anom_noise = anom_noise.btn_av(ix,:);
    elseif iSetType == 3
      disp('  get_rates : unc iSetType = 3')
      %% this is harder :  spectra and for individual timesteps scale by sqrt(average counts x 16)
      anom_noise = load('btn_avg.mat');
      if driver.i16daytimestep > 362
        iJunkTimeStep = 362;
      else
        iJunkTimeStep = driver.i16daytimestep;
      end
    end

    anom_noise = squeeze(anom_noise.btn_avg(ix,iJunkTimeStep,:));
    lala = find(isnan(anom_noise) | isinf(anom_noise));
    anom_noise(lala) = facx;
    if (length(intersect(badones(:,1),driver.iibin)) == 1) & (length(intersect(badones(:,2), driver.i16daytimestep)) == 1)
      disp('bad one bad one ... x2 for noise bad voodoo daddy')
      anom_noise = anom_noise * 2;
    end

    switch driver.rateset.ocb_set
      case {'obs','tracegas'}
        driver.rateset.rates = (real(anom.avg16_btanom(driver.i16daytimestep,:)))';        
        driver.rateset.unc_rates = 0.01*ones(size(driver.rateset.rates)); %%% till june 21, 2019
        if (iSetType == 1)    
          driver.rateset.unc_rates = facx*ones(size(driver.rateset.rates)); %%% after june 21, 2019, yay use facx = 0.0025
        elseif (iSetType >= 2)
          driver.rateset.unc_rates = anom_noise;                            %%% after june 24, 2019
        end
      case {'cal'}
        driver.rateset.rates = (real(anom.avg16_btanom(driver.i16daytimestep,:)))';
        driver.rateset.unc_rates = 0.01*ones(size(driver.rateset.rates)); %%% till june 21, 2019
        if (iSetType == 1)    
          driver.rateset.unc_rates = facx*ones(size(driver.rateset.rates)); %%% after june 21, 2019, yay use facx = 0.0025
        elseif (iSetType >= 2)
          driver.rateset.unc_rates = anom_noise;                            %%% after june 24, 2019
        end
      otherwise
        error('anom for obs + cal only')
    end  
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
  
  if iDebugRatesUseNWP > 0
    iz = (iy-1)*72 + ix;
    fprintf(1,'  ---> get_rates.m iDebugRatesUseNWP=%2i ix/iy= %2i/%2i iz/iibin=%4i/%4i \n',[iDebugRatesUseNWP ix iy iz driver.iibin])
    if iz ~= driver.iibin
      error('driver.iibin ~= iz')
    end
  end
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

