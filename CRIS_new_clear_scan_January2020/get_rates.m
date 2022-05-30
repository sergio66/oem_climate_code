function driver = get_rates(driver)

if driver.i16daytimestep > 0
  junk = load(driver.rateset.datafile,'avg16_rtime');
  [yy,mm,dd] = tai2utcSergio(junk.avg16_rtime(driver.i16daytimestep));
  fprintf(1,'    <<<<< get_rates  fname=%s timestep=%3i %4i/%2i/%2i \n',driver.rateset.datafile,driver.i16daytimestep,yy,mm,dd)
end
  
ix = driver.iibin;
if driver.i16daytimestep < 0
  %% usual data
  load(driver.rateset.datafile)
  switch driver.rateset.ocb_set
    case 'bias'
     driver.rateset.rates = real(squeeze(b_bias(ix,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_bias(ix,:,2))');
    case 'cal'
     driver.rateset.rates = real(squeeze(b_cal(ix,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_cal(ix,:,2))');
    case {'obs','tracegas'}
     driver.rateset.rates = real(squeeze(b_obs(ix,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_obs(ix,:,2))');
  end

elseif driver.i16daytimestep > 0
  anom = load(driver.rateset.datafile);
  [mmjunk,nnjunk] = size(anom.avg16_btanom);
  if driver.i16daytimestep > mmjunk
    fprintf(1,'oh oh looking for timestep %3i but there are only %3i data points in anomdata file for latbin %2i \n',driver.i16daytimestep,mmjunk,ix)
    driver.rateset.rates = zeros(1305,1);
    driver.rateset.unc_rates = 0.01*ones(size(driver.rateset.rates));
    return
  else

    iSetType = 1;  % only use facx                              fir AIRS or CRIS
    iSetType = 2;  % use noise_16day_avg_mission.mat plus facx  MOSTLY AIRS ONLY
    iSetType = 3;  % use btn_avg.mat                 plus facx  MOSTLY AIRS ONLY
    iSetType = 4;  % use btn_avg.mat                 plus facx  for CRIS

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
    elseif iSetType == 4
      disp('  get_rates : unc iSetType = 4')
      if driver.i16daytimestep > 362
        iJunkTimeStep = 362;
      else
        iJunkTimeStep = driver.i16daytimestep;
      end
      %% this is harder :  spectra and for individual timesteps scale by sqrt(average counts x 16)
      %% note the relevant code (MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfsclust_make_lats40_16dayavg_noseasonal_16daybdry_cris.m) already accounts for average counts x 16 ,, so we have to appropriately divide by NeDT given by LLS
      snpp_noise = load('/home/sergio/MATLABCODE/oem_pkg_run/CRIS_new_clear_scan_January2020//snpp_cris_mean_noise_per40_lats.mat');
      counts_noise = load(['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/BDRY16INFO_CLR_CLO/bdry16info_' num2str(driver.iibin) '.mat']);
      anom_noise = snpp_noise.nedt(driver.iibin,:)/sqrt(counts_noise.p16.count(iJunkTimeStep));      
    end

    if iSetType < 4
      anom_noise = squeeze(anom_noise.btn_avg(ix,iJunkTimeStep,:));
    end
    lala = find(isnan(anom_noise) | isinf(anom_noise));
    anom_noise(lala) = facx;
    badones = [-9999 -9999];
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

%%%%%% this was just trying to use AIRS noise for CRIS, comment this out since now we do have the noise!!!
%%%%%% driver.rateset.unc_rates = driver.rateset.unc_rates(1:1305,:);
%%%%%% load closest_airs2645_to_cris1305_chan.mat
%%%%%% driver.rateset.unc_rates = driver.rateset.unc_rates(closest_airs2645_to_cris1305_chan,:)/4;  %% cris noise about 1/4 AIRS/IASI

%% this is actually pretty good, if you stick to OEM params from AIRS "build_cov_matrices.m"
%% disp('TESTING driver.rateset.unc_rates = 0.01')
%% driver.rateset.unc_rates = 0.01 * ones(size(driver.rateset.unc_rates));
