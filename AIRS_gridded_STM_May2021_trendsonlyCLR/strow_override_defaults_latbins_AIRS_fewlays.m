function [driver,aux,topts] = strow_override_defaults_latbins_AIRS_fewlays(driver0,iNlays_retrieve,topts0);

topts = topts0;
driver = driver0;

narginS = nargin;

check_settings  %% note this logic goes through defaults settings, and sees if some need to be overriden according to topts
check_driver    %% note this logic looks at some "miaow" settings; if driver has them things are fine ... else them them from "miaow"

aux.invtype = settings.invtype;

driver.topts0 = topts;
driver.topts  = topts;

%---------------------------------------------------------------------------
% Which latitude bin
if driver.ia_OorC_DataSet_Quantile(1) <= 1
  ix = driver.iibin;
else
  ix = driver.anomalylocation;
end

%---------------------------------------------------------------------------
% Fitting [obs][cal][bias], pick one
if settings.ocb_set == -1
  driver.rateset.ocb_set  = 'bias';
  disp('bias')
elseif settings.ocb_set == 0
  driver.rateset.ocb_set  = 'obs';
  disp('obs')
elseif settings.ocb_set == +1
  driver.rateset.ocb_set  = 'cal';
  disp('cal')
elseif abs(settings.ocb_set) == +2
  driver.rateset.ocb_set  = 'anom';
  driver.rateset.ocb_set  = 'obs';
  disp('anomaly')
elseif abs(settings.ocb_set) > 2
  settings.ocb_set
  error('incorrect settings.ocb_set')
end
%---------------------------------------------------------------------------
% Raw rate data file        
%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3, -3 for CLIMCAPS L3
%%         = +5 for ERA5, +2 for MERRA2
%%         = +6 for CMIP6, -6 for AMIP6
%%             1 for N, 2 for D
[driver,settings] = set_driver_rateset_datafile(driver,settings);

fprintf(1,'[settings.ocb_set settings.descORasc driver.i16daytimestep settings.dataset] = %3i %3i %3i %3i \n',[settings.ocb_set settings.descORasc driver.i16daytimestep settings.dataset])
fprintf(1,' <<< driver.rateset.datafile >>> = %s \n',driver.rateset.datafile)

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../oem_pkg/Test/all_lagcor.mat';
driver.rateset.ncfile   = driver.rateset.datafile;

%driver
% Get rate data, do Q/A elsewhere
%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3, -3 for CLIMCAPS L3
%%         = +5 for ERA5, +2 for MERRA2
%%         = +6 for CMIP6, -6 for AMIP6
%%             1 for N, 2 for D
driver = get_rates(driver,settings,settings.iNoiseType);  %% this gets spectral rates (driver.rateset.rates), and uncertainty (driver.rateset.unc_rates)
if driver.removeEmisTrend > 0
  emiseffect = get_emissivity_trends(driver);
  plot(1:2645,driver.rateset.rates,1:2645,emiseffect); plotaxis2; hl = legend('rates','emiss trend','location','best','fontsize',10);
  driver.rateset.rates = driver.rateset.rates - emiseffect;
end

%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
[driver,iVersJac,iOldORNew,iXJac,topts] = set_driver_jacfile(driver,settings,topts);

set_the_jacobians  %% sets structure "jac" and m_ts_jac

aux.m_ts_jac = m_ts_jac;
aux.f        = jac.f;
aux.spres    = jac.spres;
aux.stemp    = jac.stemp;
aux.nlays    = jac.nlays; %% 100
aux.plays    = jac.plays; %% 1:100
aux.ptemp    = jac.ptemp; %% 1:100
aux.gas_1    = jac.gas_1; %% 1:100
aux.gas_3    = jac.gas_3; %% 1:100
aux.navg     = jac.navg;  %% N
aux.pavg     = jac.pavg;  %% 1 : N
aux.tavg     = jac.tavg;  %% 1 : N
aux.qavg     = jac.qavg;  %% 1 : N
aux.oavg     = jac.oavg;  %% 1 : N
aux.trop_P   = jac.trop_P;
aux.trop_ind = jac.trop_ind;

f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Good channel set
lenrates = length(driver.rateset.rates);

iChSet = 2; %% new chans
iChSet = 1; %% old chans (default)
iChSet = 4; %% new chans + Tonga (high alt)
iChSet = 5; %% new chans + laserlines
iChSet = 3; %% new chans, but no CFC11
iChSet = topts.iChSet;

ch = find_the_oem_channels(f,lenrates,settings.numchan,settings.chan_LW_SW,iChSet);

driver.topts.iChSet = iChSet;
driver.jacobian.chanset = ch;
%---------------------------------------------------------------------------
% Apriori file
%driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
%xb = load(driver.oem.apriori_filename,'apriori');
%xb = xb.apriori;

%xb = zeros(296,1);
xb = zeros(driver.jacobian.wvjaclays_offset + iNlays_retrieve*3,1);

if (abs(settings.set_tracegas) ~= 1) & (settings.set_tracegas ~= 2)
  settings.set_tracegas
  error('incorrect setting for overriding xb tracegas values');
end

if settings.iFixTG_NoFit(1) > 0
  disp('oh oh you wanna get rid of a trace gas')
  junk = settings.iFixTG_NoFit;
  badjunk = find(junk < 1 | junk > max(driver.jacobian.scalar_i) - 1);
  if length(badjunk) > 0 
    junk
    driver.jacobian.scalar_i
    error('can only thrown out a few of first 1 .. driver.jacobian.scalar_i - 1 gases!!!!')
  else
    numthrow = length(junk);
    fprintf(1,'throwing out %3i trace gases from the fit \n',numthrow);
%{
    aux.jacobian.scalar_i_allgases = 1:driver.jacobian.scalar_i;
    aux.jacobian.water_i_allgases  = driver.jacobian.water_i;
    aux.jacobian.temp_i_allgases   = driver.jacobian.temp_i;
    aux.jacobian.ozone_i_allgases  = driver.jacobian.ozone_i;

    driver.jacobian.scalar_i = 1:driver.jacobian.scalar_i - numthrow;
    driver.jacobian.water_i  = driver.jacobian.water_i - numthrow;
    driver.jacobian.temp_i  = driver.jacobian.temp_i - numthrow;
    driver.jacobian.ozone_i  = driver.jacobian.ozone_i - numthrow;
%}
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if settings.iFixTz_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
  iFixTz_NoFit = +1;    %%% LARABBEE LIKES THIS TURNED OFF ie keep spectra as is, just read in ERA anom and proceed
end
iZeroTVers = 0; %%% use my fit to sarta calcs, as a proxy to ERA T anomalies
iZeroTVers = 1; %%% use raw ERA T anomalies, and do the averaging here on the fly
iZeroTVers = 2; %%% use raw ERA T anomalies as saved in era_ptempanom.mat (see compare_era_anomaly_from_fit_and_model.m)

set_zeroT_nofit

%%%%%%%%%%%%%%%%%%%%%%%%%

if settings.iFixO3_NoFit >= 0 & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
  iFixO3_NoFit = settings.iFixO3_NoFit;    %%% LARABBEE LIKES THIS TURNED OFF ie keep spectra as is, just read in ERA anom and proceed
end

iZeroO3Vers = 0; %%% use my fit to sarta calcs, as a proxy to ERA T anomalies
iZeroO3Vers = 1; %%% use raw ERA T anomalies, and do the averaging here on the fly
iZeroO3Vers = 2; %%% use raw ERA T anomalies as saved in era_ptempanom.mat (see compare_era_anomaly_from_fit_and_model.m)

set_zeroO3_nofit

%%%%%%%%%%%%%%%%%%%%%%%%%

if settings.iFixWV_NoFit >= 0 & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
  iFixWV_NoFit = settings.iFixWV_NoFit;    %%% LARABBEE LIKES THIS TURNED OFF ie keep spectra as is, just read in ERA anom and proceed
end

iZeroWVVers = 0; %%% use my fit to sarta calcs, as a proxy to ERA T anomalies
iZeroWVVers = 1; %%% use raw ERA T anomalies, and do the averaging here on the fly
iZeroWVVers = 2; %%% use raw ERA T anomalies as saved in era_ptempanom.mat (see compare_era_anomaly_from_fit_and_model.m)

set_zeroWV_nofit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is where topts.set_era5_cmip6_airsL3 is used
set_apriori_ERA5_MERRA2_or_AIRSL3_MLS_geophysical   %% can set a priori to the MLS, or to the ERA5, or to the MERRA2, or to the AIRS L3, rates

iAdjLowerAtmWVfrac = 0;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!
iAdjLowerAtmWVfrac = 1;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!
iAdjLowerAtmWVfrac = topts.iAdjLowerAtmWVfrac;      %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!
driver.co2adj_ESRL = -9999;
if driver.ia_OorC_DataSet_Quantile(1) == 0
  %% these are obs, so use CO2 trends!!
  set_CO2_CH4_N2O_ESRL                                %% can set CO2/CH4/N2O to ESRL rates, can also set low atm dWV/dt using Isaac Held delta(RH)=0
%elseif driver.ia_OorC_DataSet_Quantile(1) == 1 & driver.ia_OorC_DataSet_Quantile(4) == 5
%  %% these are ERA5 simulations with CO2, so use CO2 trends!!
%  set_CO2_CH4_N2O_ESRL                                %% can set CO2/CH4/N2O to ESRL rates, can also set low atm dWV/dt using Isaac Held delta(RH)=0
%elseif driver.ia_OorC_DataSet_Quantile(1) == 1 & driver.ia_OorC_DataSet_Quantile(4) ~= 5
%  disp('using simulated spectra which are not ERA5 sims, so do not have CO2 in them')
elseif driver.ia_OorC_DataSet_Quantile(1) == 1
  %% these are ERA5 simulations with CO2, so use CO2 trends!!
  set_CO2_CH4_N2O_ESRL                                %% can set CO2/CH4/N2O to ESRL rates, can also set low atm dWV/dt using Isaac Held delta(RH)=0
end

%{
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
disp('WV xb WV xb WV xb')
  xb(driver.jacobian.water_i) = +0.01/2;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    
iSergioCO2 = -1;  %% assume ESRL CO2/CH4 rates
iSergioCO2 = +1;  %% fit for CO2/CH4 rates
iSergioCO2 = settings.iSergioCO2;
if iSergioCO2 > 0 & settings.ocb_set == 1
  disp('iSergioCO2 = +1 so RETRIEVE trace gases from ERA calc!!!!')
  xb(1:6) = 0;  %%% sergio, float the CO2 rates !!!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
elseif iSergioCO2 > 0 & settings.ocb_set == 0
  disp('iSergioCO2 = +1 so RETRIEVE trace gases from obs!!!!')
  xb(1) = co2x;
  xb(2) = n2ox;
  xb(3) = ch4x;
  xb(1:6) = 0;  %%% sergio, float the CO2 rates !!!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end

[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.iibin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%

set_Tz_O3z_noFit

%%%%%%%%%%%%%%%%%%%%%%%%%
% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';

%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
driver.oem.xb = xb;  %% note this is un-normalized xb

%---------------------------------------------------------------------------
if abs(settings.offsetrates) ~= 1
  settings.offsetrates
  error('incorrect setting for overriding obs rates by constant');
elseif settings.offsetrates > 0 & driver.i16daytimestep < 0
   disp('Offset linear rates (obs or cal or bias) by constant to get sensitivity to AIRS BT drift')
   driver.rateset.rates = driver.rateset.rates + 0.01;
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 0 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (obs) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 1 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (cal) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
end

%------------------------------------------------------------------------
if abs(settings.addco2jacs) ~= 1
  settings.addco2jacs
  error('incorrect setting for adding co2jacs to ERA calcrates');
end
if strcmp(driver.rateset.ocb_set,'cal') & settings.addco2jacs > 0
  disp('Offset rates to add in CO2 jacs to ERA rates, to pretend there is CO2')
  jac = load('co2_kcartajac.mat');
  %  haha = driver.rateset.rates;
  %  baba = jac.co2jac(ix,:)';
  %  whos haha baba
  driver.rateset.rates = driver.rateset.rates + jac.co2jac(ix,:)';
end

%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors or add to above
if driver.i16daytimestep < 0
  if driver.topts.dataset == 6
    driver.rateset.unc_rates = driver.rateset.unc_rates/4;  %% 6 years of data ==> large unc
  elseif driver.topts.dataset == 8
    driver.rateset.unc_rates = driver.rateset.unc_rates/6;  %% 6 years of data ==> large unc
  end
  nc_cor = nc_rates(driver);
  driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;    %% THIS IS AN ARRAY
  %driver.rateset.unc_rates = driver.rateset.unc_rates *sqrt(1/8); %% this accounts for counting .... 
end

%driver.rateset.unc_rates = 0.001 * ones(size(driver.rateset.unc_rates));   %% THIS IS AN ARRAY

%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>
if settings.obs_corr_matrix > 0
  addpath /home/sergio/MATLABCODE
  thecov = load('/home/sergio/MATLABCODE/oem_pkg_run/Simulate_Calcs/thecov_clear');
  [f2645,i2645] = map_2834_to_2645;
  junk = thecov.thecov; junk = junk(i2645,i2645);

  %% see https://www.mathworks.com/help/finance/corr2cov.html
  %%junk = corr2cov(driver.rateset.unc_rates,junk);
  %% new, fast
  junk0 = junk;
  junk0 = diag(driver.rateset.unc_rates) * junk0 * diag(driver.rateset.unc_rates);

  %% new,slow
  %[mm,nn] = size(junk);
  %for ii = 1 : mm
  %  for jj = 1 : nn
  %    junk(ii,jj) = junk(ii,jj) * driver.rateset.unc_rates(ii) * driver.rateset.unc_rates(jj);
  %  end
  %end
  %sum(sum(junk-junk0))/nn/nn
  junk = junk0;  
  fprintf(1,'numchans, rank, mean rank, condition number of obs cov matrix = %6i %6i %8.6e %8.6e \n',nn,rank(junk),rank(junk)/length(junk),cond(junk))

  %% old
  %plot(f2645,driver.rateset.unc_rates.^2,'r',thecov.fairs(i2645),diag(junk))
  %%keyboard_nowindow
  %slope = (2645+1);
  %xind = 1:2645; diagind=slope*(xind-1)+1;
  %junk(diagind) = driver.rateset.unc_rates.^2;

  driver.rateset.unc_rates = junk;  %% THIS IS A SQUARE MATRTIX
end
%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>

% Modify with estimated error in freq + regress errors 
%driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates.*nc_cor;
%load corr_errors
%load 15yr_corr_errors_2314chans

%---------------------------------------------------------------------------
% Do rate Q/A (empty for now)
%---------------------------------------------------------------------------

if settings.resetnorm2one == +1
  %% we already cleared jac
  %aux

  disp(' ')
  disp('settings.resetnorm2one == +1')
  disp('resetting all jabian norms to 1 ONE UNO UN MOJA')
  disp('resetting all jabian norms to 1 ONE UNO UN MOJA')
  disp('resetting all jabian norms to 1 ONE UNO UN MOJA')
  disp('resetting all jabian norms to 1 ONE UNO UN MOJA')
  disp(' ')

  boo1 = driver.qrenorm;
  boo2 = aux.xb; 
  %whos boo1 boo2
  %[driver.qrenorm'  aux.xb]

  for ii = 1 : length(boo1)
    aux.m_ts_jac(:,ii) = aux.m_ts_jac(:,ii)/driver.qrenorm(ii);
  end;
  aux.xb = aux.xb.*driver.qrenorm';
  driver.qrenorm = ones(size(driver.qrenorm));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
build_cov_matrices  %% iLatX controls the "width" of the polar regions, which have different covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(13); colormap jet; imagesc(log10(abs(aux.m_ts_jac'))); caxis([-12 0]); colorbar 
%disp('here 3'); pause

%driver.jacobian.thin = aux.m_ts_jac;
%keyboard_nowindow

driver.iaSequential = topts.iaSequential;
driver.dataset      = topts.dataset;

driver.topts = topts;
