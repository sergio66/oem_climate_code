function nc_cor =  nc_rates(driver);
% Correct rate uncertainties with 1-lag autocorrelation from ncfile

nc_cor = ones(size(driver.rateset.rates));
disp('LAGCOR = 1 LAGCOR = 1 LAGCOR = 1LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1 LAGCOR = 1')
return;

load(driver.rateset.ncfile)
iibin = driver.iibin;
switch driver.rateset.ocb_set
  case 'bias'
    %nc = (1+lagcor_bias_anom(iibin,:))./(1-lagcor_bias_anom(iibin,:));
    nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
    nc = sqrt(nc);
  case 'cal'
    %nc = (1+lagcor_cal_anom(iibin,:))./(1-lagcor_cal_anom(iibin,:));
    nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
    nc = sqrt(nc);
  case {'obs','tracegas'}
    nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
    nc = sqrt(nc);
end
nc_cor = real(nc');

