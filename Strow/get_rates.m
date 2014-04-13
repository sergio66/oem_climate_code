function driver = get_rates(driver)

ix = driver.iibin;
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

