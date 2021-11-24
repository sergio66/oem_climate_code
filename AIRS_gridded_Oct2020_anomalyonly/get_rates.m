function driver = get_rates(driver)

if driver.i16daytimestep > 0
  %% usual data
  load(driver.rateset.datafile)
  switch driver.rateset.ocb_set
    case 'obs'
     driver.rateset.rates = real(bt_anom(:,driver.latlon.timestep));
     driver.rateset.unc_rates = 0.001 * ones(size(driver.rateset.rates));
    case 'cal'
     driver.rateset.rates = real(squeeze(b_cal(ix,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_cal(ix,:,2))');
    case {'bias'}
     driver.rateset.rates = real(squeeze(b_bias(ix,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_bias(ix,:,2))');
  end

elseif driver.i16daytimestep < 0
  error('this directory is for anoms, not trends')
end
