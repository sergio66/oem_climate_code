if exist('iFixTz_NoFit','var') & strcmp(driver.rateset.ocb_set,'obs')

  disp(' >>>>> need to account for cal_T_rates by eg subbing in T(z,t) from fits to cal or raw ERA anoms >>>>> ')
  xb(driver.jacobian.temp_i) = cal_T_rates;

  aux.FixTz_NoFit = cal_T_rates;
  aux.orig_water_i = driver.jacobian.water_i;
  aux.orig_temp_i = driver.jacobian.temp_i;
  aux.orig_ozone_i = driver.jacobian.ozone_i;

  noT = setdiff(1:length(xb),driver.jacobian.temp_i);

  driver.jacobian.ozone_i =  driver.jacobian.temp_i;
  driver.jacobian.temp_i = [];

  m_ts_jac = m_ts_jac(:,noT);
  aux.m_ts_jac    = m_ts_jac;

  xb = xb(noT);
  driver.rateset.rates = driver.rateset.rates - spectra_due_to_T_jac;   %% LARRABEE DOES NOT WANT THIS
  driver.qrenorm = driver.qrenorm(noT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('iFixO3_NoFit','var') & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))

  disp(' >>>>> need to account for cal_O3_rates by eg subbing in O3(z,t) from fits to cal or raw ERA anoms >>>>> ')
  xb(driver.jacobian.ozone_i) = cal_O3_rates;

  aux.FixO3_NoFit = cal_O3_rates;
  aux.orig_water_i = driver.jacobian.water_i;
  aux.orig_temp_i = driver.jacobian.temp_i;
  aux.orig_ozone_i = driver.jacobian.ozone_i;

  noO3 = setdiff(1:length(xb),driver.jacobian.ozone_i);

  driver.jacobian.ozone_i =  [];

  m_ts_jac = m_ts_jac(:,noO3);
  aux.m_ts_jac    = m_ts_jac;

  xb = xb(noO3);
  driver.rateset.rates = driver.rateset.rates - spectra_due_to_O3_jac;   %% LARRABEE DOES NOT WANT THIS
  driver.qrenorm = driver.qrenorm(noO3);
end
