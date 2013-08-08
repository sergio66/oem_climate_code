clear all
load AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat

load /home/sergio/MATLABCODE/oem_pkg/Test/all_lagcor.mat
for ii = 1 : 36
  nc = (1+lagcor_obs_anom(ii,:))./(1-lagcor_obs_anom(ii,:));
  nc = sqrt(nc);
  nc_obs(ii,:) = nc;
  nc = (1+lagcor_bias_anom(ii,:))./(1-lagcor_bias_anom(ii,:));
  nc = sqrt(nc);
  nc_bias(ii,:) = nc;
end

load latbins.mat;
iibin = 1 : 36;
driver.rateset.obs_rates = real(squeeze(b_obs(iibin,:,2))');
driver.rateset.bias_rates = real(squeeze(b_bias(iibin,:,2))');

driver.rateset.obs_rates_err = real(squeeze(b_err_obs(iibin,:,2))');
driver.rateset.bias_rates_err = real(squeeze(b_err_bias(iibin,:,2))');

iibin = find(abs(save_lat) <= 30);
  ff = instr_chans('airs'); g = dogoodchan;

figure(1); 
  plot(ff(g),nanmean(driver.rateset.obs_rates(g,iibin)'),ff(g),nanstd(driver.rateset.obs_rates(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid
figure(2); 
  plot(ff(g),nanmean(driver.rateset.obs_rates(g,iibin)'),ff(g),nanmean(driver.rateset.obs_rates_err(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid
figure(3); 
  z = driver.rateset.obs_rates_err .* nc_obs';
  plot(ff(g),nanmean(driver.rateset.obs_rates(g,iibin)'),ff(g),nanmean(z(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid

figure(4); 
  plot(ff(g),nanmean(driver.rateset.bias_rates(g,iibin)'),ff(g),nanstd(driver.rateset.bias_rates(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid
figure(5); 
  plot(ff(g),nanmean(driver.rateset.bias_rates(g,iibin)'),ff(g),nanmean(driver.rateset.bias_rates_err(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid
figure(6); 
  z = driver.rateset.bias_rates_err .* nc_bias';
  plot(ff(g),nanmean(driver.rateset.bias_rates(g,iibin)'),ff(g),nanmean(z(g,iibin)'),'r')
  axis([625 2800 -0.1 +0.1]); grid
