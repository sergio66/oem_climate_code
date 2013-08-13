clear all

ff = instr_chans('airs'); g = dogoodchan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% we have 36 latbins serial correlations
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
iiBin36 = 1 : 36;
iiBin36 = find(abs(save_lat) <= 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% these are -90:-60 -60:-30 -30:+30 +30:+60 +60:+90
load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Sep17_2012_span_09_2002_08_2012.mat
load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_day_save_lat_Sep17_2012.mat

iiBin = 1 : 5;
driver.rateset.obs_rates = real(squeeze(b_obs(iiBin,:,2))');
driver.rateset.bias_rates = real(squeeze(b_bias(iiBin,:,2))');

driver.rateset.obs_rates_err = real(squeeze(b_err_obs(iiBin,:,2))');
driver.rateset.bias_rates_err = real(squeeze(b_err_bias(iiBin,:,2))');

iiBin = 3;
figure(1); 
  plot(ff(g),driver.rateset.obs_rates(g,iiBin),ff(g),driver.rateset.obs_rates_err(g,iiBin),'r')
  axis([625 2800 -0.1 +0.1]); grid
%figure(2); 
%  plot(ff(g),driver.rateset.obs_rates(g,iiBin),ff(g),nanmean(driver.rateset.obs_rates_err(g,iiBin)'),'r')
%  axis([625 2800 -0.1 +0.1]); grid
figure(3); 
  z = driver.rateset.obs_rates_err(:,iiBin) .* (nanmean(nc_obs(iiBin36,:)))';
  plot(ff(g),driver.rateset.obs_rates(g,iiBin),ff(g),z(g),'r')
  axis([625 2800 -0.1 +0.1]); grid

figure(4); 
  plot(ff(g),driver.rateset.bias_rates(g,iiBin),ff(g),driver.rateset.bias_rates_err(g,iiBin),'r')
  axis([625 2800 -0.1 +0.1]); grid
%figure(5); 
%  plot(ff(g),nanmean(driver.rateset.bias_rates(g,iiBin)'),ff(g),nanmean(driver.rateset.bias_rates_err(g,iiBin)'),'r')
%  axis([625 2800 -0.1 +0.1]); grid
figure(6); 
  z = driver.rateset.bias_rates_err(:,iiBin) .* (nanmean(nc_bias(iiBin36,:)))';
  plot(ff(g),driver.rateset.bias_rates(g,iiBin),ff(g),z(g),'r')
  axis([625 2800 -0.1 +0.1]); grid
