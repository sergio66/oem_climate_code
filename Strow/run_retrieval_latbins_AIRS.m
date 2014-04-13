%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
driver.debug = false;
%---------------------------------------------------------------------------
% Need oem_pkg
addpath ../../oem_pkg
%---------------------------------------------------------------------------
% Apriori file
driver.oem.apriori_filename = 'apriori_zero';
%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
driver.jacobian.varname = 'M_TS_jac_all';
%---------------------------------------------------------------------------
% Lag-1 correlation file
driver.rateset.ncfile   = 'all_lagcor.mat';
%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.sarta_error = 0.0; 
driver.oem.sarta_error = 0.0;
%---------------------------------------------------------------------------
% Perform OEM fit
% Generallly true, but if you want to do LLS only, set this to false
driver.oem.dofit = true;
%---------------------------------------------------------------------------
% how many times are we looping
driver.oem.nloop = 1;
%---------------------------------------------------------------------------
% Do one of [obs][cal][biases] 
driver.rateset.ocb_set  = 'obs';
driver.jacobian.numlays    = 97;
%---------------------------------------------------------------------------
driver = strow_override_defaults_latbins_AIRS(driver);  %% this is YOUR settings
ix = driver.iibin







%---------------------------------------------------------------------------
% Get rate data
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
% Q/A rates
% bad = find(driver.rateset.rates > driver.rateset.max);
% if length(bad) > 0
%   disp('resetting some bad input (max)');
%   driver.rateset.rates(bad) = driver.rateset.max;
% end
% bad = find(driver.rateset.rates < driver.rateset.min);
% if length(bad) > 0
%   disp('resetting some bad input (min)');

% Get jacobians
jac = load(driver.jacobian.filename);
aux.m_ts_jac = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Get the observed rate errors, which ultimately are a combination of
%  (a) get_rates.m : driver.rateset.unc_rates = real(squeeze(b_err_obs(ix,:,2))');
%  (b) nc_rates.m  : ncerrors = real(nc' .* driver.rateset.unc_rates); where nc is 1, or comes from a file
aux.ncerrors = nc_rates(driver);

% Observation errors : can adjust them with a scalar 
%ncerrors = aux.ncerrors * driver.oem.adjust_spectral_errorbars;

%se_errors.ncerrors = ncerrors;




% xset are the apriori rates, in use units eg ppm/yr, K/yr
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;
[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.ix);
end
xb = xb./driver.qrenorm';

% Form structure needed by rodgers.m
aux.xb       = xb;
aux.ncerrors = ncerrors;
%---------------------------------------------------------------------------
% Load in freq corrections
load Data/dbt_10year  % alldbt
driver.rateset.rates = driver.rateset.rates-alldbt(ix,:)'/10;
%---------------------------------------------------------------------------
% Modify with estimated error in freq = 0.01K
driver.rateset.unc_rates = ones(2378,1)*0.005;
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%% Save retrieval output
save(driver.filename,'-struct','driver');

%% Close debug file
if driver.debug
  writelog('close')
end


% PLOTTING
% 
aa = load(driver.rateset.datafile);

figure(1);
  g  = dogoodchan; ff = instr_chans;
  g1 = driver.jacobian.chanset_used;
  plot(ff(g1),driver.rateset.rates(g1),ff(g1),squeeze(aa.b_obs(ix,g1,2)),'k',ff(g1),driver.oem.fit(g1),'r.-'); grid; 
    axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','DEFINITELY OBS','fits'); set(hl,'fontsize',10)
figure(2);
  g  = dogoodchan; ff = instr_chans;
  g1 = driver.jacobian.chanset_used;
  plot(ff(g1),driver.rateset.rates(g1),ff(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r'); grid; axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','fits'); set(hl,'fontsize',10)
% if length(driver.oem.finalrates) == 200
%   figure(3)
%   plot(driver.oem.finalrates(7:103),1:97,driver.oem.finalrates(104:200),1:97,'r')
%   set(gca,'ydir','reverse');
%   title('AIRS (b) : WV frac/yr (r) T K/yr'); grid
% end

plot_retrieval_latbins

