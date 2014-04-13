%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Need oem_pkg
addpath ../../oem_pkg
%---------------------------------------------------------------------------
% Doing debug?
driver.debug = false;
driver.debug_dir = '../Debug';
%---------------------------------------------------------------------------
% Perform OEM fit?
driver.oem.dofit = true;
%---------------------------------------------------------------------------
% Open debug file if desired
if driver.debug
  writelog('open');
end;
%---------------------------------------------------------------------------
driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
driver.jacobian.varname = 'M_TS_jac_all';
%---------------------------------------------------------------------------
% Fitting [obs][cal][biases], pick one
driver.rateset.ocb_set  = 'obs';
driver.jacobian.numlays    = 97;
%---------------------------------------------------------------------------
% Lag-1 correlation file
driver.rateset.ncfile   = 'all_lagcor.mat';
%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
%---------------------------------------------------------------------------
% Oem loops?  Just one if linear.
driver.oem.nloop = 1;
%---------------------------------------------------------------------------
% Get rate data
driver = get_rates(driver);
%---------------------------------------------------------------------------
% Q/A rates
% bad = find(driver.rateset.rates > driver.rateset.max);
% if length(bad) > 0
%   disp('resetting some bad input (max)');
%   driver.rateset.rates(bad) = driver.rateset.max;
% end
% bad = find(driver.rateset.rates < driver.rateset.min);
% if length(bad) > 0
%   disp('resetting some bad input (min)');

%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors
nc_cor = nc_rates(driver);
driver.rateset.unc_rates = nc_cor.*driver.rateset.unc_rates;
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
driver = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
ix = driver.iibin
%---------------------------------------------------------------------------
% Get jacobians
jac             = load(driver.jacobian.filename);
aux.m_ts_jac    = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Load in apriori
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;
[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.ix);
end
xb = xb./driver.qrenorm';
%---------------------------------------------------------------------------
% Form structure needed by rodgers.m
aux.xb       = xb;
%---------------------------------------------------------------------------
% Load in freq corrections
load Data/dbt_10year  % alldbt
driver.rateset.rates = driver.rateset.rates-alldbt(ix,:)'/10;
%---------------------------------------------------------------------------
% Modify with estimated error in freq = 0.01K
driver.rateset.unc_rates = ones(2378,1)*0.005;
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Save retrieval output
driver.filename = ['../Output/test' int2str(driver.iibin)];
save(driver.filename,'-struct','driver');
% Close debug file
if driver.debug
  writelog('close')
end
%---------------------------------------------------------------------------
% Some simple output
fprintf('Scalar Retrievals from OEM\n')
fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
fprintf(1,'O3    (%%)     %5.3f  +- %5.3f \n',100*driver.oem.finalrates(2),100*driver.oem.finalsigs(2));
fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
%---------------------------------------------------------------------------
% Plot Results
addpath Plotutils
plot_retrieval_latbins

