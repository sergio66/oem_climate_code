%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Select the latitude bin
driver.iibin = JOB;
ix = driver.iibin;
%---------------------------------------------------------------------------
% Need oem_pkg
addpath ../../oem_pkg
%---------------------------------------------------------------------------
% Doing debug?
driver.debug = false;
driver.debug_dir = '../Debug';

% Open debug file if desired
if driver.debug
  writelog('open');
end;
%---------------------------------------------------------------------------
% Perform OEM fit?
driver.oem.dofit = true;

% Oem loops?  Just one if linear.
driver.oem.nloop = 1;
%---------------------------------------------------------------------------
% Raw rate data file
driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';

% Fitting [obs][cal][bias], pick one
driver.rateset.ocb_set  = 'obs';

% Good channel set
load /asl/s1/rates/clear/good_chanset.mat 
driver.jacobian.chanset = chanset;

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../../oem_pkg/Test/all_lagcor.mat';

% Get rate data, do Q/A elsewhere
driver = get_rates(driver);
%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:6;
driver.jacobian.water_i  = 7:103;
driver.jacobian.temp_i   = 104:200;
driver.jacobian.numlays  = 97;

% Get jacobians
jac             = load(driver.jacobian.filename);
aux.m_ts_jac    = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Apriori file
driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;
[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.ix);
end
% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';
%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
driver = strow_override_defaults_latbins_AIRS(driver);
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

