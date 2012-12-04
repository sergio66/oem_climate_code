function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

load AIRS_MATFILES/strow_stmNov2011_dobs20.mat
driver.jacobian.chanset = dobs20.jacobian.chanset;

% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +0.25;
driver.rateset.min = -0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orig 6_97_97 stuff
driver.rateset.ncfile    = 'all_lagcor.mat';
driver.jacobian.filename = 'M_TS_jac_all.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   % how many gases we want to retrieve profiles
                                  % must be at least 1 (for water)
driver.jacobian.numlays     = 97;

% TEST the SYNTHETIC CASE
% a) diag case, no damping at all .. get back tracegas rates, and biases are zeroish
%driver.rateset.datafile  = '../MakeERA_ratespectra/syntheticERArates.mat';
%driver.oem.adjust_spectral_errorbars = 0.001;
%driver.oem.lambda           = 1e-11;
%driver.oem.diag_only        = 1;

% b) some damping, and gives perfect retrievals
driver.rateset.datafile  = '../MakeERA_ratespectra/syntheticERArates.mat';
driver.oem.adjust_spectral_errorbars = 0.001;
driver.oem.lambda           = 1e-11;
driver.oem.diag_only        = 1;
driver.oem.diag_only        = false;
driver.oem.lambda_Q1        = 0.1;
driver.oem.lambda_temp      = 0.001;
driver.oem.lambda_qst       = [0.1 0.1 0.1 0.1 0.10 0.1]*1e-2;
driver.oem.lambda           = 1e-1;
driver.oem.spectralcov_filename = ['../MakeERA_ratespectra/cov_spectra_yy_2009_latbin_' num2str(ix,'%02d') '.mat'];
% TEST the SYNTHETIC CASE

fprintf(1,'%s \n',driver.rateset.datafile);

% g = dogoodchan;
% driver.jacobian.chanset = g;
% driver.jacobian.chanset = driver.jacobian.chanset(1:50:end);
