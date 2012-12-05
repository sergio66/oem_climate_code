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

% orig 6_97_97 stuff
driver.rateset.ncfile    = 'all_lagcor.mat';
driver.jacobian.filename = 'M_TS_jac_all.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   % how many gases we want to retrieve profiles
                                  % must be at least 1 (for water)
driver.jacobian.numlays     = 97;

%########################################################################
% larrabee used these in the file he gave me, but I get bad trace gas rates
% they give decent WV/T rates
%driver.oem.cov_filename     = 'cov_lls.mat'
%driver.oem.diag_only        = 0;
%driver.oem.lambda           = 1;
%driver.oem.lambda_qst       = 0.1;
%driver.oem.lambda_Q1        = 10;
%driver.oem.lambda_temp      = 1;
% ABOVE IS STROW, BELOW IS MINE %%%%%%%%%%%
% ABOVE IS STROW, BELOW IS MINE %%%%%%%%%%%

% copied rates from /home/sergio/MATLABCODE/RATES_TARO/MAT/
% driver.rateset.datafile  = 'AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Aug28_2012_span_07_2007_07_2012.mat';
% driver.rateset.datafile  = 'AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Sep11_2012_span_07_2007_07_2012.mat';
% driver.rateset.datafile  = 'AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_day_avgL1Brates_robust_Sep17_2012_span_09_2002_08_2012.mat';
% driver.rateset.datafile  = 'AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Sep17_2012_span_09_2002_08_2012.mat';

% good!!!!%%%
driver.rateset.datafile  = 'AIRS_MATFILES/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
% good!!!!%%%

% orig
%driver.oem.lambda_qst       = 0.1;
%driver.oem.lambda_Q1        = 100;
%driver.oem.lambda_temp      = 50;

% only make ** diags **
%driver.oem.lambda_qst       = [ones(1,6)]*0.1;
%driver.oem.lambda_Q1        = [ones(1,10) ones(1,30) ones(1,57)]*100;
%driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*1 ones(1,57)*1]*50;

% can individually ** tweak blocks ** [and thus duplicate driver.oem.lambda_temp = Z]
%driver.oem.lambda_qst       = [10 10 10 10 10 10                              ones(1,1)*(-9999)]*15;
%driver.oem.lambda_Q1        = [ones(1,10)*10.0 ones(1,30)*5.0 ones(1,57)*0.05 ones(1,1)*(-9999)]*1000;
%driver.oem.lambda_temp      = [ones(1,10)*10.0 ones(1,30)*5.0 ones(1,57)*0.05 ones(1,1)*(-9999)]*500;

% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat
%driver.oem.lambda_qst       = 0.1;
%driver.oem.lambda_Q1        = 1;
%driver.oem.lambda_temp      = 5;
%driver.oem.lambda_qst       = [0.1 0.1 0.1 0.1 0.10 0.1]*1e-1;
%driver.oem.lambda           = 1e-4;
%driver.oem.adjust_spectral_errorbars = 0.5;
% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat

%driver.oem.diag_only        = false;
%driver.oem.lambda_Q1        = 0.1;
%driver.oem.lambda_temp      = 0.01;
%driver.oem.lambda_qst       = [0.1 0.1 0.1 0.1 0.10 0.1]*1e-2;
%driver.oem.lambda           = 1e-2;
%driver.oem.adjust_spectral_errorbars = 0.5;

driver.oem.diag_only        = false;
driver.oem.lambda_Q1        = 5;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_qst       = [0.1 0.1 0.1 0.1 0.10 0.1]*1e-2;
driver.oem.lambda           = 1e+4;
driver.oem.adjust_spectral_errorbars = 0.5;

driver.oem.lambda_qst       = [1 1 1 1 1 1                                 ones(1,1)*(-9999)]*1e-4;
driver.oem.lambda_Q1        = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*10;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*1;
driver.oem.lambda           = 1e-2;

% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV
driver.oem.lambda_qst       = [1 1 1 1 1 1                                 ones(1,1)*(-9999)]*1e-4;
driver.oem.lambda_Q1        = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*50;
driver.oem.lambda_Q1        = [ones(1,40)*10 ones(1,30)*15 ones(1,27)*20 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*1.7 ones(1,27)*2.5 ones(1,1)*(-9999)]*5;
driver.oem.lambda           = 1e-2;
% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV


driver.oem.spectralcov_filename = ['../MakeERA_ratespectra/cov_spectra_yy_2009_latbin_' num2str(ix,'%02d') '.mat'];

fprintf(1,'%s \n',driver.rateset.datafile);

