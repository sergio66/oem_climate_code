function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

% good!!!!%%%
driver.rateset.datafile  = '/asl/s1/sergio/Rate_data/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
% good!!!!%%%

load /asl/s1/sergio/Rate_data/strow_stmNov2011_dobs20.mat
driver.jacobian.chanset = dobs20.jacobian.chanset;

% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +0.25;
driver.rateset.min = -0.25;

% orig 6_97_97 stuff
driver.rateset.ncfile    = 'all_lagcor.mat';
%driver.rateset.ocb_set  = 'bias';
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


% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV
driver.oem.lambda_qst       = [1 1 1 1 1 1                                 ones(1,1)*(-9999)]*1e-4;
driver.oem.lambda_Q1        = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*50;
driver.oem.lambda_Q1        = [ones(1,40)*10 ones(1,30)*15 ones(1,27)*20 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*1.7 ones(1,27)*2.5 ones(1,1)*(-9999)]*5;
driver.oem.lambda           = 1e-2;
driver.jacobian.filename = 'M_TS_jac_all.mat';
% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV

%driver.rateset.datafile  = ...
%  '/home/sergio/MATLABCODE/RATES_TARO/MAT/xoverocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Dec13_2012_span_09_2002_08_2012.mat';
%  driver.oem.spectralcov_filename = ['../MakeERA_ratespectra/cov_spectra_yy_2009_latbin_' num2str(ix,'%02d') '.mat'];
%  driver.oem.adjust_spectral_errorbars = 1.0;  %% notice how I need to increase this when I use spectral_cov
%  driver.oem.lambda           = 1e2;   %% pretty good
%  driver.oem.lambda_temp      = 1e2;   %% pretty good
%  driver.oem.lambda_Q1        = 1e4;   %% pretty good
%  driver.oem.lambda           = 2.5e+2;   %% very good
%  driver.oem.lambda_temp      = 2000;     %% very good
%  driver.oem.lambda_Q1        = 250000.0; %% very good

%driver.jacobian.chanset = driver.jacobian.chanset(driver.jacobian.chanset < 1406);  %% this is roughly 1305 cm-1
%driver.jacobian.chanset = driver.jacobian.chanset(driver.jacobian.chanset < 904);  %% this is roughly 961 cm-1

driver.jacobian.filename = '/asl/s1/sergio/Rate_data/Aug2013/all_kcarta_HARTMANNCO2_jacs_NOV02_2012_HARTMANN_6_97_97.mat';
driver.jacobian.filename = 'M_TS_jac_all.mat';
driver.jacobian.filename = '/asl/s1/sergio/Rate_data/Aug2013/all_kcarta_UMBCCO2_jacs_NOV02_2012_UMBCCO2_6_97_97.mat';

  %% but if smoothVSregularization = 'c'
  %% note these are PHYSICAL units eg sigma_temp_stratVALUE is in KELVIN (or KELVIN/YR)
  %% so eg Andy Tangborn finds after normalization driver.oem.sigma.temp_strat_VALUE  = 4;
  %% is a "good value" which means we default set it to 4*0.01 (0.01 is qrenorm for T(z))
  %% strat temp
    driver.oem.sigma.temp_strat_VALUE  = 4.0*0.01;       %% sigsqr
    driver.oem.sigma.temp_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.temp_strat_BOTLAY = 49;  %% stop layer
  %% trop temp
    driver.oem.sigma.temp_trop_VALUE  = 0.3*0.01;       %% sigsqr
    driver.oem.sigma.temp_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.temp_trop_TOPLAY = driver.oem.sigma.temp_strat_BOTLAY + 1; %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% strat GAS 1
    driver.oem.sigma.hum_strat_VALUE  = 1.5*0.01;       %% sigsqr
    driver.oem.sigma.hum_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.hum_strat_BOTLAY = 49;  %% stop layer
  %% trop GAS 1
    driver.oem.sigma.hum_trop_VALUE  = 1.0*0.01;       %% sigsqr
    driver.oem.sigma.hum_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.hum_trop_TOPLAY = driver.oem.sigma.hum_strat_BOTLAY + 1; %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% length_correlation for T(z) and WV(z) .. note this is in terms of INDEX
    driver.oem.sigma.l_c = 2.4;
  %% QST 1 ..6 values
    driver.oem.sigma.qst(1) = 4*2.20;  %% co2
    driver.oem.sigma.qst(2) = 1*0.01;  %% o3
    driver.oem.sigma.qst(3) = 4*1.00;  %% n2o
    driver.oem.sigma.qst(4) = 1*5.00;  %% ch4 
    driver.oem.sigma.qst(5) = 1*1.00;  %% cfc11
    driver.oem.sigma.qst(6) = 1*0.10;  %% Stemp

fprintf(1,'%s \n',driver.rateset.datafile);

