function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
%% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +1;
driver.rateset.min = -1;

driver.oem.apriori_filename = 'apriori_zero';
y = instr_chans('iasi');
ind = find(y <= 2205);
ind = (1:4:length(ind));
ind = ind(ind > 200);

ind = find((y >= 725 & y <= 800) | (y >= 930 & y <= 980) | ...
           (y >= 1250 & y <= 1500) | (y >= 2100 & y <= 2250)); 
ind = ind(1:2:length(ind));
ind = find((y >= 725 & y <= 800) | (y >= 930 & y <= 980) | ...
           (y >= 1250 & y <= 1650) | (y >= 2100 & y <= 2250)); 
ind = ind(1:2:length(ind)); indA = ind;

driver.jacobian.chanset = (1:4230);
driver.jacobian.chanset = (1:4:8200);
dfs = '/strowdata1/shared/sergio/MATLABCODE/DFS_OPTIMUMCHANS_RODGERS/';
c1 = load([dfs 'iasi_optchannelZ.mat']); %% o3
c2 = load([dfs 'iasi_optchannelM.mat']); %% ch4
c3 = load([dfs 'iasi_optchannelW.mat']); %% h2o
c4 = load([dfs 'iasi_optchannelT.mat']); %% t
c5 = unique([c1.iaFreqX; c2.iaFreqX; c3.iaFreqX; c4.iaFreqX]);
driver.jacobian.chanset = c5(c5 > 422);  %% so we start at 750 cm-1
%driver.jacobian.chanset = c5(c5 > 262);  %% so we start at 710 cm-1
%blonk = find(y < 980 | y > 1100);
%driver.jacobian.chanset = intersect(driver.jacobian.chanset,blonk);
driver.jacobian.chanset = [(300:5:4230)  582:592];
driver.jacobian.chanset = [[282:592] [1342:1820] [2642:3040]];
driver.jacobian.chanset = [[282:840] [1342:1820] [2662:3240]];
driver.jacobian.chanset = [[282:840] [1342:1820] [2602:3240]];
driver.jacobian.chanset = [[282:780] [1342:1740] [2602:3240]];
%driver.jacobian.chanset = c4.iaFreqX;
%  driver.jacobian.chanset = c4.iaFreqX(c4.iaFreqX > 422);  %%start at 750 cm-1
%  driver.jacobian.chanset = c4.iaFreqX(c4.iaFreqX > 262);  %%start at 710 cm-1
driver.jacobian.chanset = sort(unique(driver.jacobian.chanset));

indB = find(y <= 2205);
indC = find(y <= 2205 | y >= 2405);

driver.jacobian.chanset = indA;
driver.jacobian.chanset = indB;
driver.jacobian.chanset = indC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% larrabee used these in the file he gave me, but I get bad trace gas rates
driver.oem.cov_filename     = 'cov_lls.mat'
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;

%%% these give decent CO2, N2O rates!!!!! (AIRS)
driver.oem.lambda           = 0.01 * 1;
driver.oem.lambda_qst       = 0.01 * 1;
driver.oem.lambda_Q1        = 0.01 * 1;
driver.oem.lambda_temp      = 0.01 * 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 50;
driver.oem.lambda_temp      = 50;
driver.oem.lambda_Q1        = 0.1;
driver.oem.lambda_temp      = 0.1;

%% used for LLS in Aug 2012, for Utah meeting
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_Q1        = 50;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_Q1        = 100;
driver.oem.lambda_temp      = 50;

driver.rateset.datafile =  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/avg_strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_iasi_lays_spanday01_avgL1Brates_robust_Sep11_2012_iasiB1_span_07_2007_07_2012.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% orig 6_97_97 stuff
driver.jacobian.filename = ...
  '/home/sergio/MATLABCODE/BUFFER_Fit_pkg/Aux_jacs_IASI_General/JACS/Sep11_2012/all_kcarta_6_97_97_jacsB1.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1       1       1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays    = 97;
driver.oem.lambda_qst       = [0.1 0.2 0.3 0.4 0.5 0.6]*0.10*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*2 ones(1,57)*3]*100;
driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*5 ones(1,57)*1]*100;

driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1     = 0.1;
driver.oem.lambda_temp      = 0.1;

driver.oem.lambda_Q1     = 10;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_Q1     = 50;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_Q1     = 100;
driver.oem.lambda_temp      = 50;

driver.jacobian.chanset = indB;

% can individually tweak
driver.oem.lambda_qst       = [10 10 10 10 10 10                              ones(1,1)*(-9999)]*15;
driver.oem.lambda_Q1        = [ones(1,10)*10.0 ones(1,30)*5.0 ones(1,57)*0.05 ones(1,1)*(-9999)]*1000;
driver.oem.lambda_temp      = [ones(1,10)*10.0 ones(1,30)*5.0 ones(1,57)*0.05 ones(1,1)*(-9999)]*500;

%% yuk
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 5;
driver.oem.lambda_temp      = 5;
driver.oem.lambda_qst       = [0.1 0.1 0.1 0.1 0.10 0.1]*1e-1;
driver.oem.lambda           = 0;
driver.oem.adjust_spectral_errorbars = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%
