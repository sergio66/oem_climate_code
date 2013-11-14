function driver = override_defaults(driver,ix);


% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

% good!!!!%%%
driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';


%driver.rateset.datafile = '/asl/s1/rates/clear/Oct2013_MERRA/xoverocean__lays_spanday01_avgL1Brates_robust_Oct5_2013_span_09_2002_08_2012.mat';
% good!!!!%%%

%load /asl/s1/rates/clear/strow_stmNov2011_dobs20.mat
load /asl/s1/rates/clear/good_chanset.mat 
driver.jacobian.chanset = chanset;
%driver.jacobian.chanset = dobs20.jacobian.chanset;
% Remove some channels 1864 
%k=find(driver.jacobian.chanset >=500)
%driver.jacobian.chanset=driver.jacobian.chanset(k); 

% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +0.25;
driver.rateset.min = -0.25;

% orig 6_97_97 stuff
driver.rateset.ncfile    = 'all_lagcor.mat';
driver.rateset.ocb_set = 'obs'; 
%driver.rateset.ocb_set  = 'bias';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   % how many gases we want to retrieve profiles
                                  % must be at least 1 (for water)
driver.jacobian.numlays     = 97;

aux_stuff.xb(1)=0.0; 
%########################################################################
% larrabee used these in the file he gave me, but I get bad trace gas rates

% they give decent WV/T rates

% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV
driver.oem.lambda_qst       = [1 1 1 1 1 1                                 ones(1,1)*(-9999)]*1e-4;
driver.oem.lambda_Q1        = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*50;
driver.oem.lambda_Q1        = [ones(1,40)*10 ones(1,30)*15 ones(1,27)*20 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*0.9 ones(1,27)*0.8 ones(1,1)*(-9999)]*5;
driver.oem.lambda_temp      = [ones(1,40)*1.0 ones(1,30)*1.7 ones(1,27)*2.5 ones(1,1)*(-9999)]*5;
driver.oem.lambda           = 1e-2;
driver.jacobian.filename = 'M_TS_jac_all.mat';
% WORKS WELL WITH Nov02_2012_span_09_2002_08_2012.mat gives good T WV

%driver.rateset.adjust = 1;
%driver.ratest.adjust_index = 1; 
%driver.rateset.adjust_values = 2.2; 


%driver.jacobian.filename = '/asl/s1/rates/clear/Aug2013/all_kcarta_HARTMANNCO2_jacs_NOV02_2012_HARTMANN_6_97_97.mat';
%driver.jacobian.filename = 'M_TS_jac_all.mat';
%driver.jacobian.filename = '/asl/s1/rates/clear/Aug2013/all_kcarta_UMBCCO2_jacs_NOV02_2012_UMBCCO2_6_97_97.mat';

smoothVSregularization='c'; 
  %% but if smoothVSregularization = 'c'
  %% note these are PHYSICAL units eg sigma_temp_stratVALUE is in KELVIN (or KELVIN/YR)
  %% so eg Andy Tangborn finds after normalization driver.oem.sigma.temp_strat_VALUE  = 4;
  %% is a "good value" which means we default set it to 4*0.01 (0.01 is qrenorm for T(z))

  %% Latitude dependent standard deviations. 
     % Stratosphere temperature
     TS(1:36)=2.0;  TS(36)=3.0;  TS(35)=3.0; TS(34)=3.0; TS(33)=3.0; TS(32)=2.0; TS(31)=2.0; TS(30)=2.0; TS(29)=2.0; TS(28)=2.0; TS(27)=2.0; TS(26)=2.5; TS(25)=2.5; TS(24)=2.5; TS(23)=2.5; TS(22)=2.5; TS(21)=2.5; TS(20)=2.5; TS(19)=2.5; TS(18)=3.5; TS(17)=3.5; TS(16)=3.5; TS(15)=3.5; TS(14)=3.5; TS(13)=3.5; TS(12)=3.5; TS(11)=2.5; TS(10)=2.5; TS(9)=2.5; TS(8)=2.5; TS(7)=2.5; TS(6)=2.1; TS(5)=1.5; TS(4)=1.5; TS(3)=1.2; TS(2)=1.0; TS(1)=2.0;              
     % Upper troposphere temperature
     TUT(1:36)=2.0; TUT(36)=3.0; TUT(35)=3.0; TUT(34)=3.0; TUT(33)=1.5; TUT(32)=0.8; TUT(31)=0.8; TUT(30)=0.8; TUT(29)=0.8; TUT(28)=0.5; TUT(27)=0.4; TUT(26)=0.5; TUT(25)=0.5; TUT(24)=0.5; TUT(23)=0.5; TUT(22)=0.5; TUT(21)=0.5; TUT(20)=0.5; TUT(19)=0.5; TUT(18)=0.9; TUT(17)=0.9; TUT(16)=0.9; TUT(15)=0.9; TUT(14)=0.9; TUT(13)=0.9; TUT(12)=0.9; TUT(11)=0.4; TUT(10)=0.9; TUT(9)=1.8; TUT(8)=2.8; TUT(7)=2.8; TUT(6)=2.8; TUT(5)=2.8; TUT(4)=2.8; TUT(3)=2.8; TUT(2)=2.8; TUT(1)=0.0008;   
     % Lower troposphere temperature
     TLT(1:36)=2.0; TLT(36)=3.0; TLT(35)=5.0; TLT(34)=1.0; TLT(33)=0.8; TLT(32)=0.5; TLT(31)=0.3; TLT(30)=0.3; TLT(29)=0.3; TLT(28)=0.3; TLT(27)=0.2; TLT(26)=0.2; TLT(25)=0.12; TLT(24)=0.1; TLT(23)=0.1; TLT(22)=0.1; TLT(21)=0.1; TLT(20)=0.1; TLT(19)=0.1; TLT(18)=0.25; TLT(17)=0.25; TLT(16)=0.25; TLT(15)=0.25; TLT(14)=0.25;  TLT(13)=0.25; TLT(12)=0.25; TLT(11)=0.25; TLT(10)=0.25; TLT(9)=0.3; TLT(8)=0.5; TLT(7)=0.5; TLT(6)=0.5; TLT(5)=0.5; TLT(4)=0.5; TLT(3)=1.5; TLT(2)=0.8; TLT(1)=1.8;  
     % Stratosphere water vapor
     WVS(1:36)=1.0; WVS(36)=0.5; WVS(35)=0.5; WVS(34)=0.5; WVS(33)=1.0; WVS(32)=1.0; WVS(31)=2.0; WVS(30)=1.0; WVS(29)=1.0; WVS(28)=0.6; WVS(27)=0.5; WVS(26)=0.4; WVS(25)=0.3; WVS(24)=0.3; WVS(23)=0.3; WVS(22)=0.3; WVS(21)=0.3; WVS(20)=0.3; WVS(19)=0.3; WVS(18)=0.3; WVS(17)=0.3; WVS(16)=0.3; WVS(15)=0.3; WVS(14)=0.3; WVS(13)=0.3; WVS(12)=0.3; WVS(11)=0.3; WVS(10)=0.3; WVS(9)=0.3; WVS(8)=0.3; WVS(7)=0.3; WVS(6)=0.25; WVS(5)=0.25; WVS(4)=0.25; WVS(3)=.25; WVS(2)=0.25; WVS(1)=0.25;   
     % Troposphere water Vapor
     WVT(1:36)=1.0; WVT(36)=0.3; WVT(35)=0.3; WVT(34)=0.3; WVT(33)=0.25; WVT(32)=0.1; WVT(31)=0.1; WVT(30)=0.2; WVT(29)=0.2; WVT(28)=0.3; WVT(27)=0.3; WVT(26)=0.2; WVT(25)=0.15; WVT(24)=0.18; WVT(23)=0.18; WVT(22)=0.18; WVT(21)=0.18; WVT(20)=0.18; WVT(19)=0.25; WVT(18)=0.12; WVT(17)=0.08; WVT(16)=0.05; WVT(15)=0.05; WVT(14)=0.05; WVT(13)=0.05; WVT(12)=0.04; WVT(11)=0.04; WVT(10)=0.04; WVT(9)=0.20; WVT(8)=0.2; WVT(7)=0.14; WVT(6)=0.14; WVT(5)=0.14; WVT(4)=0.14; WVT(3)=0.14; WVT(2)=0.08; WVT(1)=0.20;     
     % Trace Gases
     CO2(1:36)=1.0; CO2(36)=1.0; CO2(35)=1.0; CO2(34)=1.0; CO2(33)=1.0; CO2(32)=5.0; CO2(31)=5.0; CO2(30)=5.0; CO2(29)=5.0; CO2(28)=5.0; CO2(27)=5.0; CO2(26)=5.0; CO2(25)=5.0; CO2(24)=5.0; CO2(23)=5.0; CO2(22)=5.0; CO2(21)=5.0; CO2(20)=5.0; CO2(19)=5.0; CO2(18)=5.0; CO2(17)=5.0; CO2(16)=5.0; CO2(15)=5.0; CO2(14)=5.0; CO2(13)=5.0; CO2(12)=5.0; CO2(11)=5.0; CO2(10)=5.0; CO2(9)=5.0; CO2(8)=5.0; CO2(7)=5.0; CO2(6)=5.0; CO2(5)=5.0; CO2(4)=5.0; CO2(3)=5.0; CO2(2)=5.0; CO2(1)=1.0;  
     O3(1:36)=1.0;  O3(36)=1.0;  O3(35)=1.0;  O3(34)=1.0; O3(33)=1.0; O3(32)=2.0; O3(31)=2.0; O3(30)=2.0; O3(29)=2.0; O3(28)=2.0; O3(27)=2.0; O3(26)=2.0;O3(25)=2.0;O3(24)=2.0; O3(23)=2.0; O3(22)=2.0; O3(21)=2.0; O3(20)=2.0; O3(19)=2.0; O3(18)=2.0; O3(17)=2.0; O3(16)=2.0; O3(15)=2.0; O3(14)=2.0; O3(13)=2.0; O3(12)=2.0; O3(11)=2.0; O3(10)=2.0; O3(9)=2.0; O3(8)=2.0; O3(7)=2.0; O3(6)=2.0; O3(5)=2.0; O3(4)=2.0; O3(3)=2.0; O3(2)=2.0; O3(1)=2.0;  
     N2O(1:36)=1.0; N2O(36)=1.0; N2O(35)=1.0; N2O(34)=1.0; N2O(33)=1.0; N20(32)=2.0; N20(31)=2.0; N2O(30)=2.0; N2O(29)=2.0; N2O(28)=2.0; N2O(27)=2.0; N2O(26)=2.0; N2O(25)=2.0;N2O(24)=2.0; N2O(23)=2.0; N2O(22)=2.0; N2O(21)=2.0;  N2O(20)=2.0; N2O(19)=2.0; N20(18)=2.0; N2O(17)=2.0; N2O(16)=2.0; N2O(15)=2.0; N2O(14)=2.0; N2O(13)=2.0; N2O(12)=2.0; N2O(11)=2.0; N2O(10)=2.0; N2O(9)=2.0; N2O(8)=2.0; N2O(7)=2.0; N2O(6)=2.0; N2O(5)=2.0; N2O(4)=2.0; N2O(3)=2.0; N2O(2)=2.0; N2O(1)=2.0;    
     CH4(1:36)=1.0; CH5(36)=1.0; CH5(35)=1.0; CH5(34)=1.0; CH5(33)=1.0; CH5(32)=2.0; CH5(31)=1.0; CH5(30)=1.0; CH5(29)=1.0; CH5(28)=1.0; CH5(27)=1.0; CH5(26)=1.0; CH5(25)=1.0; CH5(24)=1.0; CH5(23)=1.0; CH5(22)=1.0; CH5(21)=1.0; CH5(20)=1.0; CH5(19)=1.0; CH5(18)=1.0; CH5(17)=1.0; CH5(16)=1.0; CH5(15)=1.0; CH5(14)=1.0; CH5(13)=1.0; CH5(12)=1.0; CH5(11)=1.0; CH5(10)=1.0; CH5(9)=1.0; CH5(8)=1.0; CH5(7)=1.0; CH5(6)=1.0; CH5(5)=1.0; CH5(4)=1.0; CH5(3)=1.0; CH5(2)=1.0; CH5(1)=1.0;  
     CFC11(1:36)=1.0; 
     % Surface temperature  
     STEMP(1:36)=1.0; STEMP(36)=20.0; STEMP(35)=0.1; STEMP(34)=0.1; STEMP(33)=0.1; STEMP(32) = 0.1; STEMP(31)=0.1; STEMP(30)=0.1; STEMP(29)=3.0; STEMP(28)=0.1; STEMP(27)=0.1; STEMP(26)=0.1; STEMP(25)=0.1; STEMP(24)=0.1; STEMP(23)=0.1; STEMP(22)=0.1; STEMP(21)=3.0;  STEMP(20)=3.0; STEMP(19)=3.0; STEMP(18)=3.0; STEMP(17)=3.0; STEMP(16)=0.2; STEMP(15)=0.2; STEMP(14)=0.2; STEMP(13)=0.2; STEMP(12)=0.2; STEMP(11)=0.2; STEMP(10)=0.2; STEMP(9)=0.2; STEMP(8)=0.2; STEMP(7)=0.2; STEMP(6)=0.2; STEMP(5)=0.2; STEMP(4)=0.2; STEMP(3)=0.2; STEMP(2)=0.2; STEMP(1)=1.2;  
     % Correlation length (in grid units) 
     C_L(1:36) = 1.5; C_L(36)=2.1; C_L(35)=1.8; C_L(34)=1.8; C_L(33)=1.8; C_L(32) =2.2; C_L(31)=2.2; C_L(30)=2.2;C_L(29)=2.2; C_L(28)=2.2; C_L(27)=2.2;   C_L(26)=2.2; C_L(25)=2.2; C_L(24)=2.2; C_L(23)=2.2; C_L(22)=2.2; C_L(21)=2.2; C_L(20)=2.2; C_L(19)=2.2; C_L(18)=2.2;  C_L(17)=2.2; C_L(16)=2.2; C_L(15)=2.2; C_L(14)=2.2; C_L(13)=2.2; C_L(12)=2.2; C_L(11)=2.2;  C_L(10)=2.2; C_L(9)=2.2; C_L(8)=2.2; C_L(7)=2.2; C_L(6)=2.2; C_L(5)=2.2; C_L(4)=2.2; C_L(3)=2.2; C_L(2)=2.2; C_L(1)=2.2;  
    
  %% strat temp
    driver.oem.sigma.temp_strat_VALUE  = TS(ix)*0.01;       %% sigsqr
    driver.oem.sigma.temp_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.temp_strat_BOTLAY = 49;  %% stop layer
  %% upper trop temp 
    driver.oem.sigma.temp_upper_trop_VALUE = TUT(ix)*0.01; 
    driver.oem.sigma.temp_upper_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.temp_upper_trop_BOTLAY = 79; 
  %% lower trop temp
    driver.oem.sigma.temp_trop_VALUE  = TLT(ix)*0.01;       %% sigsqr
    driver.oem.sigma.temp_trop_TOPLAY = 80;  %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.temp_trop_TOPLAY = driver.oem.sigma.temp_upper_trop_BOTLAY + 1; %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% strat GAS 1
    driver.oem.sigma.hum_strat_VALUE  = WVS(ix)*0.01;       %% sigsqr
    driver.oem.sigma.hum_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.hum_strat_BOTLAY = 49;  %% stop layer
  %% trop GAS 1
    driver.oem.sigma.hum_trop_VALUE  = WVT(ix)*0.01;       %% sigsqr
    driver.oem.sigma.hum_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.hum_trop_TOPLAY = driver.oem.sigma.hum_strat_BOTLAY + 1; %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% length_correlation for T(z) and WV(z) .. note this is in terms of INDEX
    driver.oem.sigma.l_c = C_L(ix);
  %% QST 1 ..6 values
    driver.oem.sigma.qst(1) = CO2(ix)*2.20;  %% co2
    driver.oem.sigma.qst(2) = O3(ix)*0.01;  %% o3
    driver.oem.sigma.qst(3) = N2O(ix)*1.00;  %% n2o
    driver.oem.sigma.qst(4) = CH4(ix)*5.00;  %% ch4 
    driver.oem.sigma.qst(5) = CFC11(ix)*1.00;  %% cfc11
    driver.oem.sigma.qst(6) = STEMP(ix)*0.10;  %% Stemp


%driver.oem.regularizationVScovariances = 'ERA'; 
driver.oem.regularizationVScovariances = 'C'; 

driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 0.1;
driver.oem.lambda_temp      = 10;
fprintf(1,'%s \n',driver.rateset.datafile);

