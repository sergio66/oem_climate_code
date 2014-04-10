function driver = strow_override_defaults_latbins_AIRS(driver,ix);

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

% 10-year rate file
driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
% 10-year Merra calc rates?
%driver.rateset.datafile = '/asl/s1/rates/clear/Oct2013_MERRA/xoverocean__lays_spanday01_avgL1Brates_robust_Oct5_2013_span_09_2002_08_2012.mat';

% Good channel set
load /asl/s1/rates/clear/good_chanset.mat 
driver.jacobian.chanset = chanset;

% Remove channels outside the max/min B(T) values
driver.rateset.max = 320;
driver.rateset.min = 180;
% Remove channels outside the max/min dB(T)/dt values
driver.rateset.max = +0.25;
driver.rateset.min = -0.25;

% Rate errors
driver.rateset.ncfile    = 'all_lagcor.mat';
% Fits observed rates
driver.rateset.ocb_set = 'obs'; 
% Fit bias rates (vs ERA, etc.)
%driver.rateset.ocb_set  = 'bias';

% Jacobian definitions
driver.jacobian.qstnames   = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [  1    1    1     1     1       1];
% How many gas profiles to fit (at least 1for water vapor)
driver.jacobian.numQlays   = 1;
% Number of layers in Jacobians
driver.jacobian.numlays     = 97;
% Jacobian file, (why no full path?)
driver.jacobian.filename = 'M_TS_jac_all.mat';

% What is this?
aux_stuff.xb(1)=0.0; 

driver.oem.regularizationVScovariances = 'C'; 
% driver.oem.regularizationVScovariances = 'ERA'; 

% If want to use Andy's old way of setting up covariance
% andy_cov_build

% If want to do old empirical regularization instead (OLD)
% empirical_reg

% Build covariance matrices; fixed, water, temperature

% num_fixed = length(find(driver.jacobian.qstYesOrNo == 1));
% if driver.jacobian.numQlays ~= 1
%    disp('Error, only T and Q profiles together for now:')
% end

% Using pmat (= 97) here, will kinda assume hard-coded for now
pmat_size = driver.jacobian.numlays;
pmat_size = 97;

load qrenorm
fnorm = qrenorm(1:6);
wnorm = qrenorm(7:103);
tnorm = qrenorm(104:200);

% Make sure re-scale to Jacobians before squaring cov matrix
for i=1:pmat_size
   for j=1:pmat_size
      mat_od(i,j) = abs(i-j);
   end
end

% Relative off-diagonal
l_c = 2.4;
mat_od = exp(-mat_od.^2./(1*l_c^2));

% % Sample values for c structure, use co3lev.m to build
% %   now using normal physical units
% c.trans1 = 20;
% c.trans2 = 73;
% c.width1 = 1/2;
% c.width2 = 1/3;
% c.lev1 = 1;
% c.lev2 = 3;
% c.lev3 = 4;

% Ignore stuff above, for now use constant profile variances
% and conentrate on getting the off-diagonals correct

% Temperature level uncertainties, then scaled and squared
tunc     = ones(1,pmat_size)*0.02;
t_sigma = (tunc./tnorm).^2;
% Make cov matrix
tmat = (t_sigma'*t_sigma).*mat_od;

% Temperature level uncertainties, then scaled and squared
wunc     = ones(1,pmat_size)*0.005;
w_sigma = (wunc./wnorm).^2;
% Make cov matrix
wmat = (w_sigma'*w_sigma).*mat_od;

%fmat = zeros(6,6);
%            CO2(ppm) O3(frac) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    
%fmatd = [5/2.2     0.02       2      0.2      0.8        0.01];
%     CO2(ppm) O3(frac) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    
fmatd = [2     0.02       2      0.2      0.8        0.01];
fmat  = diag(fmatd.*fnorm); 

driver.oem.fmat = fmat;
driver.oem.wmat = wmat;
driver.oem.tmat = tmat;

