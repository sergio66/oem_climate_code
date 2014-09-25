function [driver,aux] = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
% Which latitude bin
ix = driver.iibin;
%---------------------------------------------------------------------------
% Fitting [obs][cal][bias], pick one
driver.rateset.ocb_set  = 'obs';
%---------------------------------------------------------------------------
% Raw rate data file
%driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
driver.rateset.datafile  = '/asl/s1/rates/Clear/lls_robust_all_fitout_gdays.mat';
driver.rateset.datafile = '/asl/s1/rates/Cloud/Apr2014/xoverocean_gsx_1day_ctr2378_era_cld_5deg_lays_spanday01_avgL1Brates_robust_Apr16_2014_span_09_2002_08_2012.mat';

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../../oem_pkg/Test/all_lagcor.mat';
% Get rate data, do Q/A elsewhere
driver = get_rates(driver);
%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
driver.jacobian.filename = '/asl/s1/rates/Cloud/Apr2014/sarta_M_TS_jac_all_6_6_97_97_cld.mat';
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:12;
driver.jacobian.water_i  = 13:109;
driver.jacobian.temp_i   = 110:206;
driver.jacobian.numlays  = 97;

% Get jacobians
jac             = load(driver.jacobian.filename);
aux.m_ts_jac    = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Good channel set
load /asl/s1/rates/Clear/good_chanset.mat 
driver.jacobian.chanset = chanset;

% Remove ozone channels
%k = find(f(chanset) < 970 |  f(chanset) > 1117 );
k = find(f(chanset) < 970 |  f(chanset) > 1200 );
driver.jacobian.chanset = chanset(k);
%---------------------------------------------------------------------------
% Apriori file
driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;
xb = zeros(1,206);
%xb(1) = 0;  % Set CO2 apriori to zero (now at 2)
[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.iibin);
end
% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';
%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
% Convert radiance rates into bt rate
% deriv = drdbt(f,rad2bt(f,driver.rateset.r));
% driver.rateset.rates = driver.rateset.rates./(1000*deriv);
% driver.rateset.unc_rates = driver.rateset.unc_rates./(1000*deriv);
%---------------------------------------------------------------------------
% Load in freq corrections
load Data/dbt_10year  % alldbt
driver.rateset.rates = driver.rateset.rates-alldbt(ix,:)'/10;

% % Offset rates to get sensitivity to AIRS BT drift
% if driver.rateset.ocb_set  == 'obs';
%    driver.rateset.rates = driver.rateset.rates + 0.1;
% end
%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors or add to above
%nc_cor = nc_rates(driver);
% Modify with estimated error in freq + regress errors 
%driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates.*nc_cor;
driver.rateset.unc_rates = ones(2378,1)*0.001;
%---------------------------------------------------------------------------
% Do rate Q/A (empty for now)
%---------------------------------------------------------------------------
% Build covariance matrices; fixed, water, temperature
pmat_size = driver.jacobian.numlays;

% Normalization depends on parameter
fnorm = driver.qrenorm(driver.jacobian.scalar_i);
wnorm = driver.qrenorm(driver.jacobian.water_i);
tnorm = driver.qrenorm(driver.jacobian.temp_i);

% Make sure re-scale to Jacobians before squaring cov matrix
for i=1:pmat_size
   for j=1:pmat_size
      mat_od(i,j) = abs(i-j);
   end
end

% Defines tropopause index trpi
trop_index

% Relative off-diagonal
% l_c = 2.4;
l_c = 1.25;
mat_od = exp(-mat_od.^2./(1*l_c^2));


%%% clear sky
%{
for i=1:36
   ct(i).trans1 = trpi(i);
   ct(i).trans2 = trpi(i)+10;
   ct(i).lev1 = 0.005*4;
   ct(i).lev2 = 0.005*4;
%    ct(i).lev1 = 0.002;
%    ct(i).lev2 = 0.002;
   ct(i).lev3 = ct(i).lev2;
   ct(i).width1 = 1/5;
   ct(i).width2 = ct(i).width1;

   cw(i).trans1 = trpi(i);
   cw(i).trans2 = trpi(i)+10;
%    cw(i).lev1 = 0.01;
%    cw(i).lev2 = 0.005;
   cw(i).lev1 = 0.005*4;
   cw(i).lev2 = 0.005*4;
   cw(i).lev3 = cw(i).lev2;
   cw(i).width1 = 1/5;
   cw(i).width2 = cw(i).width1;
end
%}

for i=1:36
   ct(i).trans1 = trpi(i);
   ct(i).trans2 = trpi(i)+10;
   ct(i).lev1 = 0.02*10*1;
   ct(i).lev2 = 0.01*10*2;    %% for clouds
   ct(i).lev3 = ct(i).lev2;
   ct(i).width1 = 1/5*1/2;    %% for clouds
   ct(i).width1 = 1/5;        %% for clouds
   ct(i).width2 = ct(i).width1;

   cw(i).trans1 = trpi(i);
   cw(i).trans2 = trpi(i)+10;
   cw(i).lev1 = 0.01*10*1;
   cw(i).lev2 = 0.005*10*2;   %% for clouds
   cw(i).lev2 = 0.005;
   cw(i).lev3 = cw(i).lev2;
   cw(i).width1 = 1/5*1/2;    %% for clouds
   cw(i).width1 = 1/5;
   cw(i).width2 = cw(i).width1;
end

% Temperature level uncertainties, then scaled and squared
%tunc     = ones(1,pmat_size)*0.01;
tunc = cov2lev(ct(ix));
t_sigma = (tunc./tnorm);%.^2;
% Make cov matrix
tmat = (t_sigma'*t_sigma).*mat_od;
driver.oem.tunc = tunc;

% Water level uncertainties, then scaled and squared
%wunc     = ones(1,pmat_size)*0.02;
wunc = cov2lev(cw(ix));
w_sigma = (wunc./wnorm);%.^2;
% Make cov matrix
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) O3(frac) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    
%fmat_orgi = [5/2.2     0.02       2      0.2      0.8        0.01];
fmatd = [2     0.1       2      1      1        0.1];

%% put in clouds
% 'CNG1'    'CNG2'    'CTP1'    'CTP2'    'CSZ1'    'CSZ2' 
fmatdCLD = [0.01 0.01 0.01 0.01 0.01 0.01];

fmatd = [fmatd fmatdCLD];
fmat  = diag(fmatd./fnorm); 
driver.oem.cov = blkdiag(fmat,wmat,tmat);

%---------------------------------------------------------------------
% Empirical regularization parameters and switches
driver.oem.reg_type = 'reg_and_cov'; % 'reg_and_cov','cov','reg' are other choices
% Separate reg weights for water, temperature profiles
driver.oem.alpha_water = 200;
driver.oem.alpha_temp = 10;


