function [driver,aux] = strow_override_defaults_latbins_AIRS(driver);

%%%%%%%%%% look at /home/sergio/MATLABCODE/oem_pkg_run_Aug2014/Strow/overrides_history.m
%%%%%%%%%% look at /home/sergio/MATLABCODE/oem_pkg_run_Aug2014/Strow/overrides_history.m
%%%%%%%%%% look at /home/sergio/MATLABCODE/oem_pkg_run_Aug2014/Strow/overrides_history.m
%%%%%%%%%% look at /home/sergio/MATLABCODE/oem_pkg_run_Aug2014/Strow/overrides_history.m
%%%%%%%%%% look at /home/sergio/MATLABCODE/oem_pkg_run_Aug2014/Strow/overrides_history.m

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

%% see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS
%% see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS/SARTA_CLOUD_JACS
driver.jacobian.filename = '/asl/s1/rates/Cloud/Apr2014/sarta_M_TS_jac_all_6_4_97_97_cld_replaceT.mat';

driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:10;
driver.jacobian.water_i  = 11:107;
driver.jacobian.temp_i   = 108:204;
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
xb = zeros(1,204);

xb(1) = 0;  % Set CO2 apriori to zero (now at 2)
%xb(1) = 2.2;  % Set CO2 apriori to zero (now at 2)

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
driver.rateset.unc_rates = ones(2378,1)*0.001;
%driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates;%.*nc_cor;
%driver.rateset.unc_rates = driver.rateset.unc_rates;%.*nc_cor;

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

ct(ix).stren = 0.005;
ct(ix).trans1 = trpi(ix);
%ct(ix).trans2 = trpi(ix)+10;
ct(ix).lev1 = ct(ix).stren;
ct(ix).lev2 = ct(ix).stren/2;
ct(ix).lev3 = ct(ix).stren/4;
ct(ix).width1 = 1/5;
%ct(ix).width2 = ct(ix).width1;

cw(ix).stren = 0.001;
cw(ix).trans1 = trpi(ix);
%cw(ix).trans2 = trpi(ix)+10;
cw(ix).lev1 = cw(ix).stren;
cw(ix).lev2 = cw(ix).stren/2;
cw(ix).lev3 = cw(ix).stren/4;
cw(ix).width1 = 1/5;
%cw(ix).width2 = cw(ix).width1;

% Temperature level uncertainties, then scaled and squared
tunc = cov2lev(ct(ix));
t_sigma = (tunc./tnorm);
% Make cov matrix
tmat = (t_sigma'*t_sigma).*mat_od;
driver.oem.tunc = tunc;

% Water level uncertainties, then scaled and squared
wunc = cov2lev(cw(ix));
w_sigma = (wunc./wnorm);
% Make cov matrix
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) O3(frac) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    
%fmat_orgi = [5/2.2     0.02       2      0.2      0.8        0.01];
fmatd = [2     0.1       2      1      1        0.1];

%% put in clouds
% 'CNG1'    'CNG2'    'CTP1'    'CTP2' 
fmatdCLD = [0.01 0.01 0.01 0.01];
fmatd = [fmatd fmatdCLD];

fmatd = [2 0.1 1 10 1 0.1 1/20 1/20 1/20 1/20]*1;

fmatd = [0.00001 0.1 1 10 1 0.1 1/20 1/20 1/20 1/20]*1;
fmatd = [0.00001 0.1 1 10 1 0.1 1/20 1/20 1/20 1/20]*0.1;
fmatd = [2       0.1 1 10 1 0.1 1/20 1/20 1/20 1/20]*1;

fmat  = diag(fmatd./fnorm); 
driver.oem.cov = blkdiag(fmat,wmat,tmat);

%---------------------------------------------------------------------
% Empirical regularization parameters and switches
driver.oem.reg_type = 'reg_and_cov'; % 'reg_and_cov','cov','reg' are other choices
% Separate reg weights for water, temperature profiles
driver.oem.alpha_water = 5*10;
driver.oem.alpha_temp = 1*10;

