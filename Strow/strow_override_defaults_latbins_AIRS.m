function driver = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
% Which latitude bin
ix = driver.iibin;
regress_rates = driver.rateset.unc_rates;
%---------------------------------------------------------------------------
% Load in freq corrections
load Data/dbt_10year  % alldbt
driver.rateset.rates = driver.rateset.rates-alldbt(ix,:)'/10;

% Modify rates with lag-1 correlation errors or add to above
nc_cor = nc_rates(driver);
% Modify with estimated error in freq + regress errors 
driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates.*nc_cor;
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
l_c = 0.4;
mat_od = exp(-mat_od.^2./(1*l_c^2));

% Override a priori values 

xb(1)=2.2; 

% Attempt to change l_c in strat, unstable as Andy said
% Maybe should change Jacobians and have fewer strat layers
% l_c_trop = 2.4;
% l_c_strat = 2.4;
% mat_od_strat = mat_od(1:trpi(ix),1:trpi(ix));
% mat_od_trop  = mat_od(trpi(ix)+1:end,trpi(ix)+1:end);
% mat_od_strat = exp(-mat_od_strat.^2./(1*l_c_strat^2));
% mat_od_trop = exp(-mat_od_trop.^2./(1*l_c_trop^2));
% No off-diagonal between strat and trop
% mat_od = blkdiag(mat_od_strat,mat_od_trop);

for i=1:36
   ct(i).trans1 = trpi(i);
   ct(i).trans2 = trpi(i)+10;
   ct(i).lev1 = 0.02;
   ct(i).lev2 = 0.01;
   ct(i).lev3 = ct(i).lev2;
   ct(i).width1 = 1/5;
   ct(i).width2 = ct(i).width1;

   cw(i).trans1 = trpi(i);
   cw(i).trans2 = trpi(i)+10;
   cw(i).lev1 = 0.01;
   cw(i).lev2 = 0.005;
   cw(i).lev3 = cw(i).lev2;
   cw(i).width1 = 1/5;
   cw(i).width2 = cw(i).width1;
end

% Temperature level uncertainties, then scaled and squared
%tunc     = ones(1,pmat_size)*0.01;
tunc = cov2lev(ct(ix));
t_sigma = (tunc./tnorm).^2;
% Make cov matrix
tmat = (t_sigma'*t_sigma).*mat_od;
driver.oem.tunc = tunc;

% Water level uncertainties, then scaled and squared
%wunc     = ones(1,pmat_size)*0.02;
wunc = cov2lev(cw(ix));
w_sigma = (wunc./wnorm).^2;
% Make cov matrix
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) O3(frac) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    
%fmat_orgi = [5/2.2     0.02       2      0.2      0.8        0.01];
fmatd = [2     0.1       2      10      1        0.1];
fmat  = diag(fmatd.*fnorm); 
driver.oem.cov = blkdiag(fmat,wmat,tmat);
