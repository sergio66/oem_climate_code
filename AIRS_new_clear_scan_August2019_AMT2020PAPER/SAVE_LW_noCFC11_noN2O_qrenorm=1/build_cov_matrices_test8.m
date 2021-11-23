% Build covariance matrices; fixed, water, temperature
pmat_size = driver.jacobian.numlays;

% Normalization depends on parameter
fnorm = driver.qrenorm(driver.jacobian.scalar_i);
wnorm = driver.qrenorm(driver.jacobian.water_i);
tnorm = driver.qrenorm(driver.jacobian.temp_i);
oznorm = driver.qrenorm(driver.jacobian.ozone_i);

% Make sure re-scale to Jacobians before squaring cov matrix
for i=1:pmat_size
  for j=1:pmat_size
    mat_od(i,j) = abs(i-j);
  end
end

% Defines tropopause index trpi
trop_index

%%%%%%%%%%%%%%%%%%%%%%%%%
%%         lc   ct.lev1  ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2 coz_wide  alpha_T  alpha_w  alpha_oz
%% THIS IS FOR ERA CALS 
cov_set = [2.4  0.002    0.002     1/5       0.005*4  0.005*4     1/5      0.005*4   0.005*4    1/5        10       200     200];
cov_set = [1.0  0.005*4  0.005*4   1/5       0.005*0.1  0.005*0.1 1/5      0.005*4   0.005*4    1/5        10       200     200]; %great era calc rates, O3 rates too1 large
cov_set = [1.0  0.005*4  0.005*4   1/5       0.005*0.1  0.005*0.1 1/5      0.001/2   0.001/1    1/2        10       200     20]; %great era calc rates, O3 getting better
cov_set = [1.0  0.005*4  0.005*4   1/5       0.005*0.1  0.005*0.1 1/5      0.001/2   0.001/0.75 1/2        10       200     20]; %great era calc rates, O3 getting even better

%% THIS IS FOR OBS and CALS, pretty good baseline start, convert_strowrates2oemrates_random_13_year_no_nucal_v2.mat
cov_set = [1.0  0.005*1  0.005*1   1/2       0.005*0.1  0.005*0.1 1/2      0.001/2   0.001/0.75 1/2        10       200     20]; 
cov_set = [1.0  0.005/2  0.005/2   1/2       0.005*0.1  0.005*0.1 1/2      0.001/2   0.001/0.75 1/2        10       20      20]; %pretty good for obs
cov_set = [1.0  0.005/2  0.005/2   1/2       0.005/20     0.005/20  1/2      0.001/2   0.001/0.75 1/2        10       2     20]; %pretty good for obs
cov_set = [1.0  0.005/2  0.005/2   1/2       0.005/25     0.005/25  1/2      0.001/2   0.001/0.75 1/2        10       2     20]; %pretty good for obs  YEAH YEAH YEAH

%Sergio above

% pretty nice
%cov_set = [1.0  0.005/2  0.005/2       1/2       0.005/4     0.005/4  1/2      0.01/2   0.01/2   1/2        5E4     8E7    2E6]; %pretty good for obs  YEAH YEAH YEAH

% good below for sept 2016
%cov_set = [1.0  0.005        0.005       1/2       0.005/4     0.005/4  1/2      0.01/4   0.01/4   1/2        1E6     8E7    2E5]; %pretty good for obs  YEAH YEAH YEAH

cov_set = [1.0  0.005        0.005       1/2       0.005/4     0.005/4  1/2      0.01/4   0.01/4   1/2        1E5     8E6    2E4]; %pretty good for obs  YEAH YEAH YEAH



%%%%%%%%%%%%%%%%%%%%%%%%%

% Relative off-diagonal
l_c = cov_set(1);;
mat_od = exp(-mat_od.^2./(1*l_c^2));

   ct(ix).trans1 = trpi(ix);
   ct(ix).trans2 = trpi(ix)+10;
   ct(ix).lev1 = cov_set(2);
   ct(ix).lev2 = cov_set(3);
   ct(ix).lev3 = ct(ix).lev2;
   ct(ix).width1 = cov_set(4);
   ct(ix).width2 = ct(ix).width1;

   cw(ix).trans1 = trpi(ix);
   cw(ix).trans2 = trpi(ix)+10;
   cw(ix).trans1 = cw(ix).trans1 + 10;  %% new move towards groud
   cw(ix).trans2 = cw(ix).trans2 + 10;  %% new move towards ground
   cw(ix).lev1 = cov_set(5);
   cw(ix).lev2 = cov_set(6);
   cw(ix).lev3 = cw(ix).lev2;
   cw(ix).width1 = cov_set(7);
   cw(ix).width2 = cw(ix).width1;

   coz(ix).trans1 = trpi(ix);
   coz(ix).trans2 = trpi(ix)+10;
   coz(ix).trans1 = coz(ix).trans1 - 40; %% new move towards TOA
   coz(ix).trans2 = coz(ix).trans2 - 40; %% new move towards TOA
   coz(ix).lev1 = cov_set(8);
   coz(ix).lev2 = cov_set(9);
   coz(ix).lev3 = coz(ix).lev2;
   coz(ix).width1 = cov_set(10);
   coz(ix).width2 = coz(ix).width1;

% Temperature level uncertainties, then scaled and squared
%tunc     = ones(1,pmat_size)*0.01;
tunc = cov2lev(ct(ix));
t_sigma = (tunc./tnorm);%.^2;
tmat = (t_sigma'*t_sigma).*mat_od;
driver.oem.tunc = tunc;

% Water level uncertainties, then scaled and squared
%wunc     = ones(1,pmat_size)*0.02;
wunc = cov2lev(cw(ix));
w_sigma = (wunc./wnorm);%.^2;
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

% ozone level uncertainties, then scaled and squared
%ozunc     = ones(1,pmat_size)*0.02;
ozunc = cov2lev(coz(ix));
oz_sigma = (ozunc./oznorm);%.^2;
ozmat = (oz_sigma'*oz_sigma).*mat_od;
driver.oem.ozunc = ozunc;

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) N2O(ppb) CH4(ppb) CFC11(ppt) Tsurf(K)    cng1 cng2 cpsz1 cpsz2]   
%fmat_orgi = [5/2.2       2      0.2      0.8        0.01      0.01  0.01  0.0001  0.0001];
fmatd      = [2         0.1        1      1          0.1       0.01  0.01  0.0001  0.0001]; %% works quite well for ERA CALCS

fmatd      = [2         0.1        1      1          0.1       0.01  0.01  0.01  0.01] * 2; %%pretty good!!! or OBS
                                                                                            %
% fmatd      = [1E-8         0.1        1      1          0.01       0.01  0.01  0.01  0.01] * 4; %%


% Goog b elow
%fmatd      = [0.3         20        20      1          0.01       0.01  0.01  0.01  0.01] * 4; %%


% For calcs
fmatd      = [1E-8      20     20      1          0.01       0.01  0.01  0.01  0.01] * 4; %%
 %Sergio'last above

%fmatd      = [2         0.1       1E7     1          0.1       0.01  0.01  0.01  0.01] * 4; %%

%fmatd      = [1e-8         0.1        1      1          0.1       0.01  0.01  0.01  0.01] * 4; %% 

%if length(driver.jacobian.scalar_i) ~ length(fmatd)
%  fmatd = fmatd(driver.jacobian.scalar_i);
%end

fmat  = diag(fmatd./fnorm); 
driver.oem.cov = blkdiag(fmat,wmat,tmat,ozmat);
%---------------------------------------------------------------------
% Empirical regularization parameters and switches
driver.oem.reg_type = 'reg_and_cov'; % 'reg_and_cov','cov','reg' are other choices
% Separate reg weights for water, temperature profiles
driver.oem.alpha_temp =  cov_set(11);
driver.oem.alpha_water = cov_set(12);
driver.oem.alpha_ozone = cov_set(13);
