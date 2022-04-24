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
if iNlays_retrieve < 90
  trpi = floor(trpi/(100/iNlays_retrieve));
end

%%%%%%%%%%%%%%%%%%%%
%% see ../Strow/strow_override_defaults_latbins_AIRS.m

%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLABCODE/oem_pkg_run/AIRS_AllSky_97_O3jacs_Apr2017/strow_build_cov_matrix.m
%%         lc   ct.lev1  ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide alpha_T  alpha_w  alpha_oz
%%               TROP    STRAT               TROP     STRAT
%% if lc is tiny, mat_od ~ diagnol; if lc is large, mat_od is very broad
%%
%% if cX_wide is tiny, transition is very smooth and gradual ... if cX_wide is large, almost step function transition from trop to strat
%% reg ~ Scov + S_tikonov ~ (sig_q)^-1 + alpha L1'L1       where sig_q is uncertainty in param (eg K^2 for temperature)
%%   if (sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization
%%      (sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov    regularization
%%   if sig_q -> 0   then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%      sig_q -> INF then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%   if alpha -> 0   then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%      alpha -> INF then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%
%%               sigT_t    sigT_s                sigWV_t   sigWV_s             sigO3_t   sigO3_s
%%         lc   ct.lev1  ct.lev2   ct_wide     cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2    coz_wide  alpha_T  alpha_w  alpha_oz

%%         lc   ct.lev1   ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2 coz_wide  alpha_T  alpha_w  alpha_oz

if driver.i16daytimestep > 0
  %% worked great for anomalies!!!
  %%% topts.obs_corr_matrix = -1, 10,20 lays, topts.invtype = 1,3
  cov_set = [1.0  0.005/2   0.005/2   1/2       0.005/25     0.005/25   1/2      0.001/2   0.001/0.75     1/2        1E-5     1E-5  1E-5]; %pretty good for obs  YEAH YEAH YEAH
  cov_set = [1.0  0.05/2    0.05/2    1/2       0.15/25      0.15/25    1/2      0.05/2    0.05/0.75      1/2        1E-1     1E-1  1E-1]; %works pretty well for topts.obs_corr_matrix = -1  and 10 lays, 
                                                                 %OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly less wiggly)  <<<<< used as DEFAULT starting for all the SAVE_blah dirs, qrenorm ~= 1 >>>>
  
elseif driver.i16daytimestep < 0
  %% trying to get goooooood trends
  %%         lc   ct.lev1   ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2 coz_wide  alpha_T  alpha_w  alpha_oz
  cov_set = [1.0  0.05/2/2    0.05/2/2    1/2       0.15/25/2      0.15/25/2    1/2      0.05/2/2    0.05/0.75/2      1/2        2*1E-1     2*1E-1  2*1E-1]; 
  cov_set = [1.0  0.05/2*1000.0  0.05/2*1000.0  1/2       0.15/25*1000.0     0.15/25*1000.0   1/2      0.05/2*1000.0   0.05/0.75*1000.0    1/2        20*1E-3     20*1E-1  20*1E-1];  
  cov_set = [1.0  0.05/2*1000.0  0.05/2*1000.0  1/2       0.15/25*1000.0     0.15/25*1000.0   1/2      0.15/25*1000.0  0.15/25*1000.0      1/2        20*1E-2     20*1E-2  20*1E-2];  
  cov_set = [1.0  0.05/2*1.0E4  0.05/2*1.0E4    1/2       0.15/25*1.0E4     0.15/25*1.0E4     1/2      0.15/25*1.0E4   0.15/25*1.0E4       1/2        20*1E-3     20*1E-3  20*1E-3];  
  cov_set = [1.0  0.05          0.05            1/2       0.15/25           0.15/25           1/2      0.15/25         0.15/25             1/2        20*1E-3     20*1E-3  20*1E-3];  

  %% from ../AIRS_new_clear_scan_August2019_AMT2020PAPER/build_cov_matrices.m
  cov_set = [1.0  0.05/2/2    0.05/2/2    1/2       0.15/25/2      0.15/25/2    1/2      0.05/2/2    0.05/0.75/2      1/2        2*1E-1     2*1E-1  2*1E-1];
  cov_set = [1.0  0.05/2/20    0.05/2/20   1/2       0.15/25/20      0.15/25/20    1/2      0.05/2/20    0.05/0.75/20      1/2        20*1E-1     20*1E-1  20*1E-1];
  cov_set = [1.0  0.05/2*2    0.05/2*2    1/2       0.15/25*2      0.15/25*2    1/2      0.05/2*2    0.05/0.75*2      1/2        2*1E-1     2*1E-1  2*1E-1];

  cov_set = [1.0  0.05/2*2    0.05/2*2    1/2       0.15/25*2      0.15/25*2    1/2      0.15/25*2    0.15/25*2      1/2        2*1E-1     2*1E-1  2*1E-1];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cov_set = [1.0  0.05          0.05            1/2       0.15/25           0.15/25           1/2      0.15/25         0.15/25             1/2        20*1E-6     20*1E-6  20*1E-6];  %worked gret AIRS STM, 20 layers though WV O3 mebbe a little too tight, Q00, *8-16!! not so good for very cloudy eg Q04
  cov_set = [1.0  0.1          0.1            1/2       0.15/10           0.15/10           1/2      0.15/10         0.15/10             1/2        20*1E-8     20*1E-8  20*1E-8];  %worked GOOD for iDebugRatesUseNWP = 3X,%X (ERA5 and AIRS L3reconstructed rates, though I think T(z) rates could be a little looser)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  cov_set = [1.0  0.2          0.2            1/2       0.15/1           0.15/1           1/2      0.15/1         0.15/1             1/2        20*1E-9     20*1E-9  20*1E-9];  %even better, getting closer to AIRS L3, but maybe too loose
%  cov_set = [1.0  0.2          0.2            1/2       0.15/1           0.15/1           1/2      0.15/1         0.15/1             1/2        10*1E-9     10*1E-9  10*1E-9];  % works well for iQuant == 8 to 16, lousy for iQuant 4
%  cov_set = [1.0  0.1          0.1            1/2       0.15/2           0.15/2           1/2      0.15/2         0.15/2             1/2        20*1E-8     20*1E-8  20*1E-8];  %not bad
  %cov_set = [1.0  0.025/2       0.025/2         1/2       0.15/200           0.15/200           1/2      0.15/200         0.15/200             1/2        20*1E-4     20*1E-4  20*1E-4];  %33 layers, tighten things up

  %% worked GOOD for iDebugRatesUseNWP = 32,52,62 (ERA5 and AIRS L3reconstructed rates, though I think T(z) rates could be a little looser), 
  %% trend for T ~ 0.01 K/yr, WVfrac and O3frac ~ 1e-3/yr while uncertainties for T ~ 0.04K/yr, WVfrac and O3frac ~ 4e-3/yr -- so unc ~ x4 the retrieved value, till Oct 15, 2021
  cov_set = [1.0  0.2          0.2            1/2       0.15/10           0.15/10           1/2      0.15/10         0.15/10             1/2        20*1E-8     20*1E-8  20*1E-8];    %% some spikes when fmatd is changed so CO2 sticks to 2.2 ppm
  cov_set = [1.0  0.1          0.1            1/2       0.15/20           0.15/20           1/2      0.15/20         0.15/20             1/2        20*1E-6     20*1E-6  20*1E-6];    %% some spikes when fmatd is changed so CO2 sticks to 2.2 ppm

  %% increase WV,O3 uncertainty in covariances
  %% trend for T ~ 0.01 K/yr, WVfrac and O3frac ~ 1e-3/yr while uncertainties for T ~ 0.04K/yr, WVfrac and O3frac ~ 3e-2/yr -- so unc ~ x40 the retrieved value OOOOOOOOOOOOPS, Oct 16, 2021
%  cov_set = [1.0  0.2          0.2            1/2       0.15           0.15           1/2      0.15         0.15             1/2        20*1E-12     20*1E-12  20*1E-12];  

  %% reduce T,WV,O3 uncertainty in covariances, makes WV DOF too small?????
  %% trend for T ~ 0.01 K/yr, WVfrac and O3frac ~ 1e-3/yr while uncertainties for T ~ 0.03K/yr, WVfrac and O3frac ~ 1e-3/yr -- so unc ~ x3 the retrieved value, actually pretty good, Oct 16, 2021 
  %%                                                                                     and Feb 4, 2022 (very good simulated ERA5 results)
  %% works well with sticking CO2 to 2.2 ppm/yr, no spikes
  cov_set = [1.0  0.025         0.05            1/2       0.15/50           0.15/50           1/2      0.15/50         0.15/50             1/2        20*1E-7     20*1E-7  20*1E-7];  %% ok   excellent simulated ERA5 spectral rates Feb 4, 2022
  cov_set = [1.0  0.05*1        0.05*1          1/2       0.15/50*1         0.15/50*1         1/2      0.15/50*1       0.15/50*1           1/2        20*1E-7     20*1E-7  20*1E-7];  %% very excellent simulated ERA5 spectral rates Feb 4, 2022 till Feb 15, 2022

  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02              0.02              1/2      0.02            0.02                1/2        20*1E-7     20*1E-7  20*1E-7];  %% try x100 unc, Feb 16 2022-Apr7,2022 :  great JPLMay 2022 talk!!!
  cov_set = [1.0  0.05*3        0.05*3          1/2       0.15/50*3         0.15/50*3         1/2      0.15/50*3       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% try x3   unc, Feb 16 2022-Apr7,2022 :  great but maybe still constricts WV/O3

  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% try strat only x3   unc
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E+2     20*1E+2  20*1E+2];  %% try strat only x3   unc

%%% after AIRS STM 2021, testing quantile 16
%%%  cov_set = [1.0  0.05          0.05            1/2       0.15/10           0.15/10           1/2      0.15/10         0.15/10             1/2        20*1E-6     20*1E-6  20*1E-6];  %
%%%  cov_set = [1.0  0.05          0.05            1/2       0.15/05           0.15/05           1/2      0.15/05         0.15/05             1/2        20*1E-6     20*1E-6  20*1E-6];  %

end

fprintf(1,'cov_set lc = %8.6f \n',cov_set(1));
junk = cov_set([2 3 4  11]); fprintf(1,'      T  : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk)
junk = cov_set([5 6 7  12]); fprintf(1,'      WV : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk);
junk = exp10(junk(1:2))-1;   fprintf(1,'           percent sig_trop  sig_strat  = %8.6e %8.6e \n',junk*100);
junk = cov_set([8 9 10 13]); fprintf(1,'      O3 : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk)
junk = exp10(junk(1:2))-1;   fprintf(1,'           percent sig_trop  sig_strat  = %8.6e %8.6e \n',junk*100);
%%%%%%%%%%%%%%%%%%%%%%%%%

% Relative off-diagonal
l_c = cov_set(1);;
mat_od = exp(-mat_od.^2./(1*l_c^2));

   iOffX = 10;
   ct(ix).trans1 = trpi(ix);
   ct(ix).trans2 = trpi(ix) + -floor(iOffX/(100/iNlays_retrieve));; %% offset by 10 (when using 100 AIRS layers) or 2 (when using 20 FAT layers)

   ct(ix).trans1 = ct(ix).trans1-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)
   ct(ix).trans2 = ct(ix).trans2-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)

   ct(ix).lev1 = cov_set(2);
   ct(ix).lev2 = cov_set(3);
   ct(ix).lev3 = ct(ix).lev2;
   ct(ix).width1 = cov_set(4);
   ct(ix).width2 = ct(ix).width1;

   iOffX = 10;
%   cw(ix).trans1 = trpi(ix);
%   cw(ix).trans2 = trpi(ix) + -floor(iOffX/(100/iNlays_retrieve));       %% offset by 10 (when using 100 AIRS layers) or 2 (when using 20 FAT layers)
   cw(ix).trans1 = trpi(ix) + -floor(iOffX/(100/iNlays_retrieve));  %% new move towards ground (+1) away from gnd/to TOA (-)
   cw(ix).trans2 = trpi(ix) + -floor(iOffX/(100/iNlays_retrieve));  %% new move towards ground (+1) away from gnd/to TOA (-)

   cw(ix).trans1 = cw(ix).trans1-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)
   cw(ix).trans2 = cw(ix).trans2-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)

   cw(ix).lev1 = cov_set(5);
   cw(ix).lev2 = cov_set(6);
   cw(ix).lev3 = cw(ix).lev2;
   cw(ix).width1 = cov_set(7);
   cw(ix).width2 = cw(ix).width1;

   iOffX = 10;
   coz(ix).trans1 = trpi(ix);
   coz(ix).trans2 = trpi(ix) + floor(iOffX/(100/iNlays_retrieve));       %% offset by 10 (when using 100 AIRS layers) or 2 (when using 20 FAT layers)
   coz(ix).trans1 = coz(ix).trans1-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)
   coz(ix).trans2 = coz(ix).trans2-floor(4*iOffX/(100/iNlays_retrieve)); %% new move towards ground (+1) away from gnd/to TOA (-)
   coz(ix).lev1 = cov_set(8);
   coz(ix).lev2 = cov_set(9);
   coz(ix).lev3 = coz(ix).lev2;
   coz(ix).width1 = cov_set(10);
   coz(ix).width2 = coz(ix).width1;

% Temperature level uncertainties, then scaled and squared
if ~exist('iFixTz_NoFit')
  tunc = cov2lev(ct(ix),driver.jacobian.numlays);
  t_sigma = (tunc./tnorm);
  tmat = (t_sigma'*t_sigma).*mat_od;
  driver.oem.tunc = tunc;
elseif exist('iFixTz_NoFit','var') & ~strcmp(driver.rateset.ocb_set,'obs')
  tunc = cov2lev(ct(ix),driver.jacobian.numlays);
  t_sigma = (tunc./tnorm);
  tmat = (t_sigma'*t_sigma).*mat_od;
  driver.oem.tunc = tunc;
elseif exist('iFixTz_NoFit','var') & strcmp(driver.rateset.ocb_set,'obs')
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixTz_NoFit > 0
    disp('making sure T is not fitted!!!!')

    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e10;  %% THIS COMBO WORKS, CO2, CFC,N2O, CH4 good but still see residual dT/dt
    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e16;  %% THIS COMBO WORKS, much smaller dT/dt, CO2 good, CFC,O3,WV bad
    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e12;  %% 

    driver.oem.tunc = []; cov_set(11) = 1.0e10;
  end
end

% ozone level uncertainties, then scaled and squared
if ~exist('iFixO3_NoFit')
  ozunc = cov2lev(coz(ix),driver.jacobian.numlays);
  oz_sigma = (ozunc./oznorm);
  ozmat = (oz_sigma'*oz_sigma).*mat_od;
  driver.oem.ozunc = ozunc;
elseif exist('iFixO3_NoFit','var') & ~strcmp(driver.rateset.ocb_set,'obs') & ~strcmp(driver.rateset.ocb_set,'cal')
  ozunc = cov2lev(coz(ix),driver.jacobian.numlays);
  oz_sigma = (ozunc./oznorm);
  ozmat = (oz_sigma'*oz_sigma).*mat_od;
  driver.oem.ozunc = ozunc;
elseif exist('iFixO3_NoFit','var') & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixO3_NoFit >= 0
    disp('making sure O3 is not fitted!!!!')

    %driver.oem.o3unc = ones(size(o3unc)) * 1.0e-16;  cov_set(11) = 1.0e10;  %% THIS COMBO WORKS, CO2, CFC,N2O, CH4 good but still see residual dT/dt
    %driver.oem.o3unc = ones(size(o3unc)) * 1.0e-16;  cov_set(11) = 1.0e16;  %% THIS COMBO WORKS, much smaller dT/dt, CO2 good, CFC,O3,WV bad
    %driver.oem.o3unc = ones(size(o3unc)) * 1.0e-16;  cov_set(11) = 1.0e12;  %% 

    driver.oem.o3unc = []; cov_set(13) = 1.0e10;
  end
end

% Water level uncertainties, then scaled and squared
wunc = cov2lev(cw(ix),driver.jacobian.numlays);
w_sigma = (wunc./wnorm);
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

%keyboard_nowindow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) N2O(ppb) CH4(ppb) CFC11(ppt) CFC12 (ppt) Tsurf(K)    cng1 cng2 cpsz1 cpsz2]   
if settings.co2lays == 1
  fmatd = [2     0.1       2      1     1            0.1/3]*3;   %% scanty test time ONLY done Tue Jul 16 15:45:49 EDT 2019
  fmatd = [2     0.1       2      0.25  0.1          0.1];       %% try this, CFC --> -1 ppm/yr, pretty good CFC!!      : files 8/30/19  *newchans_newkcartajacs_finitediffV3_5_tweakCFCunc.mat 
  fmatd = [2     1.0       5      2     2            0.1];       %% try loosening  so all gases have unc ~ 1.0 year growth : files 8/30/19  *newchans_newkcartajacs_finitediffV3_5_allTraceGasUnc_1yearrate.mat
  fmatd = [2     1.0       5      2     2            0.1]/10;    %% try tightening so all gases have unc ~ 0.1 year growth : files 9/03/19  *newchans_newkcartajacs_finitediffV3_5_allTraceGasUnc_0p1yearrate.mat
  fmatd = [2     0.1       2      0.02    0.02       0.1];       %% 99.99999999999999999999% of the time, works pretty well but CFC11,CFC12 rates are too high by x2  -3.4 ppm/yr (July/Aug 2019) >>> ORIG, BEST
  fmatd = [2     0.1       2      1     1            1];         %% 99.99999999999999999999% of the time, works pretty well but CFC11,CFC12 rates are too high by x2  -3.4 ppm/yr (July/Aug 2019) >>> ORIG, BEST
elseif settings.co2lays == 3
  fmatd = [2 2 2 0.1       2      0.02     0.02            0.1];
end
%fmatd = ones(size(fmatd))*1e-4;
fmatd(6) = 0.1;

if settings.set_tracegas == 1 & settings.co2lays == 1 & driver.i16daytimestep < 0
  %fmatd(1:3) = fmatd(1:3)*1e-7;   %% xb(1:3) (Co2/N2o/Ch4) so we should not change those values .. gives large spectral bias, great gas rates eg CO2=2.2, terrible T rates  ******
  %fmatd(1:3) = fmatd(1:3)*1e-5;   %% xb(1:3) (Co2/N2o/Ch4) so we should not change those values .. gives large spectral bias, great gas rates eg CO2=2.2, terrible T rates  ******

  %fmatd(1:3) = fmatd(1:3)*1e-1;    %% gives small spectral bias, lousy SARTA trace gas rates eg CO2=1.0, semi ok KCARTA trace gas rates   *** >>>  WORKS FOR AIRS STM Oct 2020 and May 2021 and AUG 2021<<< ****  
  %fmatd(1:3) = fmatd(1:3)*1.0;     %% gives small spectral bias, lousy SARTA trace gas rates eg CO2=1.0, semi ok KCARTA trace gas rates   *** >>>  WORKS FOR AIRS STM Oct 2020 and May 2021 <<< ****


  %fmatd(1:3) = fmatd(1:3)*1e-2;    %% gives small spectral bias, TESTING SEPT 2021, works prety well except in tropics where CO2 becomes 1.9 ppm/yr
  fmatd(1:3) = fmatd(1:3)*1e-3;     %% so add in another factor of 10 on Oct 22, 2021

elseif settings.set_tracegas == 1 & settings.co2lays == 3 & driver.i16daytimestep < 0
  fmatd(1:5) = fmatd(1:5)*0.0000001;   %% have put in xb(1:3) so we should not change those values .. recall 1,2,3 = CO2/N2O/CH4 and 4/5 are cld1,cld2
  fmatd(1:5) = fmatd(1:5)*1e-3;    %% xb(1:3) (Co2/N2o/Ch4) so we should not change those values .. gives small spectral bias, lousy SARTA trace gas rates eg CO2=1.0, semi ok KCARTA trace gas rates  ******
end

%% NEW Aug 16, 2020
%if driver.i16daytimestep < 0
%  fmatd = 1*[2     0.1       2      1     1            0.1];       %% 99.99999999999999999999% of the time, works pretty well but CFC11,CFC12 rates are too high by x2  -3.4 ppm/yr (July/Aug 2019) >>> ORIG, BEST
%end

fmatd
%keyboard_nowindow

%{
%% before Oct 21, 2021
if settings.iFixTG_NoFit(1) > 0
  disp('setting uncertainties for some tracegases to be 0, setting a priori xb and jacs for these gases also to be 0')
  xb(settings.iFixTG_NoFit) = 0.0;
  aux.xb(settings.iFixTG_NoFit) = 0.0;
  driver.oem.xb(settings.iFixTG_NoFit) = 0.0;
  fmatd(settings.iFixTG_NoFit) = fmatd(settings.iFixTG_NoFit) * 0.0000000000000001; %% this way we do not change values
  
  m_ts_jac(:,settings.iFixTG_NoFit) = eps;
  aux.m_ts_jac(:,settings.iFixTG_NoFit) = eps;
end
%}
%% after Oct 21, 2021
if settings.iFixTG_NoFit(1) > 0
  disp('setting uncertainties for some tracegases to be 0, NOT setting a priori xb, or jac, for these gases also to be 0')
  %xb(settings.iFixTG_NoFit) = 0.0;
  %aux.xb(settings.iFixTG_NoFit) = 0.0;
  %driver.oem.xb(settings.iFixTG_NoFit) = 0.0;
  fmatd(settings.iFixTG_NoFit) = fmatd(settings.iFixTG_NoFit) * 0.0000000000000001; %% this way we do not change values
  
  %m_ts_jac(:,settings.iFixTG_NoFit) = eps;
  %aux.m_ts_jac(:,settings.iFixTG_NoFit) = eps;
end

fmat  = diag(fmatd./fnorm); 
if exist('tmat','var') & exist('ozmat','var')
  driver.oem.cov = blkdiag(fmat,wmat,tmat,ozmat);
elseif exist('ozmat','var')
  driver.oem.cov = blkdiag(fmat,wmat,ozmat);
elseif exist('tmat','var')
  driver.oem.cov = blkdiag(fmat,wmat,tmat);
end

if topts.tie_sst_lowestlayer > 0 & exist('tmat','var')
  %% now tie together surface temp with lowest layers
  wah1 = driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.scalar_i(end));
  wah2 = driver.oem.cov(driver.jacobian.temp_i(end),driver.jacobian.temp_i(end));
  %driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.temp_i(end-1):driver.jacobian.temp_i(end)) = sqrt(abs(wah1*wah2));
  %driver.oem.cov(driver.jacobian.temp_i(end-1):driver.jacobian.temp_i(end),driver.jacobian.scalar_i(end)) = sqrt(abs(wah1*wah2));
  driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.temp_i(end):driver.jacobian.temp_i(end)) = sqrt(abs(wah1*wah2));
  driver.oem.cov(driver.jacobian.temp_i(end):driver.jacobian.temp_i(end),driver.jacobian.scalar_i(end)) = sqrt(abs(wah1*wah2));
end

%---------------------------------------------------------------------
% Empirical regularization parameters and switches
%% from oem_run/rodgers.m
%% switch driver.oem.reg_type
%%   case 'reg_and_cov'
%%     r = rcov + rc;
%%   case 'reg'
%%     r = rc;
%%   case 'cov'
%%     r = rcov;
%%   otherwise
%%     disp('Incorrect choice driver.oem.reg_type')
%% end
%%

driver.oem.reg_type = 'cov';         % 'reg_and_cov','cov','reg' are other choices
driver.oem.reg_type = 'reg';         % 'reg_and_cov','cov','reg' are other choices
driver.oem.reg_type = 'reg_and_cov'; % 'reg_and_cov','cov','reg' are other choices DEFAULT

% Separate reg weights for water, temperature, ozone profiles
driver.oem.alpha_water = cov_set(12);

if ~exist('iFixTz_NoFit')
  driver.oem.alpha_temp =  cov_set(11);
elseif exist('iFixTz_NoFit','var') & ~strcmp(driver.rateset.ocb_set,'obs')
  driver.oem.alpha_temp =  cov_set(11);
elseif exist('iFixTz_NoFit','var') & strcmp(driver.rateset.ocb_set,'obs')
  disp('oops iFixTz_NoFit exists, for driver.rateset.ocb_set = obs so not setting driver.oem.alpha_temp')
end

if ~exist('iFixO3_NoFit')
  driver.oem.alpha_ozone =  cov_set(13);
elseif exist('iFixO3_NoFit','var') & ~strcmp(driver.rateset.ocb_set,'obs') & ~strcmp(driver.rateset.ocb_set,'cal')
  driver.oem.alpha_ozone =  cov_set(13);
elseif exist('iFixO3_NoFit','var') & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
  disp('oops iFixO3_NoFit exists, for driver.rateset.ocb_set = obs so not setting driver.oem.alpha_ozone')
end

driver.oem.cov_set = cov_set;
driver.oem.fmat    = fmat;
