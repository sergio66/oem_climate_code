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
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E-1     20*1E-1  20*1E-1];  %% try strat only x3   unc, too higly damped??
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E-2     20*1E-2  20*1E-2];  %% try strat only x3   unc

