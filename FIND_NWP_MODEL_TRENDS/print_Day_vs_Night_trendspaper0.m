dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN';
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN_Temp';

fprintf(1,'will be saving plots/figs etc to %s \n',dir0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compare_SKT_trends_Day_vs_Night.m
figure(20); aslprint_asis([dir0 '/skt_Day_versus_Night_4panel.pdf']);
figure(21); aslprint_asis([dir0 '/skt_Day_minus_Night_4panel.pdf']);
figure(22); aslprint_asis([dir0 '/skt_Day_Night_avg_6panel.pdf']);
figure(23); aslprint_asis([dir0 '/skt_Day_Night_meantrend_over6.pdf']);
figure(24); aslprint_asis([dir0 '/skt_Day_Night_stddevtrend_over6.pdf']);
figure(27); aslprint_asis([dir0 '/skt_Day_Night_meantrend_umbc.pdf']);
figure(28); aslprint_asis([dir0 '/skt_Day_Night_stddevtrend_umbc.pdf']);

%% compare_colWV_trends_Day_vs_Night.m
figure(30); aslprint_asis([dir0 '/colwv_Day_versus_Night_4panel.pdf']);
figure(31); aslprint_asis([dir0 '/colwv_Day_minus_Night_4panel.pdf']);
figure(32); aslprint_asis([dir0 '/colwv_Day_Night_avg_5panel.pdf']);
figure(33); aslprint_asis([dir0 '/colwv_Day_Night_meantrend_over5.pdf']);
figure(34); aslprint_asis([dir0 '/colwv_Day_Night_stddevtrend_over5.pdf']);
figure(37); aslprint_asis([dir0 '/colwv_Day_Night_meantrend_umbc.pdf']);
figure(38); aslprint_asis([dir0 '/colwv_Day_Night_stddevtrend_umbc.pdf']);

%% compare_T_trends_Day_vs_Night.m
figure(40); aslprint_asis([dir0 '/tz_Day_versus_Night_4panel.pdf']);
figure(41); aslprint_asis([dir0 '/tz_Day_minus_Night_4panel.pdf']);
figure(42); aslprint_asis([dir0 '/tz_Day_Night_avg_5panel.pdf']);
figure(43); aslprint_asis([dir0 '/tz_Day_Night_meantrend_over5.pdf']);
figure(44); aslprint_asis([dir0 '/tz_Day_Night_stddevtrend_over5.pdf']);
figure(47); aslprint_asis([dir0 '/tz_Day_Night_meantrend_umbc.pdf']);
figure(48); aslprint_asis([dir0 '/tz_Day_Night_stddevtrend_umbc.pdf']);

%% compare_WV_trends_Day_vs_Night.m
figure(50); aslprint_asis([dir0 '/wvz_Day_versus_Night_4panel.pdf']);
figure(51); aslprint_asis([dir0 '/wvz_Day_minus_Night_4panel.pdf']);
figure(52); aslprint_asis([dir0 '/wvz_Day_Night_avg_5panel.pdf']);
figure(53); aslprint_asis([dir0 '/wvz_Day_Night_meantrend_over5.pdf']);
figure(54); aslprint_asis([dir0 '/wvz_Day_Night_stddevtrend_over5.pdf']);
figure(57); aslprint_asis([dir0 '/wvz_Day_Night_meantrend_umbc.pdf']);
figure(58); aslprint_asis([dir0 '/wvz_Day_Night_stddevtrend_umbc.pdf']);

%% compare_RH_trends_Day_vs_Night.m
figure(60); aslprint_asis([dir0 '/rh_Day_versus_Night_4panel.pdf']);
figure(61); aslprint_asis([dir0 '/rh_Day_minus_Night_4panel.pdf']);
figure(62); aslprint_asis([dir0 '/rh_Day_Night_avg_5panel.pdf']);
figure(63); aslprint_asis([dir0 '/rh_Day_Night_meantrend_over5.pdf']);
figure(64); aslprint_asis([dir0 '/rh_Day_Night_stddevtrend_over5.pdf']);
figure(67); aslprint_asis([dir0 '/rh_Day_Night_meantrend_umbc.pdf']);
figure(68); aslprint_asis([dir0 '/rh_Day_Night_stddevtrend_umbc.pdf']);

%% compare_RHsurf_trends_Day_vs_Night.m
figure(70); aslprint_asis([dir0 '/rhsurf_Day_versus_Night_4panel.pdf']);
figure(71); aslprint_asis([dir0 '/rhsurf_Day_minus_Night_4panel.pdf']);
figure(72); aslprint_asis([dir0 '/rhsurf_Day_Night_avg_5panel.pdf']);
figure(73); aslprint_asis([dir0 '/rhsurf_Day_Night_meantrend_over5.pdf']);
figure(74); aslprint_asis([dir0 '/rhsurf_Day_Night_stddevtrend_over5.pdf']);
figure(77); aslprint_asis([dir0 '/rhsurf_Day_Night_meantrend_umbc.pdf']);
figure(78); aslprint_asis([dir0 '/rhsurf_Day_Night_stddevtrend_umbc.pdf']);

figure(81); aslprint_asis([dir0 '/dST_dt_land_ocean_zonal_allmodels.pdf']);
figure(84); aslprint_asis([dir0 '/dST_dt_ocean_zonal_allmodels.pdf']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(90); aslprint_asis([dir0 '/plot_avg_over_6_models_mean_std_skt_trends.pdf']);
figure(91); aslprint_asis([dir0 '/plot_avg_over_5_models_mean_std_mmw_trends.pdf'])
figure(92); aslprint_asis([dir0 '/mmw_trends_5results.pdf'])
figure(93); aslprint_asis([dir0 '/plot_avg_over_5_models_mean_std_T_trends.pdf']);
figure(94); aslprint_asis([dir0 '/plot_avg_over_5_models_mean_std_WVfrac_trends.pdf'])
figure(95); aslprint_asis([dir0 '/plot_avg_over_5_models_mean_std_RH_trends.pdf'])
figure(96); aslprint_asis([dir0 '/plot_avg_over_5_models_mean_std_RHSURF_trends.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103); aslprint_asis([dir0 '/ilr_trends_zonal.pdf']);
figure(104); aslprint_asis([dir0 '/ilr_trends_umbc.pdf']);



