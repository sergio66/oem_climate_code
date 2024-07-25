dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN';
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN_Temp';

fprintf(1,'will be saving plots/figs etc to %s \n',dir0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('printing compare_SKT_trends_Day_vs_Night.m');
figure(20); sergioprintfig([dir0 '/skt_Day_versus_Night_4panel']);
figure(21); sergioprintfig([dir0 '/skt_Day_minus_Night_4panel']);
figure(22); sergioprintfig([dir0 '/skt_Day_Night_avg_6panel']);
figure(23); sergioprintfig([dir0 '/skt_Day_Night_meantrend_over6']);
figure(24); sergioprintfig([dir0 '/skt_Day_Night_stddevtrend_over6']);
figure(27); sergioprintfig([dir0 '/skt_Day_Night_meantrend_umbc']);
figure(28); sergioprintfig([dir0 '/skt_Day_Night_stddevtrend_umbc']);

disp('printing compare_colWV_trends_Day_vs_Night.m');
figure(30); sergioprintfig([dir0 '/colwv_Day_versus_Night_4panel']);
figure(31); sergioprintfig([dir0 '/colwv_Day_minus_Night_4panel']);
figure(32); sergioprintfig([dir0 '/colwv_Day_Night_avg_5panel']);
figure(33); sergioprintfig([dir0 '/colwv_Day_Night_meantrend_over5']);
figure(34); sergioprintfig([dir0 '/colwv_Day_Night_stddevtrend_over5']);
figure(37); sergioprintfig([dir0 '/colwv_Day_Night_meantrend_umbc']);
figure(38); sergioprintfig([dir0 '/colwv_Day_Night_stddevtrend_umbc']);

disp('printing compare_T_trends_Day_vs_Night.m');
figure(40); sergioprintfig([dir0 '/tz_Day_versus_Night_4panel']);
figure(41); sergioprintfig([dir0 '/tz_Day_minus_Night_4panel']);
figure(42); sergioprintfig([dir0 '/tz_Day_Night_avg_5panel']);
figure(43); sergioprintfig([dir0 '/tz_Day_Night_meantrend_over5']);
figure(44); sergioprintfig([dir0 '/tz_Day_Night_stddevtrend_over5']);
figure(47); sergioprintfig([dir0 '/tz_Day_Night_meantrend_umbc']);
figure(48); sergioprintfig([dir0 '/tz_Day_Night_stddevtrend_umbc']);

disp('printing compare_WV_trends_Day_vs_Night.m');
figure(50); sergioprintfig([dir0 '/wvz_Day_versus_Night_4panel']);
figure(51); sergioprintfig([dir0 '/wvz_Day_minus_Night_4panel']);
figure(52); sergioprintfig([dir0 '/wvz_Day_Night_avg_5panel']);
figure(53); sergioprintfig([dir0 '/wvz_Day_Night_meantrend_over5']);
figure(54); sergioprintfig([dir0 '/wvz_Day_Night_stddevtrend_over5']);
figure(57); sergioprintfig([dir0 '/wvz_Day_Night_meantrend_umbc']);
figure(58); sergioprintfig([dir0 '/wvz_Day_Night_stddevtrend_umbc']);

disp('printing compare_RH_trends_Day_vs_Night.m');
figure(60); sergioprintfig([dir0 '/rh_Day_versus_Night_4panel']);
figure(61); sergioprintfig([dir0 '/rh_Day_minus_Night_4panel']);
figure(62); sergioprintfig([dir0 '/rh_Day_Night_avg_5panel']);
figure(63); sergioprintfig([dir0 '/rh_Day_Night_meantrend_over5']);
figure(64); sergioprintfig([dir0 '/rh_Day_Night_stddevtrend_over5']);
figure(67); sergioprintfig([dir0 '/rh_Day_Night_meantrend_umbc']);
figure(68); sergioprintfig([dir0 '/rh_Day_Night_stddevtrend_umbc']);

disp('printing compare_RHsurf_trends_Day_vs_Night.m');
figure(70); sergioprintfig([dir0 '/rhsurf_Day_versus_Night_4panel']);
figure(71); sergioprintfig([dir0 '/rhsurf_Day_minus_Night_4panel']);
figure(72); sergioprintfig([dir0 '/rhsurf_Day_Night_avg_5panel']);
figure(73); sergioprintfig([dir0 '/rhsurf_Day_Night_meantrend_over5']);
figure(74); sergioprintfig([dir0 '/rhsurf_Day_Night_stddevtrend_over5']);
figure(77); sergioprintfig([dir0 '/rhsurf_Day_Night_meantrend_umbc']);
figure(78); sergioprintfig([dir0 '/rhsurf_Day_Night_stddevtrend_umbc']);

disp('printing misc.m');
figure(81); sergioprintfig([dir0 '/dST_dt_land_ocean_zonal_allmodels']);
figure(84); sergioprintfig([dir0 '/dST_dt_ocean_zonal_allmodels']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(90); sergioprintfig([dir0 '/plot_avg_over_6_models_mean_std_skt_trends']);
figure(91); sergioprintfig([dir0 '/plot_avg_over_5_models_mean_std_mmw_trends'])
figure(92); sergioprintfig([dir0 '/mmw_trends_5results'])
figure(93); sergioprintfig([dir0 '/plot_avg_over_5_models_mean_std_T_trends']);
figure(94); sergioprintfig([dir0 '/plot_avg_over_5_models_mean_std_WVfrac_trends'])
figure(95); sergioprintfig([dir0 '/plot_avg_over_5_models_mean_std_RH_trends'])
figure(96); sergioprintfig([dir0 '/plot_avg_over_5_models_mean_std_RHSURF_trends'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103); sergioprintfig([dir0 '/ilr_trends_zonal']);
figure(104); sergioprintfig([dir0 '/ilr_trends_umbc']);



