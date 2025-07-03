dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN/';
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN_Temp/';
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/';

fprintf(1,'will be saving plots/figs etc to %s \n',dir0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('printing compare_SKT_trends_Day_vs_Night.m');
figure(20); sergioprintfig([dir0 '/nosmooth_skt_Day_versus_Night_4panel']);         %%%% XXXX need to save data    fig B1 (appendix)
figure(21); sergioprintfig([dir0 '/nosmooth_skt_Day_minus_Night_4panel']);
figure(22); sergioprintfig([dir0 '/nosmooth_skt_Day_Night_avg_6panel']);
figure(23); sergioprintfig([dir0 '/nosmooth_skt_Day_Night_meantrend_over6']);
figure(24); sergioprintfig([dir0 '/nosmooth_skt_Day_Night_stddevtrend_over6']);
figure(27); sergioprintfig([dir0 '/nosmooth_skt_Day_Night_meantrend_umbc']);
figure(28); sergioprintfig([dir0 '/nosmooth_skt_Day_Night_stddevtrend_umbc']);

disp('printing compare_colWV_trends_Day_vs_Night.m');
figure(30); sergioprintfig([dir0 '/nosmooth_colwv_Day_versus_Night_4panel']);
figure(31); sergioprintfig([dir0 '/nosmooth_colwv_Day_minus_Night_4panel']);
figure(32); sergioprintfig([dir0 '/nosmooth_colwv_Day_Night_avg_5panel']);
figure(33); sergioprintfig([dir0 '/nosmooth_colwv_Day_Night_meantrend_over5']);
figure(34); sergioprintfig([dir0 '/nosmooth_colwv_Day_Night_stddevtrend_over5']);
figure(37); sergioprintfig([dir0 '/nosmooth_colwv_Day_Night_meantrend_umbc']);
figure(38); sergioprintfig([dir0 '/nosmooth_colwv_Day_Night_stddevtrend_umbc']);

disp('printing compare_T_trends_Day_vs_Night.m');
figure(40); sergioprintfig([dir0 '/nosmooth_tz_Day_versus_Night_4panel']);
figure(41); sergioprintfig([dir0 '/nosmooth_tz_Day_minus_Night_4panel']);
figure(42); sergioprintfig([dir0 '/nosmooth_tz_Day_Night_avg_5panel']);               %%%% XXXX need to save data    fig11
figure(43); sergioprintfig([dir0 '/nosmooth_tz_Day_Night_meantrend_over5']);
figure(44); sergioprintfig([dir0 '/nosmooth_tz_Day_Night_stddevtrend_over5']);
figure(47); sergioprintfig([dir0 '/nosmooth_tz_Day_Night_meantrend_umbc']);
figure(48); sergioprintfig([dir0 '/nosmooth_tz_Day_Night_stddevtrend_umbc']);

disp('printing compare_WV_trends_Day_vs_Night.m');
figure(50); sergioprintfig([dir0 '/nosmooth_wvz_Day_versus_Night_4panel']);
figure(51); sergioprintfig([dir0 '/nosmooth_wvz_Day_minus_Night_4panel']);
figure(52); sergioprintfig([dir0 '/nosmooth_wvz_Day_Night_avg_5panel']);              %%%% XXXX need to save data     fig12
figure(53); sergioprintfig([dir0 '/nosmooth_wvz_Day_Night_meantrend_over5']);
figure(54); sergioprintfig([dir0 '/nosmooth_wvz_Day_Night_stddevtrend_over5']);
figure(57); sergioprintfig([dir0 '/nosmooth_wvz_Day_Night_meantrend_umbc']);
figure(58); sergioprintfig([dir0 '/nosmooth_wvz_Day_Night_stddevtrend_umbc']);

disp('printing compare_RH_trends_Day_vs_Night.m');
figure(60); sergioprintfig([dir0 '/nosmooth_rh_Day_versus_Night_4panel']);
figure(61); sergioprintfig([dir0 '/nosmooth_rh_Day_minus_Night_4panel']);
figure(62); sergioprintfig([dir0 '/nosmooth_rh_Day_Night_avg_5panel']);
figure(63); sergioprintfig([dir0 '/nosmooth_rh_Day_Night_meantrend_over5']);
figure(64); sergioprintfig([dir0 '/nosmooth_rh_Day_Night_stddevtrend_over5']);
figure(67); sergioprintfig([dir0 '/nosmooth_rh_Day_Night_meantrend_umbc']);
figure(68); sergioprintfig([dir0 '/nosmooth_rh_Day_Night_stddevtrend_umbc']);

disp('printing compare_RHsurf_trends_Day_vs_Night.m');
figure(70); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_versus_Night_4panel']);
figure(71); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_minus_Night_4panel']);
figure(72); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_Night_avg_5panel']);
figure(73); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_Night_meantrend_over5']);
figure(74); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_Night_stddevtrend_over5']);
figure(77); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_Night_meantrend_umbc']);
figure(78); sergioprintfig([dir0 '/nosmooth_rhsurf_Day_Night_stddevtrend_umbc']);

disp('printing misc.m');
figure(81); sergioprintfig([dir0 '/nosmooth_dST_dt_land_ocean_zonal_allmodels']);
figure(84); sergioprintfig([dir0 '/nosmooth_dST_dt_ocean_zonal_allmodels']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(90); sergioprintfig([dir0 '/nosmooth_plot_avg_over_6_models_mean_std_skt_trends']);
figure(91); sergioprintfig([dir0 '/nosmooth_plot_avg_over_5_models_mean_std_mmw_trends'])
figure(92); sergioprintfig([dir0 '/nosmooth_mmw_trends_5results'])                                  %%%% XXXX need to save data   fig8
figure(93); sergioprintfig([dir0 '/nosmooth_plot_avg_over_5_models_mean_std_T_trends']);
figure(94); sergioprintfig([dir0 '/nosmooth_plot_avg_over_5_models_mean_std_WVfrac_trends'])
figure(95); sergioprintfig([dir0 '/nosmooth_plot_avg_over_5_models_mean_std_RH_trends'])
figure(96); sergioprintfig([dir0 '/nosmooth_plot_avg_over_5_models_mean_std_RHSURF_trends'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103); sergioprintfig([dir0 '/nosmooth_ilr_trends_zonal']);
figure(104); sergioprintfig([dir0 '/nosmooth_ilr_trends_umbc']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSaveJGR = -1; 
%% << iV and iV6 are from get_umbc_day_night_name.m >>

if iSaveJGR > 0 & iV == 3 & iV3 == 66
  figure(52);
  figure(152); clf
  fig9_WVtrends_noMLS = fig12_WVtrends;
    iFig = 152; 
        profile_plots_2x1x2tiledlayout_tall(fig9_WVtrends_noMLS.rlat,fig9_WVtrends_noMLS.plays100,fig9_WVtrends_noMLS.umbc,fig9_WVtrends_noMLS.airsL3,fig9_WVtrends_noMLS.climcaps,fig9_WVtrends_noMLS.merra2,fig9_WVtrends_noMLS.era5,iFig,fig9_WVtrends_noMLS.plotoptions2x1x2);
  
  %% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig9.mat fig9_WVtrends_noMLS

elseif iSaveJGR > 0 & iV == 3 & iV3 == 6
  %% from compare_SKT_trends_Day_vs_Night.m
  figure(20);
  figure(120); clf
    iFig = 120;
    aslmap_3x4tiledlayout(figB1_sktDNtrends.umbc_D,figB1_sktDNtrends.airsL3_D,figB1_sktDNtrends.climcaps_D,figB1_sktDNtrends.era5_D,...
                          figB1_sktDNtrends.umbc_N,figB1_sktDNtrends.airsL3_N,figB1_sktDNtrends.climcaps_N,figB1_sktDNtrends.era5_N,...
                          figB1_sktDNtrends.umbc_X,figB1_sktDNtrends.airsL3_X,figB1_sktDNtrends.climcaps_X,figB1_sktDNtrends.era5_X,...
                          iFig,figB1_sktDNtrends.plotoptions);
  figB1_sktDNtrends.comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_trends_Day_vs_Night.m  --> compare_SKT_trends_Day_vs_Night.m';
  %% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/figB1.mat figB1_sktDNtrends
  
  %% from compare_T_trends_Day_vs_Night
  figure(42)
  figure(142); clf
    iFig = 142; 
        profile_plots_2x1x2tiledlayout_tall(fig11_Ttrends.rlat,fig11_Ttrends.plays100,fig11_Ttrends.umbc,fig11_Ttrends.airsL3,fig11_Ttrends.climcaps,fig11_Ttrends.merra2,fig11_Ttrends.era5,iFig,fig11_Ttrends.plotoptions2x1x2);
  %% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig11.mat fig11_Ttrends
  
  figure(52);
  figure(152); clf
    iFig = 152; 
        profile_plots_2x1x2tiledlayout_tall(fig12_WVtrends.rlat,fig12_WVtrends.plays100,fig12_WVtrends.umbc,fig12_WVtrends.airsL3,fig12_WVtrends.climcaps,fig12_WVtrends.merra2,fig12_WVtrends.era5,iFig,fig12_WVtrends.plotoptions2x1x2);
  %% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig12.mat fig12_WVtrends
  
  %% from compare_colWV_trends_Day_vs_Night.m
  figure(92)
  figure(192); clf
  plot(fig10_colwvtrend.rlat,fig10_colwvtrend.umbc,'k',fig10_colwvtrend.rlat,fig10_colwvtrend.airsL3,'b',fig10_colwvtrend.rlat,fig10_colwvtrend.climcaps,'g',...
       fig10_colwvtrend.rlat,fig10_colwvtrend.era5,'r',fig10_colwvtrend.rlat,fig10_colwvtrend.merra2,'m','linewidth',2);
  ylim([-1 +2]*0.05);
  xlim([-1 +1]*90); 
  plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
  xlabel('Latitude [deg]'); ylabel('d mmw/dt [mm/yr]');
  fig10_colwvtrend.comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_trends_Day_vs_Night.m  --> compare_colWV_trends_Day_vs_Night.m';
  %% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig10.mat fig10_colwvtrend
end
