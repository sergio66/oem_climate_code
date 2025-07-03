%% cd  /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/
%% do_the_plots_ecmwf_or_era_16days_tile_v2
%%   plot_QX_vs_uniformclear_filter
%%
%% stop at Fig 08,09,11 and save
%% save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig2.mat data_fig2

%{
cp -a /umbc/xfs2/strow/asl/s1/sergio/home/git/MATLABCODE_Git/PLOTTER/simplemap.m .
cp -a /umbc/xfs2/strow/asl/s1/sergio/home/git/MATLABCODE_Git/PLOTMISC/shadedErrorBar.m .
%}

load fig2.mat

  figure(1); clf; simplemap(data_fig2.uniform_clear_180x380.bt1231); colormap jet
  figure(2); clf; simplemap(data_fig2.tileQ90_180x380.bt1231); colormap jet
  figure(2); cx = caxis;
  figure(1); caxis(cx);

  figure(3);
  plot(data_fig2.rlat180,data_fig2.zonal.uniform_clear_mean,'b',data_fig2.rlat180,data_fig2.zonal.Q90_mean,'r','linewidth',2);
  hold on
  shadedErrorBar(data_fig2.rlat180,data_fig2.zonal.uniform_clear_mean,data_fig2.zonal.uniform_clear_std,'b',0.3);
  shadedErrorBar(data_fig2.rlat180,data_fig2.zonal.Q90_mean,data_fig2.zonal.Q90_std,'r',0.3);
  hold off
  hl = legend('QX filter','clear filter','location','south');
  xlim([-90 +90]); ylabel('BT1231 observed (K)'); xlabel('Latitude [deg]');

  figure(1); set(gca,'fontsize',18)
  figure(2); set(gca,'fontsize',18)
  figure(3); set(gca,'fontsize',18)

