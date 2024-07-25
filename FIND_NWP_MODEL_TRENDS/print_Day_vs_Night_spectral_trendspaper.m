iPrint = +1;
iPrint = -1;

dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN/';
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN_Temp/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% if iPrint > 0
%%   addpath /asl/matlib/plotutils
%% 
%%   fprintf(1,'will be saving plots/figs etc to %s \n',dir0)
%% 
%%   figure(1); aslprint([dir0 'BT1231_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
%%   figure(2); aslprint([dir0 'BT1226_BT1231_zonaltrends_mean_obs_chirp_era5_airsL3_climcapsL3.pdf']);
%%   figure(3); aslprint([dir0 'BT1226_BT1231_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
%%   figure(4); aslprint([dir0 'BT1419_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
%%   figure(5); aslprint([dir0 'BT1519_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
%% 
%%   figure(6); aslprint([dir0 'land_ocean_spectra_all_tropics_midlats_polar_nighttrends.pdf']);
%%   figure(7); aslprint([dir0 'land_ocean_spectra_all_tropics_midlats_polar_daytrends.pdf']);
%%   figure(8); aslprint([dir0 'land_ocean_spectra_all_tropics_midlats_polar_avgtrends.pdf']);  %% see also Figs/land_ocean_spectra_all_tropics_midlats_polar.pdf; do_correlations_of_spectra
%%   figure(9); aslprint([dir0 'land_ocean_spectra_all_tropics_midlats_polar_difftrends.pdf']);
%% end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iPrint > 0
  addpath /asl/matlib/plotutils

  fprintf(1,'will be saving plots/figs etc to %s \n',dir0)

  figure(1); sergioprintfig([dir0 'BT1231_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
  figure(2); set(gca,'fontsize',12); sergioprintfig([dir0 'BT1226_BT1231_zonaltrends_mean_obs_chirp_era5_airsL3_climcapsL3.pdf']);
  figure(3); sergioprintfig([dir0 'BT1226_BT1231_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
  figure(4); sergioprintfig([dir0 'BT1419_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);
  figure(5); sergioprintfig([dir0 'BT1519_trends_5panels_D_vs_N_obs_chirp_era5_airsL3_climcapsL3.pdf']);

  figure(6); sergioprintfig([dir0 'land_ocean_spectra_all_tropics_midlats_polar_nighttrends.pdf']);
  figure(7); sergioprintfig([dir0 'land_ocean_spectra_all_tropics_midlats_polar_daytrends.pdf']);
  figure(8); sergioprintfig([dir0 'land_ocean_spectra_all_tropics_midlats_polar_avgtrends.pdf']);  %% see also Figs/land_ocean_spectra_all_tropics_midlats_polar.pdf; do_correlations_of_spectra
  figure(9); sergioprintfig([dir0 'land_ocean_spectra_all_tropics_midlats_polar_difftrends.pdf']);
end

