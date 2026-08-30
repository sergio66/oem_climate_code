%{
     -rw-rw-r-- 1 sergio pi_sergio  5223 Jan 22  2026 set_zeroT_nofit.m
     -rw-rw-r-- 1 sergio pi_sergio  5566 Jan 22  2026 set_zeroWV_nofit.m
     -rw-rw-r-- 1 sergio pi_sergio  2577 Jan 22  2026 show_unc.m
     -rw-rw-r-- 1 sergio pi_sergio 16997 Jan 22  2026 strow_override_defaults_latbins_AIRS_fewlays.m
     -rw-rw-r-- 1 sergio pi_sergio  1159 Jan 22  2026 trop_index.m
     -rw-rw-r-- 1 sergio pi_sergio 19346 Jan 22  2026 set_driver_rateset_datafile.m
     -rw-rw-r-- 1 sergio pi_sergio 14498 Jan 22  2026 set_the_AMSU_jcobians.m
     -rw-rw-r-- 1 sergio pi_sergio 14651 Jan 22  2026 set_the_jacobians.m
     -rw-rw-r-- 1 sergio pi_sergio  5566 Jan 22  2026 set_zeroO3_nofit.m
     -rw-rw-r-- 1 sergio pi_sergio 16079 Jan 22  2026 set_apriori_ERA5_MERRA2_or_AIRSL3_MLS_geophysical.m
     -rw-rw-r-- 1 sergio pi_sergio  8911 Jan 22  2026 set_driver_jacfile.m
     -rw-rw-r-- 1 sergio pi_sergio  2348 Jan 22  2026 set_Tz_O3z_noFit.m
     -rw-rw-r-- 1 sergio pi_sergio   862 Jan 22  2026 nc_rates.m
     -rw-rw-r-- 1 sergio pi_sergio  2163 Jan 22  2026 precipitation_vs_skt_changes.m
     -rw-rw-r-- 1 sergio pi_sergio 11584 Jan 22  2026 see_clust_put_together_jacs_cldERA5_2022.m
     -rw-rw-r-- 1 sergio pi_sergio 18539 Jan 22  2026 set_CO2_CH4_N2O_ESRL.m
     -rw-rw-r-- 1 sergio pi_sergio  5687 Jan 22  2026 get_co2_n2o_ch4_for_strow_override.m
     -rw-rw-r-- 1 sergio pi_sergio  1960 Jan 22  2026 get_ESRL_TRACE_GAS_2002_2022.m
     -rw-rw-r-- 1 sergio pi_sergio  6701 Jan 22  2026 get_jac_fast.m
     -rw-rw-r-- 1 sergio pi_sergio  8750 Jan 22  2026 get_rates.m
     -rw-rw-r-- 1 sergio pi_sergio  5313 Jan 22  2026 jeevanjee_PBL_deltaRH_Ts.m
     -rw-rw-r-- 1 sergio pi_sergio   744 Jan 22  2026 combinejaclays.m
OOPS -rw-rw-r-- 1 sergio pi_sergio 27470 Jan 22  2026 do_the_cov_set_numbers.m  
     -rw-rw-r-- 1 sergio pi_sergio  1239 Jan 22  2026 estimate_fracWV_for_deltaRH_zero.m
     -rw-rw-r-- 1 sergio pi_sergio  6595 Jan 22  2026 find_the_oem_channels.m
     -rw-rw-r-- 1 sergio pi_sergio  2317 Jan 22  2026 find_wgtA_wgtB.m
     -rw-rw-r-- 1 sergio pi_sergio 54539 Jan 22  2026 clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
OOPS -rw-rw-r-- 1 sergio pi_sergio 36093 Jan 22  2026 build_cov_matrices.m
     -rw-rw-r-- 1 sergio pi_sergio  7251 Jan 22  2026 change_important_topts_settings.m
     -rw-rw-r-- 1 sergio pi_sergio  1495 Jan 22  2026 check_driver.m
     -rw-rw-r-- 1 sergio pi_sergio  7150 Jan 22  2026 check_settings.m
     -rw-rw-r-- 1 sergio pi_sergio  7259 Jan 22  2026 choose_goodchans_from_2645.m
     %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ls -1 /home/sergio/git/JGR_July2025/CODE
list_files = {...
'build_cov_matrices.m',...
'change_important_topts_settings.m',...
'check_driver.m',...
'check_settings.m',...
'choose_goodchans_from_2645.m',...
'clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m',...
'combinejaclays.m',...
'do_the_cov_set_numbers.m',...
'estimate_fracWV_for_deltaRH_zero.m',...
'find_the_oem_channels.m',...
'find_wgtA_wgtB.m',...
'get_co2_n2o_ch4_for_strow_override.m',...
'get_ESRL_TRACE_GAS_2002_2022.m',...
'get_jac_fast.m',...
'get_rates.m',...
'jeevanjee_PBL_deltaRH_Ts.m',...
'nc_rates.m',...
'precipitation_vs_skt_changes.m',...
'see_clust_put_together_jacs_cldERA5_2022.m',...
'set_apriori_ERA5_MERRA2_or_AIRSL3_MLS_geophysical.m',...
'set_CO2_CH4_N2O_ESRL.m',...
'set_driver_jacfile.m',...
'set_driver_rateset_datafile.m',...
'set_the_AMSU_jcobians.m',...
'set_the_jacobians.m',...
'set_Tz_O3z_noFit.m',...
'set_zeroO3_nofit.m',...
'set_zeroT_nofit.m',...
'set_zeroWV_nofit.m',...
'show_unc.m',...
'strow_override_defaults_latbins_AIRS_fewlays.m',...
'trop_index.m'}

for ii = 1 : length(list_files)
  clc
  differ = ['!diff ' list_files{ii} ' /home/sergio/git/JGR_July2025/CODE/.'];
  eval(differ)
  disp(' ')
  fprintf(1,'>>>> just diffed %s ret to continue if big changes, see OOPS above\n',list_files{ii});
  pause
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat

addpath /home/sergio/git/JGR_July2025/PLOTTER
make_jgr_fig3

figure(5); clf; plot(h.vchan,squeeze(nanmean(b_desc,1)))
figure(5); clf; plot(h.vchan,squeeze(nanmean(squeeze(nanmean(b_desc,1)),1)))
figure(6); clf; plot(h.vchan,squeeze(nanmean(squeeze(nanmean(b_err_desc,1)),1)))

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
make_jgr_fig5
ax = axis; cx = caxis;

yy = unique(Y);
figure(7); clf; pcolor(h.vchan,yy,squeeze(nanmean(b_asc,1)))
figure(7); clf; pcolor(h.vchan,yy,squeeze(nanmean(b_desc,1)))
shading flat; colorbar; colormap(usa2); caxis(cx); axis(ax);

figure(8); clf;
pcolor(h.vchan,yy,data_fig5b.l1c_trendD - squeeze(nanmean(b_desc,1)))
shading flat; colorbar; colormap(usa2); caxis(cx/50); axis(ax);
