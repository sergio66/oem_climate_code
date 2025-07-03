%% cd  /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/
%% do_the_plots_ecmwf_or_era_16days_tile_v2
%% 
%% stop at Fig25 and save
%% save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig1.mat data_fig1

jett = jet; jett(1,:) = 1;

load fig1.mat
figure(1); clf
pcolor(data_fig1.dbtt,data_fig1.latbins,data_fig1.histc); colormap(jett); colorbar; shading interp
hold on;
plot(data_fig1.Q0,data_fig1.latbins,'c','linewidth',2);
plot(data_fig1.mean,data_fig1.latbins,'rx-','linewidth',2);
plot(data_fig1.Q50,data_fig1.latbins,'linewidth',2);
plot(data_fig1.Q90,data_fig1.latbins,'color',[1 1 1]*0.5,'linewidth',2);
plot(data_fig1.Q1,data_fig1.latbins,'b','linewidth',2);
hold off
xlim([210 310])
%wah = [['    ']; ['min ']; ['mean']; num2str(quants(IA(1))','%4.2f'); ['0.90']; ['max ']];
wah = [['    ']; ['min ']; ['mean']; ['0.50']; ['0.90']; ['max ']];
hl = legend(wah,'location','west','fontsize',8);
set(gca,'ydir','normal')
ylim([-90 +90])
xlabel('BT 1231 (K)'); ylabel('Latitude')
