addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/plotutils

% cp /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m contrast_obsrates_05_10_15_20.m

if ~exist('a05')
  dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/';
  a05 = load([dir0 'iType_10_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  a10 = load([dir0 'iType_11_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  a15 = load([dir0 'iType_12_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  a20 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  
  loader = ['load ' dir0 'h2645structure.mat'];
  eval(loader)

  [h,ha,p,pa] = rtpread([dir0 'summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp']);
end
  
rlat = ones(2645,1)*p.rlat;
cosrlat = cos(rlat*pi/180);

a05.b_asc_zonal = squeeze(nanmean(a05.b_asc,1));  a05.b_desc_zonal = squeeze(nanmean(a05.b_desc,1)); 
a10.b_asc_zonal = squeeze(nanmean(a10.b_asc,1));  a10.b_desc_zonal = squeeze(nanmean(a10.b_desc,1)); 
a15.b_asc_zonal = squeeze(nanmean(a15.b_asc,1));  a15.b_desc_zonal = squeeze(nanmean(a15.b_desc,1)); 
a20.b_asc_zonal = squeeze(nanmean(a20.b_asc,1));  a20.b_desc_zonal = squeeze(nanmean(a20.b_desc,1)); 

bah = a05.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a10.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a15.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a20.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = a05.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a10.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a15.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a20.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = a05.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a10.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a15.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a20.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = a05.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a10.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a15.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a20.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = a05.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a10.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a15.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = a20.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

figure(1); clf;
plot(h.vchan,mean_BT); xlabel('Wavenumber cm-1'); ylabel('BT(K)'); 
xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);

figure(2); clf
plot(h.vchan,mean_desc_trend); xlabel('Wavenumber cm-1'); ylabel('desc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);
  ylim([-1 +1]*0.05)

figure(2); clf
subplot(211); plot(h.vchan,mean_desc_trend); xlabel('Wavenumber cm-1'); ylabel('desc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);
subplot(212); plot(h.vchan,mean_err_desc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);

figure(3); clf
subplot(211); plot(h.vchan,mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('asc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);
subplot(212); plot(h.vchan,mean_err_asc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('q01')
  q01 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat']);
  q02 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat']);
  q03 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  q04 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat']);
  q05 = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat']);
end

bah = q01.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q02.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q03.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q04.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q05.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,5) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = q01.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q02.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q03.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q04.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q05.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,5) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = q01.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q02.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q03.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q04.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q05.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,5) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = q01.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q02.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q03.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q04.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q05.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,5) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

bah = q01.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,1) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q02.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,2) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q03.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,3) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q04.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,4) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);
bah = q05.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,5) = nanmean(cosrlat.*bah,2)./nanmean(cosrlat,2);

figure(4); clf;
plot(h.vchan,q_mean_BT); xlabel('Wavenumber cm-1'); ylabel('BT(K)'); 
xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

figure(5); clf
subplot(211); plot(h.vchan,q_mean_desc_trend); xlabel('Wavenumber cm-1'); ylabel('desc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);
subplot(212); plot(h.vchan,q_mean_err_desc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

figure(6); clf
subplot(211); plot(h.vchan,q_mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('asc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);
subplot(212); plot(h.vchan,q_mean_err_asc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

figure(7); clf
plot(h.vchan,q_mean_desc_trend-q_mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('desc-asc d(BT)/dt (K/yr)'); 
xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

figure(8); clf
plot(h.vchan,q_mean_desc_trend(:,1)*ones(1,4)-q_mean_desc_trend(:,2:5)); xlabel('Wavenumber cm-1'); ylabel('differences in d(BT)/dt (K/yr)'); 
xlim([640 1640]); plotaxis2; hl = legend('Q50-Q80','Q50-Q90','Q50-Q95','Q50-Q97','location','best','fontsize',8);

figure(9); clf
plot(h.vchan,q_mean_desc_trend(:,3)*ones(1,4)-q_mean_desc_trend(:,[1 2 4 5])); xlabel('Wavenumber cm-1'); ylabel('differences in d(BT)/dt (K/yr)'); 
xlim([640 1640]); plotaxis2; hl = legend('Q90-Q50','Q90-Q80','Q90-Q95','Q90-Q97','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 05;  %% use CarbonTracker CO2 trends
  results05yrs = load(strUMBC,'results','resultsO3','resultsWV','resultsT','pavg');
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2v2.mat';    iNumYears = 10;  %% use CarbonTracker CO2 trends
  results10yrs = load(strUMBC,'results','resultsO3','resultsWV','resultsT','pavg','rlat');
a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 15;  %% use CarbonTracker CO2 trends
  results15yrs = load(strUMBC,'results','resultsO3','resultsWV','resultsT','pavg');
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 20;  %% use CarbonTracker CO2 trends ****
  results20yrs = load(strUMBC,'results','resultsO3','resultsWV','resultsT','pavg');

z11 = squeeze(nanmean(reshape(results05yrs.resultsWV,72,64,49),1))';
z12 = squeeze(nanmean(reshape(results10yrs.resultsWV,72,64,49),1))';
z21 = squeeze(nanmean(reshape(results15yrs.resultsWV,72,64,49),1))';
z22 = squeeze(nanmean(reshape(results20yrs.resultsWV,72,64,49),1))';
plotoptions.str11 = '05 yrs'; plotoptions.str12 = '10 yrs'; plotoptions.str21 = '15 yrs'; plotoptions.str22 = '20 yrs';
plotoptions.xstr = 'Latitude'; plotoptions.ystr = 'Pressure (mb)';
plotoptions.cx = [-1 +1]*0.015; plotoptions.maintitle = 'd/dt WVfrac'; plotoptions.cmap = usa2; plotoptions.yReverseDir = +1; plotoptions.yLinearOrLog = +1; plotoptions.yLimits = [100 1000]; plotoptions.xLimits = [-90 +90]; 
profile_plots_2x2tiledlayout(results10yrs.rlat,results10yrs.pavg,z11,z12,z21,z22,10,plotoptions);

z11 = squeeze(nanmean(reshape(results05yrs.resultsT,72,64,49),1))';
z12 = squeeze(nanmean(reshape(results10yrs.resultsT,72,64,49),1))';
z21 = squeeze(nanmean(reshape(results15yrs.resultsT,72,64,49),1))';
z22 = squeeze(nanmean(reshape(results20yrs.resultsT,72,64,49),1))';
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = usa2; plotoptions.yReverseDir = +1; plotoptions.yLinearOrLog = -1; plotoptions.yLimits = [10 1000]; plotoptions.xLimits = [-90 +90]; 
profile_plots_2x2tiledlayout(results10yrs.rlat,results10yrs.pavg,z11,z12,z21,z22,11,plotoptions);

z11 = squeeze(nanmean(reshape(results05yrs.resultsO3,72,64,49),1))';
z12 = squeeze(nanmean(reshape(results10yrs.resultsO3,72,64,49),1))';
z21 = squeeze(nanmean(reshape(results15yrs.resultsO3,72,64,49),1))';
z22 = squeeze(nanmean(reshape(results20yrs.resultsO3,72,64,49),1))';
plotoptions.cx = [-1 +1]*0.005; plotoptions.maintitle = 'd/dt O3frac'; plotoptions.cmap = usa2; plotoptions.yReverseDir = +1; plotoptions.yLinearOrLog = -1; plotoptions.yLimits = [1 100]; plotoptions.xLimits = [-90 +90]; 
profile_plots_2x2tiledlayout(results10yrs.rlat,results10yrs.pavg,z11,z12,z21,z22,12,plotoptions);

z11 = reshape(a05.b_desc,4608,2645); z11 = z11(:,1520); 
z12 = reshape(a10.b_desc,4608,2645); z12 = z12(:,1520); 
z21 = reshape(a15.b_desc,4608,2645); z21 = z21(:,1520); 
z22 = reshape(a20.b_desc,4608,2645); z22 = z22(:,1520); 
plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
plotoptions.cx = [-1 +1]*0.2; plotoptions.maintitle = 'dBT1231/dt'; plotoptions.cmap = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
aslmap_2x2tiledlayout(z11,z12,z21,z22,13,plotoptions);

z11 = results05yrs.results(:,6); z12 = results10yrs.results(:,6); z21 = results15yrs.results(:,6); z22 = results20yrs.results(:,6); 
plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
plotoptions.cx = [-1 +1]*0.2; plotoptions.maintitle = 'dSKT/dt'; plotoptions.cmap = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
aslmap_2x2tiledlayout(z11,z12,z21,z22,14,plotoptions);

plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Wavenumber (cm-1)';
z11 = a05.b_desc_zonal; z12 = a10.b_desc_zonal; z21 = a15.b_desc_zonal; z22 = a20.b_desc_zonal;
plotoptions.cx = [-1 +1]*0.2; plotoptions.maintitle = 'dBT/dt'; plotoptions.cmap = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
profile_plots_2x2tiledlayout(a05.h.vchan,results10yrs.rlat,z11,z12,z21,z22,15,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = -1;
if iPrint > 0
  figure(1); xlim([640 1640]); aslprint('../Figs/mean_BT_obs_desc_05_10_15_20_years.pdf');
  figure(2); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_desc_05_10_15_20_years.pdf');
  figure(3); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_asc_05_10_15_20_years.pdf');

  figure(4); xlim([640 1640]); aslprint('../Figs/mean_BT_obs_desc_20_years_Q50_Q80_Q90_Q95_Q97.pdf');
  figure(5); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_desc_20_years_Q50_Q80_Q90_Q95_Q97.pdf');
  figure(6); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_asc_20_years_Q50_Q80_Q90_Q95_Q97.pdf');

  figure(8); xlim([640 1640]); axis([640 1640 -0.005 +0.005]); aslprint('../Figs/trend_dBTdt_dsc_20_years_Q50_versus_Q80_Q90_Q95_Q97.pdf');
end
