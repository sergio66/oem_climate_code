%% cp /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m compare_05_10_15_20_spectral_rates.m

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/plotutils

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
rlat = cos(rlat*pi/180);

bah = a05.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a10.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a15.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a20.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = a05.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a10.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a15.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a20.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = a05.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a10.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a15.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a20.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = a05.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a10.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a15.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a20.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = a05.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a10.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a15.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = a20.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

figure(1); clf;
plot(h.vchan,mean_BT); xlabel('Wavenumber cm-1'); ylabel('BT(K)'); 
xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);

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
  q_mean_BT(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q02.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q03.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q04.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q05.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_BT(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = q01.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q02.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q03.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q04.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q05.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_desc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = q01.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q02.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q03.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q04.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q05.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_asc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = q01.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q02.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q03.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q04.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q05.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_desc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = q01.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q02.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q03.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q04.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = q05.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  q_mean_err_asc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

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
