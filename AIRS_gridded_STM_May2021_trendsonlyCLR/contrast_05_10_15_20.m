addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/plotutils

%% this code is originally in /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m
%% see also driver_contrast_obsrates_05_10_15_20.m

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
  ylabel('d(BT)/dt (K/yr)');
subplot(212); plot(h.vchan,mean_err_desc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);

figure(3); clf
subplot(211); plot(h.vchan,mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('asc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('05','10','15','20','location','best','fontsize',8);
  ylabel('d(BT)/dt (K/yr)');
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
  ylabel('d(BT)/dt (K/yr)');
subplot(212); plot(h.vchan,q_mean_err_desc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

figure(6); clf
subplot(211); plot(h.vchan,q_mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('asc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);
  ylabel('d(BT)/dt (K/yr)');
subplot(212); plot(h.vchan,q_mean_err_asc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

figure(10); clf
berrD = squeeze(nanmean(a20.b_err_desc,1));
berrA = squeeze(nanmean(a20.b_err_asc,1));
berr = 0.5*(berrA+berrD);
rlat = nanmean(reshape(p.rlat,72,64),1);

pcolor(h.vchan,rlat,berr);  %% (A+D)/2
pcolor(h.vchan,rlat,berrA); %% A
pcolor(h.vchan,rlat,berrD); %% D
 shading flat; colorbar; colormap(jet);title('Obs Rates Unc'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])

figure(11); clf
bsigD = squeeze(nanmean(a20.b_desc,1));
bsigA = squeeze(nanmean(a20.b_asc,1));
bsig = 0.5*(bsigA+bsigD);
rlat = nanmean(reshape(p.rlat,72,64),1);

pcolor(h.vchan,rlat,bsig); 
pcolor(h.vchan,rlat,bsigA); 
pcolor(h.vchan,rlat,bsigD); 
 shading flat; colorbar; colormap(usa2);title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640]); caxis([-1 +1]*0.1)

figure(12); clf;
pcolor(h.vchan,rlat,abs(bsig./berr)); 
pcolor(h.vchan,rlat,abs(bsigA./berrA)); 
pcolor(h.vchan,rlat,abs(bsigD./berrD)); 
 shading flat; colorbar; colormap(jet);title('Obs Rates/Unc = S/N'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640]); caxis([0 +1]*5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now change fig5 so it uses tiles from driver_compare_cldy_Qtrends_2002_2022.m

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat
moo01 = b_desc;
unc01 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat
moo02 = b_desc;
unc02 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat
moo03 = b_desc;
unc03 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat
moo04 = b_desc;
unc04 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat
moo05 = b_desc;
unc05 = b_err_desc;

load h2645structure.mat

zoo01 = permute(moo01,[3 1 2]); zoo01 = reshape(zoo01,2645,72*64);
zoo02 = permute(moo02,[3 1 2]); zoo02 = reshape(zoo02,2645,72*64);
zoo03 = permute(moo03,[3 1 2]); zoo03 = reshape(zoo03,2645,72*64);
zoo04 = permute(moo04,[3 1 2]); zoo04 = reshape(zoo04,2645,72*64);
zoo05 = permute(moo05,[3 1 2]); zoo05 = reshape(zoo05,2645,72*64);

uoo01 = permute(unc01,[3 1 2]); uoo01 = reshape(uoo01,2645,72*64);
uoo02 = permute(unc02,[3 1 2]); uoo02 = reshape(uoo02,2645,72*64);
uoo03 = permute(unc03,[3 1 2]); uoo03 = reshape(uoo03,2645,72*64);
uoo04 = permute(unc04,[3 1 2]); uoo04 = reshape(uoo04,2645,72*64);
uoo05 = permute(unc05,[3 1 2]); uoo05 = reshape(uoo05,2645,72*64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos moo* unc* zoo* uoo*
do_XX_YY_from_X_Y
cosYY = ones(2645,1)*cos(YY*pi/180);
czoo01 = nansum(cosYY.*zoo01,2)./nansum(cosYY,2);  cuoo01 = nansum(cosYY.*uoo01,2)./nansum(cosYY,2); 
czoo02 = nansum(cosYY.*zoo02,2)./nansum(cosYY,2);  cuoo02 = nansum(cosYY.*uoo02,2)./nansum(cosYY,2); 
czoo03 = nansum(cosYY.*zoo03,2)./nansum(cosYY,2);  cuoo03 = nansum(cosYY.*uoo03,2)./nansum(cosYY,2); 
czoo04 = nansum(cosYY.*zoo04,2)./nansum(cosYY,2);  cuoo04 = nansum(cosYY.*uoo04,2)./nansum(cosYY,2); 
czoo05 = nansum(cosYY.*zoo05,2)./nansum(cosYY,2);  cuoo05 = nansum(cosYY.*uoo05,2)./nansum(cosYY,2); 

figure(20); clf
subplot(211); plot(h.vchan,nanmean(zoo01,2),h.vchan,nanmean(zoo02,2),h.vchan,nanmean(zoo03,2),h.vchan,nanmean(zoo04,2),h.vchan,nanmean(zoo05,2))
subplot(212); plot(h.vchan,nanmean(uoo01,2),h.vchan,nanmean(uoo02,2),h.vchan,nanmean(uoo03,2),h.vchan,nanmean(uoo04,2),h.vchan,nanmean(uoo05,2))

figure(20); clf
ta = tiledlayout(2,1);
tafov(1) = nexttile; plot(h.vchan,czoo01,h.vchan,czoo02,h.vchan,czoo03,h.vchan,czoo04,h.vchan,czoo05); plotaxis2;
  xlim([640 1620]); ylabel('d(BT)/dt [K yr^{-1}]');  
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'location','south'); 
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.45 0.25 0.25])
  set(gca,'fontsize',14)
tafov(2) = nexttile; plot(h.vchan,cuoo01,h.vchan,cuoo02,h.vchan,cuoo03,h.vchan,cuoo04,h.vchan,cuoo05); plotaxis2;
  xlim([640 1620]); 
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',8,'location','southeast'); 
  legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.4625 0.2 0.2])
  xlabel('Wavenumber [cm^{-1}]'); ylabel('d(BT)/dt [K yr^{-1}]');
  ylim([0 0.031])
  set(gca,'fontsize',14)

ta.Padding = 'none';
ta.TileSpacing = 'compact';
ta.Padding = 'compact';
ta.TileSpacing = 'tight';
tafov(1).XTickLabel = '';

%{
data_fig5b.vchan      = h.vchan;
data_fig5b.rlat       = rlat;
data_fig5b.l1c_trendD = bsigD;
data_fig5b.l1c_uncD   = berrD;
save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig5b.mat data_fig5b

data_fig3a.vchan = h.vchan;
data_fig3a.q_mean_BT = q_mean_BT;

data_fig3b.vchan = h.vchan;
data_fig3b.q01.trend = czoo01; data_fig3b.q01.unc = cuoo01; 
data_fig3b.q02.trend = czoo02; data_fig3b.q02.unc = cuoo02; 
data_fig3b.q03.trend = czoo03; data_fig3b.q03.unc = cuoo03; 
data_fig3b.q04.trend = czoo04; data_fig3b.q04.unc = cuoo04; 
data_fig3b.q05.trend = czoo05; data_fig3b.q05.unc = cuoo05; 
save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig3.mat data_fig3b data_fig3a

%% sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/trend_dBTdt_desc_20_years_Q50_Q80_Q90_Q95_Q97_bigfont')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


  %% figure(10); xlim([640 1640]);  aslprint('../Figs/obs_trend_unc_20_years.pdf');  %% notice this is night
end
