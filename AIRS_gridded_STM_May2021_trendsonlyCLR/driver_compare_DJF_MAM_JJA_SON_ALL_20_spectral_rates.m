%% cp /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m compare_05_10_15_20_spectral_rates.m

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps

if ~exist('aALL')
  dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/';

  aDJF = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03_DJF.mat']);
  aMAM = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03_MAM.mat']);
  aJJA = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03_JJA.mat']);
  aSON = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03_SON.mat']);

  aALL = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  
  loader = ['load ' dir0 'h2645structure.mat'];
  eval(loader)

  [h,ha,p,pa] = rtpread([dir0 'summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp']);
end
  
rlat = ones(2645,1)*p.rlat;
rlat = cos(rlat*pi/180);

bah = aDJF.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aMAM.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aJJA.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aSON.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aALL.mean_BT; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_BT(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = aDJF.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aMAM.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aJJA.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aSON.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aALL.b_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_desc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = aDJF.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aMAM.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aJJA.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aSON.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aALL.b_err_desc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_desc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = aDJF.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aMAM.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aJJA.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aSON.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aALL.b_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_asc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

bah = aDJF.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,1) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aMAM.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,2) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aJJA.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,3) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aSON.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,4) = nanmean(rlat.*bah,2)./nanmean(rlat,2);
bah = aALL.b_err_asc; bah = permute(bah,[3 1 2]); bah = reshape(bah,2645,4608); 
  mean_err_asc_trend(:,5) = nanmean(rlat.*bah,2)./nanmean(rlat,2);

figure(1); clf;
plot(h.vchan,mean_BT); xlabel('Wavenumber cm-1'); ylabel('BT(K)'); 
xlim([640 1640]); plotaxis2; hl = legend('DJF','MAM','JJA','SON','ALL','location','best','fontsize',8);

figure(2); clf
subplot(211); plot(h.vchan,mean_desc_trend); xlabel('Wavenumber cm-1'); ylabel('desc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('DJF','MAM','JJA','SON','ALL','location','best','fontsize',8);
  ylim([-1 +1]*0.10);
subplot(212); plot(h.vchan,mean_err_desc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('DJF','MAM','JJA','SON','ALL','location','best','fontsize',8);
  ylim([0 +1]*0.10);

figure(3); clf
subplot(211); plot(h.vchan,mean_asc_trend); xlabel('Wavenumber cm-1'); ylabel('asc d(BT)/dt (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('DJF','MAM','JJA','SON','ALL','location','best','fontsize',8);
  ylim([-1 +1]*0.10);
subplot(212); plot(h.vchan,mean_err_asc_trend); xlabel('Wavenumber cm-1'); ylabel('unc (K/yr)'); 
  xlim([640 1640]); plotaxis2; hl = legend('DJF','MAM','JJA','SON','ALL','location','best','fontsize',8);
  ylim([0 +1]*0.10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load latB64.mat
rlon73 = -180 : 5 : +180;  rlat65 = latB2;
rlon = 0.5*(rlon73(1:end-1)+rlon73(2:end));
rlat = 0.5*(rlat65(1:end-1)+rlat65(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

iFig = 3;
iFig = iFig + 1; mooDJF = squeeze(aDJF.b_desc(:,:,1520)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooDJF,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('DJF BT1231 trend')
iFig = iFig + 1; mooMAM = squeeze(aMAM.b_desc(:,:,1520)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooMAM,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('MAM BT1231 trend')
iFig = iFig + 1; mooJJA = squeeze(aJJA.b_desc(:,:,1520)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooJJA,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('JJA BT1231 trend')
iFig = iFig + 1; mooSON = squeeze(aSON.b_desc(:,:,1520)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooSON,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('SON BT1231 trend')
iFig = iFig + 1; mooALL = squeeze(aALL.b_desc(:,:,1520)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooALL,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('ALL BT1231 trend')
aXYZ.b_desc = (aDJF.b_desc + aMAM.b_desc + aJJA.b_desc + aSON.b_desc)/4;
iFig = iFig + 1; mooALL = squeeze(aXYZ.b_desc(:,:,1520)); figure(iFig); simplemap_thickcoast(Y,X,mooALL,5); colormap(usa2); caxis([-1 +1]*0.15); title('ALL<4> BT1231 trend')

%iFig = iFig + 1; 
%  mooALL = squeeze(aXYZ.b_desc(:,:,1520))./squeeze(aALL.b_desc(:,:,1520)); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL-1,5); colormap(usa2); caxis([-1 +1]*0.5); title('1- ALL<4>/ALL BT1231 trend')
%  mooALL = squeeze(aXYZ.b_desc(:,:,1520)) - squeeze(aALL.b_desc(:,:,1520)); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL,5); colormap(usa2); caxis([-1 +1]*0.015); title('1- ALL<4>/ALL BT1231 trend')

iFig = iFig;
junk.color = 'k'; junk.iNorS = +1;
iFig = iFig + 1; mooDJF = squeeze(aDJF.b_desc(:,:,1520)); aslmap_polar(iFig,rlat65,rlon73,smoothn((reshape(mooDJF,72,64)') ,1), [-90 +90],[-180 +180],junk); colormap(usa2); caxis([-1 +1]*0.15); title('DJF BT1231 trend')
iFig = iFig + 1; mooMAM = squeeze(aMAM.b_desc(:,:,1520)); aslmap_polar(iFig,rlat65,rlon73,smoothn((reshape(mooMAM,72,64)') ,1), [-90 +90],[-180 +180],junk); colormap(usa2); caxis([-1 +1]*0.15); title('MAM BT1231 trend')
iFig = iFig + 1; mooJJA = squeeze(aJJA.b_desc(:,:,1520)); aslmap_polar(iFig,rlat65,rlon73,smoothn((reshape(mooJJA,72,64)') ,1), [-90 +90],[-180 +180],junk); colormap(usa2); caxis([-1 +1]*0.15); title('JJA BT1231 trend')
iFig = iFig + 1; mooSON = squeeze(aSON.b_desc(:,:,1520)); aslmap_polar(iFig,rlat65,rlon73,smoothn((reshape(mooSON,72,64)') ,1), [-90 +90],[-180 +180],junk); colormap(usa2); caxis([-1 +1]*0.15); title('SON BT1231 trend')
iFig = iFig + 1; mooALL = squeeze(aALL.b_desc(:,:,1520)); aslmap_polar(iFig,rlat65,rlon73,smoothn((reshape(mooALL,72,64)') ,1), [-90 +90],[-180 +180],junk); colormap(usa2); caxis([-1 +1]*0.15); title('ALL BT1231 trend')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = -1;
if iPrint > 0
  figure(1); xlim([640 1640]); aslprint('../Figs/mean_BT_obs_desc_DJF_MAM_JJA_SON_years.pdf');
  figure(2); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_desc_DJF_MAM_JJA_SON_years.pdf');
  figure(3); xlim([640 1640]); aslprint('../Figs/trend_dBTdt_asc_DJF_MAM_JJA_SON_years.pdf');
end
