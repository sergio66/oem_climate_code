iNumYears = 20;

junk = load(['ceres_trends_' num2str(iNumYears,'%02d') 'year.mat']);
  ceres_trend = junk.ceres_trend;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'era5_spectral_olr');
  era5_spectral_olr = junk.era5_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;

load latB64.mat
figure(1); clf
figure(1); plot(ceres_trend.trend_lat,ceres_trend.trend_lw,ceres_trend.trend_lat,ceres_trend.trend_lw_clr,'linewidth',6);
miaow = era5_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - era5_spectral_olr.olr0_ecRad.clr;;
  figure(1); hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
miaow = merra2_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - merra2_spectral_olr.olr0_ecRad.clr;;
  figure(1); hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
miaow = umbc_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr;;
  figure(1); hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'x-','linewidth',4);
miaow = airsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
  figure(1); hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
miaow = climcapsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - climcapsL3_spectral_olr.olr0_ecRad.clr;;
  figure(1); hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
hold off; plotaxis2; hl = legend('CERES allsky','CERES clrsky','ERA5','MERRA2','THIS WORK','AIRSL3','CLIMAPS L3','location','best','fontsize',10); 
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
axax = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat')
miaow = airsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
  figure(2); plot(ceres_trend.trend_lat,ceres_trend.trend_lw,ceres_trend.trend_lat,ceres_trend.trend_lw_clr,...
                 meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'kx-',...
                 meanvaluebin(latB2),nanmean(reshape(axax.thestats64x72_other.clrolrrate,72,64),1),meanvaluebin(latB2),nanmean(reshape(axax.thestats64x72_other.olrrate,72,64),1),'linewidth',2);
plotaxis2; legend('CERES allsky','CERES clrsky','Sergio AIRS L3','Joel AIRS L3 clr','Joel AIRS L3 allsky','location','best','fontsize',10);
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

%hmm I did not do CLIMCAPS OLR
%/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat
%axax = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat');
%axax = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
