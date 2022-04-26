moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','b_err_desc');
b_err_desc = moonoise.b_err_desc; clear moonoise;
b_err_desc = permute(b_err_desc,[3 1 2]);
b_err_desc = reshape(b_err_desc,2645,72*64);

moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','b_desc');
b_desc = moonoise.b_desc; clear moonoise;
b_desc = permute(b_desc,[3 1 2]);
b_desc = reshape(b_desc,2645,72*64);

moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','h');
f = moonoise.h.vchan; clear moonoise;

addpath ../FIND_NWP_MODEL_TRENDS/
driver_get_the_model_trends

settings.iIgnoreChans_CH4 = -1;
settings.iIgnoreChans_N2O = -1;
settings.iIgnoreChans_SO2 = -1;
%chanset = jacobian.chanset;
chanset = 1 : 2645;
[hMean17years,ha,pMean17years,pa]     = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
plotopt.iUpperWavenumLimit = 1620;
plotopt.rlon = pMean17years.rlon;
plotopt.rlat = pMean17years.rlat;
plotopt.scatterORaslmap = -1;

miaow = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_globalSSTfeedback.mat','fits');
fits = miaow.fits; clear miaow;

disp('hit ret to compare fits vs obs'); pause
[raaBadFov,indBadFov] = plot_spectral_region_chisqr(b_desc(chanset,:),0*b_desc(chanset,:),0*b_desc(chanset,:),fits(chanset,:),f(chanset,:),b_err_desc(chanset,:),-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*0.25); end; figure(20); caxis([0 1]*0.1)

disp('hit ret to compare ERA5/MERRA2 vs obs'); pause
[raaBadFov,indBadFov] = plot_spectral_region_chisqr(b_desc(chanset,:),0*b_desc(chanset,:),0*b_desc(chanset,:),era5.era5_spectral_rates(chanset,:),f(chanset,:),b_err_desc(chanset,:),-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*10); end; figure(20); caxis([0 1]*2)

disp('hit ret to compare AIRSL3/CLIMCAPSL3 vs obs'); pause
[raaBadFov,indBadFov] = plot_spectral_region_chisqr(b_desc(chanset,:),0*b_desc(chanset,:),0*b_desc(chanset,:),airsL3.airsL3_spectral_rates(chanset,:),f(chanset,:),b_err_desc(chanset,:),-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*10); end; figure(20); caxis([0 1]*2)

disp('hit ret to compare CMIP6/AMIP6 vs obs'); pause
[raaBadFov,indBadFov] = plot_spectral_region_chisqr(b_desc(chanset,:),0*b_desc(chanset,:),0*b_desc(chanset,:),cmip6.cmip6_spectral_rates(chanset,:),f(chanset,:),b_err_desc(chanset,:),-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*10); end; figure(20); caxis([0 1]*2)
