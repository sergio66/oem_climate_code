disp('see "simple_get_the_model_trends_do_feedbacks.m" ')
disp('see "compare_OLR_trend.m" ')

addpath /asl/matlib/h4tools/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

if ~exist('era5_spectral_olr')
  [~,~,p,~] = rtpread('summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp');
  %junk = load('olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat','era5_spectral_olr');
  %junk = load('/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_20_unc_factor0.25.mat','era5_spectral_olr');

  iWhich = input('Enter (0/default) ERA5 (1) MERRA2 (2) AIRS L3 (3) CLIMCAPSL3 (4) UMBC : ');
  if length(iWhich) == 0
    iWhich = 0;
  end

  fERA5   = '/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_20.mat';
  fMERRA2 = '/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_20.mat';

  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20.mat'; %% default, but things get written this all the time .. so could be dangerous/inconsistent
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20_defaultGULP_STM_Oct2023.mat';  %% Sounder Oct 2023 STM
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20_SEQN_vers10.mat';              %% trends paper

  if iWhich == 0
    junk = load('fERA5','era5_spectral_olr');  %% this is in compare_OLR_trend.m
    spectral_olr = junk.era5_spectral_olr;
    dataset_name = 'ERA5';
  elseif iWhich == 2
    junk = load('fERA5','airsL3_spectral_olr');  %% this is in compare_OLR_trend.m
    spectral_olr = junk.airsL3_spectral_olr;
    dataset_name = 'AIRS L3';
  elseif iWhich == 1
    junk = load('fMERRA2','merra2_spectral_olr');  %% this is in compare_OLR_trend.m
    spectral_olr = junk.merra2_spectral_olr;
    dataset_name = 'MERRA2';
  elseif iWhich == 3
    junk = load('fMERRA2','climcapsL3_spectral_olr');  %% this is in compare_OLR_trend.m
    spectral_olr = junk.climcapsL3_spectral_olr;
    dataset_name = 'CLIMCAPS L3';
  elseif iWhich == 4
    junk = load('fUMBC','umbc_spectral_olr','strUMBC','iNumYears');  %% this is in compare_OLR_trend.m
    fprintf(1,'UMBC : %2i year dataset = %s \n',junk.iNumYears,junk.strUMBC)
    spectral_olr = junk.umbc_spectral_olr;
    dataset_name = 'UMBC';
  end

  clear junk
  junk = load('latB64.mat');
  lat64 = 0.5*(junk.latB2(1:end-1)+junk.latB2(2:end));
  clear junk
end

olr0_bands    = spectral_olr.olr0_ecRad.bands;
skt_bands     = spectral_olr.skt_ecRad.bands;
planck_bands  = spectral_olr.planck_ecRad.bands;
lapse_bands   = spectral_olr.lapse_ecRad.bands;
o3_bands      = spectral_olr.o3_ecRad.bands;
wv_bands      = spectral_olr.wv_ecRad.bands;
t_co2_bands   = spectral_olr.ptemp_co2_ecRad.bands;               %% compute_feedbacks_generic_ecRad.m, line 352 : T + CO2 only
allpert_bands = spectral_olr.perts9999.atm_skt_ghg_ecRad.bands;   %% compute_feedbacks_generic_ecRad.m, line 156 : the whole shebang : SKT,T,g1,g2,g3,g4,g6

bands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
bandcenter = 0.5*(bands(1:end-1) + bands(2:end));

figure(1); clf
plot(bandcenter,wv_bands.clr - olr0_bands.clr);
plot(bandcenter,nanmean(skt_bands.clr'-olr0_bands.clr'),bandcenter,nanmean(wv_bands.clr'-olr0_bands.clr'),bandcenter,nanmean(planck_bands.clr'-olr0_bands.clr'),bandcenter,nanmean(lapse_bands.clr'-olr0_bands.clr'),'linewidth',2);
plotaxis2; hl = legend('SKT','WV','PLANCK','LAPSE','location','best','fontsize',10); title(dataset_name)

figure(2); clf
wahSKT = squeeze(nanmean(reshape(skt_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahSKT'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (SKT)']); colormap jet
cx = caxis; 
if min(cx) >= 0
  jett = jet; jett(1,:) = 1;
  colormap(jett);
else
  colormap(usa2);
  cx = abs(cx);
  cx = max(cx);
  caxis([-cx +cx]);
end

figure(3); clf
wahWV = squeeze(nanmean(reshape(wv_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahWV'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (WV)']); colormap jet
caxis([-1 +1]*10*1e-3); colormap(usa2);

figure(4); clf
wahPL = squeeze(nanmean(reshape(planck_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahPL'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (PLANCK)']); colormap jet

figure(5); clf
wahLA = squeeze(nanmean(reshape(lapse_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahPL'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (LAPSE)']); colormap jet

figure(6); clf
wahTCO2 = squeeze(nanmean(reshape(t_co2_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahTCO2'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (TEMP+CO2)']); colormap jet
caxis([-1 +1]*2*1e-2); colormap(usa2);

figure(7); clf
wahO3 = squeeze(nanmean(reshape(o3_bands.clr-olr0_bands.clr,16,72,64),2));
pcolor(bandcenter,lat64,wahO3'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (O3)']); colormap jet
caxis([-1 +1]*2*1e-2); colormap(usa2);

wahALL_1 = squeeze(nanmean(reshape(t_co2_bands.clr-olr0_bands.clr,16,72,64),2)) + squeeze(nanmean(reshape(wv_bands.clr-olr0_bands.clr,16,72,64),2)) + squeeze(nanmean(reshape(skt_bands.clr-olr0_bands.clr,16,72,64),2)) + ...
           squeeze(nanmean(reshape(o3_bands.clr-olr0_bands.clr,16,72,64),2));
figure(08); clf; pcolor(bandcenter,lat64,wahALL_1'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (SUM ALL v1)']); colormap jet
caxis([-1 +1]*2*1e-2); colormap(usa2);

wahALL_2 = squeeze(nanmean(reshape(allpert_bands.clr-olr0_bands.clr,16,72,64),2));
figure(09); clf; pcolor(bandcenter,lat64,wahALL_2'); shading interp; colorbar; xlabel('RRTM band center [cm-1]'); ylabel('Latitude'); title([dataset_name '\delta OLR due to \delta (ALL v2)']); colormap jet
caxis([-1 +1]*2*1e-2); colormap(usa2);

figure(10); clf
plot(lat64,nansum(wahALL_1,1),lat64,nansum(wahALL_2,1),'linewidth',2); plotaxis2; ylabel('Flux Trend'); xlabel('Latitude'); hl = legend('v1 = missing g4,6','v2 = all'); title(dataset_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCosAvg = +1;
if iCosAvg > 0
  wunk = (ones(16,1) * cos(p.rlat*pi/180));
else
  wunk = ones(size(wunk));
end

%junkLA  = (lapse_bands.clr-olr0_bands.clr) .* wunk;  
%junkPL  = (planck_bands.clr-olr0_bands.clr) .* wunk;  
junkO3   = (o3_bands.clr-olr0_bands.clr) .* wunk;  
junkWV   = (wv_bands.clr-olr0_bands.clr) .* wunk;  
junkSKT  = (skt_bands.clr-olr0_bands.clr) .* wunk;  
junkTCO2 = (t_co2_bands.clr-olr0_bands.clr) .* wunk;  
junkALL  = (allpert_bands.clr-olr0_bands.clr) .* wunk;  

%meanjunkLA  = nansum(junkLA,2)./nansum(wunk,2);
%meanjunkPL  = nansum(junkPL,2)./nansum(wunk,2);
meanjunkWV   = nansum(junkWV,2)./nansum(wunk,2);
meanjunkO3   = nansum(junkO3,2)./nansum(wunk,2);
meanjunkTCO2 = nansum(junkTCO2,2)./nansum(wunk,2);
meanjunkSKT = nansum(junkSKT,2)./nansum(wunk,2);
meanjunkALL = nansum(junkALL,2)./nansum(wunk,2);
figure(1); clf; plot(bandcenter,meanjunkSKT,'g',bandcenter,meanjunkWV,'b',bandcenter,meanjunkTCO2,'r',bandcenter,meanjunkO3,'c','linewidth',2);
           hold on; plot(bandcenter,meanjunkSKT+meanjunkWV+meanjunkTCO2+meanjunkO3,'k--',bandcenter,meanjunkALL,'kx-','linewidth',4); hold off
plotaxis2; hl = legend('SKT','WV','T+CO2','O3','ALLv1 = SUM(4)','ALLv2','location','best','fontsize',10);

%figure(1); clf; plot(bandcenter,meanjunkSKT,bandcenter,meanjunkWV,bandcenter,meanjunkT,bandcenter,meanjunkSKT+meanjunkWV+meanjunkT+meanjunkO3,'kx-','linewidth',2);
%plotaxis2; hl = legend('SKT','WV','T+CO2','ALL','location','best','fontsize',10);

punk = cos(lat64*pi/180);
fprintf(1,'%s : add (cosine weighted)  flux trend = sum(T+WV+SKT+O3+GHG) = %8.4f W/m2/yr           sum(all even CH4,N2O) = %8.4f W/m2/yr \n',dataset_name,sum(meanjunkSKT+meanjunkWV+meanjunkTCO2+meanjunkO3),sum(meanjunkALL))
