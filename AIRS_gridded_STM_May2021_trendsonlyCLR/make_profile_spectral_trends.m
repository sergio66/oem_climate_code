function nwp_spectral_trends = make_profile_spectral_trends(era,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,iERAorCMIP6,iVersJac,iDoSpectralRates);

%% inputs 
%%   era_or_cmip6 = 100 layers
%%   era5         = 100 layers
%%   airsL3       = JOEL : 24 (T,O3) or 12 (WV) layers  CHRIS : 100 (T,O3) or 66 (WV) layers
%%   results*     = 6 scalars, default T(z),WV(z),O3(z), typically 20 layers but could be eg 10 or 25 or 50 or 97
%%   iERAorCMIP6  = 1  : using ERA
%%                = 2  : using CMIP6
%%                = -1 : who cares, just do the calcs for first argument ie you have sent in dummies for others
%%   iVersJac     = 2019 or ERA5 2002-2021 
%%   iDoSpectralRates = optional +1 : do sum(jac x dX/dt .. slow) or 
%%                               -1 : ignore this step .. fast
%% 
%% note "mask" is not part of argument list, this routine computes spectral rates for all 4608 profiles
%% you apply the mask to the results externally after this routine has been called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 12
  iERAorCMIP6 = +1;      %% assume we are using ERA fields
  iVersJac    = 2021;    %% assume we are using ERA5 2002-2021 jacs
  iDoSpectralRates = +1; %% compute (sum(jac x dX/dt)
elseif nargin == 13
  iVersJac    = 2021;  %% assume we are using ERA5 2002-2021 jacs
  iDoSpectralRates = +1; %% compute (sum(jac x dX/dt)
elseif nargin == 14
  iDoSpectralRates = +1; %% compute (sum(jac x dX/dt)
end

umbc_spectral_rates     = zeros(2645,4608);
era_spectral_rates      = zeros(2645,4608);
era5_spectral_rates     = zeros(2645,4608);
airsL3_spectral_rates   = zeros(2645,4608);

iERA5orERAI = 2019; %% ERA-I 2002-2019
iERA5orERAI = 2021; %% ERA5  2002-2021
iERA5orERAI = iVersJac;
if iERA5orERAI == 2019
  [h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp'); %% line 23 of see_clust_put_together_jacs_clr.m
elseif iERA5orERAI == 2021
  [h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp'); %% line 23 of see_clust_put_together_jacs_clr_ERA5.m
end
[~,kcarta.subjac.ppmv2] = layers2ppmv(h1_4608,p1_4608,1:4608,2);
[~,kcarta.subjac.ppmv4] = layers2ppmv(h1_4608,p1_4608,1:4608,4);
[~,kcarta.subjac.ppmv6] = layers2ppmv(h1_4608,p1_4608,1:4608,6);

[mmUMBC,nnUMBC] = size(resultsWV);    
if nnUMBC ~= 20
  fprintf(1,' <<<<<<<<<<< WARNING make_profile_spectral_trends.m has length of WV,T,O3 retrievals as %3i and not 20 \n',nnUMBC);
end

%% see get_rates.m
for ii = 1 : 64
  ind = (ii-1)*72 + (1:72);
  % this is assuming I am reading in  
  strlatbin = num2str(ii,'%02d');
  xdriver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];
  load(xdriver.rateset.datafile)
  xdriver.rateset.rates(:,ind) = real(thesave.xtrend);
  xdriver.rateset.unc_rates(:,ind) = real(thesave.xtrendErr);
end
% plot(nanmean(xdriver.rateset.rates,2)); error('lskjglskjgs')

if iDoSpectralRates == +1
  display('  calling compute_jac_spectral_rates_N_put_airs24_cris66_to_100layers')
  compute_jac_spectral_rates_N_put_airs24_cris66_to_100layers
elseif iDoSpectralRates == -1
  display('  calling compute_put_airs24_cris66_to_100layers')
  compute_put_airs24_cris66_to_100layers
end

fprintf(1,'\n')

nwp_spectral_trends.era_spectral_rates     = era_spectral_rates;
nwp_spectral_trends.era_100_layertrends    = era_100_layertrends;

if iERAorCMIP6 > 0
  nwp_spectral_trends.umbc_spectral_rates = umbc_spectral_rates;
  nwp_spectral_trends.umbc_20_layertrends     = umbc_20_layertrends;

  nwp_spectral_trends.era5_spectral_rates    = era5_spectral_rates;
  nwp_spectral_trends.era5_100_layertrends   = era5_100_layertrends;  
  
  nwp_spectral_trends.airsL3_spectral_rates  = airsL3_spectral_rates;
  nwp_spectral_trends.airsL3_100_layertrends = airsL3_100_layertrends;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iERAorCMIP6 == 2
  %% this was CMIP6 so rename
  nwp_spectral_trends.cmip6_spectral_rates = nwp_spectral_trends.era_spectral_rates; nwp_spectral_trends = rmfield(nwp_spectral_trends,'era_spectral_rates');
  nwp_spectral_trends.cmip6_100_layertrends = nwp_spectral_trends.era_100_layertrends; nwp_spectral_trends = rmfield(nwp_spectral_trends,'era_100_layertrends');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see gather_gridded_retrieval_results_plots
%{
figure(9)
%% note we have removed MASK from argument list here since calcs above assume 4608 geophysiccal params

%% testing
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(fits(:,mask),2),'r',f,nanmean(umbc_spectral_rates(:,mask),2),'k')
plot(f,nanmean(rates(:,mask)-fits(:,mask),2),'b',f,nanstd(rates(:,mask)-fits(:,mask),[],2),'c',f,nanmean(rates(:,mask)-umbc_spectral_rates(:,mask),2),'r',f,nanstd(rates(:,mask)-umbc_spectral_rates(:,mask),[],2),'m')

plot(f,rates(:,mask),'b',f,airsL3_spectral_rates(:,mask),'c',f,era_spectral_rates(:,mask),'g',f,era5_spectral_rates(:,mask),'r',f,fits(:,mask),'m')
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(fits(:,mask),2),'m',f,nanmean(umbc_spectral_rates(:,mask),2),'k')
plot(f,nanmean(fits(:,mask)-umbc_spectral_rates(:,mask),2),'k',f,nanstd(fits(:,mask)-umbc_spectral_rates(:,mask),[],2),'g--')

%% this is DA PLOT
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),'c',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era_spectral_rates(:,mask),2),'g',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'r',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.fits(:,mask),2),'m','linewidth',2)
plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.rates(:,mask),2),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),'c',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era_spectral_rates(:,mask),2),'g',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'r',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),'m','linewidth',2)
  plotaxis2; hl = legend('AIRS Obs','AIRS L3','ERA','ERA5','quick test','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');
%}

disp('swapping ERA5 reconstructed spectral rates with SARTA spectral rates');
plot(f,nanmean(xdriver.rateset.rates,2),f,nanmean(nwp_spectral_trends.era5_spectral_rates - xdriver.rateset.rates,2),f,nanstd(nwp_spectral_trends.era5_spectral_rates - xdriver.rateset.rates,[],2)); xlim([640 1640]); grid
  hl = legend('from SARTA','SARTA-reconstruction mean','SARTA-reconstruction std','location','best','fontsize',10);
nwp_spectral_trends.era5_spectral_rates = xdriver.rateset.rates;
%plot(nanmean(nwp_spectral_trends.era5_spectral_rates - xdriver.rateset.rates,2)

