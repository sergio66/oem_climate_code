addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load BASE profile

iClearMem = input('clear all the memory and start from scratch (-1 NO default/ +1 YES) : ');
if length(iClearMem) == 0
  iClearMem = 0;
end
if iClearMem > 0
  clear all
  iClearMem = 1;
  load('llsmap5');
end

if ~exist('results')
  disp('WARNING : using saved results,resultsWV,resultsT')
  savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_globalSSTfeedback.mat';
  savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';

  fprintf(1,'savename = %s \n',savename);
  load(savename);  
end

if ~exist ('h') & ~exist('p')
  %% see eg SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_gather_RH_rates_AIRSL3_NWP_XMIP.m
  [hMean17years,ha,pMean17years,pa]     = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
  iLoad = 1;
    iDorA = 1;
    if iDorA > 0
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
    else
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
    end
    era_prof = load(fin);
    hTimeStep1 = era_prof.hnew_op;
    pTimeStep1 = era_prof.pnew_op;
  %% now see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/find_T_RH_trends.m
  h = hTimeStep1; p = pTimeStep1;      %% I been using this in eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_gather_gridded_retrieval_results
  h = hMean17years; p = pMean17years;  %% I think I should use this
end

%% now that we have flux calcs, can we do some plots of how flux has changed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute OLR

clear cmip6_OLR_spectral_olr

fracWV = cmip6.trend_gas_1; fracWV(101,:) = 0;
fracO3 = cmip6.trend_gas_3; fracO3(101,:) = 0;
deltaT = cmip6.trend_ptemp; deltaT(101,:) = 0;
results(:,6) = cmip6.trend_stemp;

fracWVunc = cmip6.trend_gas_1_err; fracWVunc(101,:) = 0;
fracO3unc = cmip6.trend_gas_3_err; fracO3unc(101,:) = 0;
deltaTunc = cmip6.trend_ptemp_err; deltaTunc(101,:) = 0;
resultsunc(:,6) = cmip6.trend_stemp_err;

if iClearMem > 0
  px = p;
  cmip6_OLR_spectral_olr.olr0       = compute_olr(h,px);
  cmip6_OLR_spectral_olr.olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);    
        
  px = p;
  fracJUNK = fracWV; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
  px.gas_1 = px.gas_1 .* (1 + fracJUNK);
  cmip6_OLR_spectral_olr.wv = compute_olr(h,px);
  cmip6_OLR_spectral_olr.wv_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.wv,2)); title('water vapor change')
  figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.wv_ecRad.clr); title('water vapor change')

elseif iClearMem < 0
  cmip6_OLR_spectral_olr.olr0       = cmip6_spectral_olr.olr0;
  cmip6_OLR_spectral_olr.olr0_ecRad = cmip6_spectral_olr.olr0_ecRad;

  cmip6_OLR_spectral_olr.wv         = cmip6_spectral_olr.wv;
  cmip6_OLR_spectral_olr.olr0_ecRad = cmip6_spectral_olr.wv_ecRad;
end

%% these all need to be freshly brewed
px = p;
px.gas_2 = px.gas_2*(1+2.2/400);
px.gas_6 = px.gas_6*(1+5/1860);
cmip6_OLR_spectral_olr.tracegas       = compute_olr(h,px);
cmip6_OLR_spectral_olr.tracegas_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1); 
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.tracegas,2)); title('tracegas change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.tracegas_ecRad.clr); title('tracegas change')

px = p;
fracJUNK = fracO3; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_3 = px.gas_3 .* (1 + fracJUNK);
cmip6_OLR_spectral_olr.o3 = compute_olr(h,px);
cmip6_OLR_spectral_olr.o3_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.o3,2)); title('Ozone change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.o3_ecRad.clr); title('Ozone change')

%% remember for feedbacks this was part of lapse rate so SST was also perturbed in addition to T(z)
%% so redo it cleanly
deltaTx = deltaT; deltaTx(deltaTx > 0.1) = 0.1; deltaTx(deltaTx < -0.1) = -0.1; deltaTx(isnan(deltaTx)) = 0;
plays101 = plays; plays101(101) = plays(100)+25;
figure(3); clf; pcolor(rlat,plays101,squeeze(nanmean(reshape(deltaTx,101,72,64),2))); colorbar; set(gca,'ydir','reverse'); colormap(llsmap5); caxis([-1 +1]*0.125); shading interp; set(gca,'yscale','log'); ylim([1 1000])
px = p;
px.ptemp = px.ptemp + deltaTx;
cmip6_OLR_spectral_olr.ptemp = compute_olr(h,px);   
cmip6_OLR_spectral_olr.ptemp_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.ptemp,2)); title('T(z) change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ptemp_ecRad.clr); title('T(z) change')
dfluxT = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ptemp_ecRad.clr; dfluxT(abs(dfluxT) > 2) = NaN;
figure(2); plot(px.rlat,dfluxT)

%% remember for feedbacks SST used global SST
%% so redo it cleanly
px = p;
indSST    = results(:,6)';
px.stemp = px.stemp + indSST;
cmip6_OLR_spectral_olr.skt = compute_olr(h,px);   
cmip6_OLR_spectral_olr.skt_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.skt,2)); title('SST change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.skt_ecRad.clr); title('SST change')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = p;
px.gas_2 = px.gas_2*(1+2.2/400);
px.gas_6 = px.gas_6*(1+5/1860);

indSST    = results(:,6)';
px.stemp = px.stemp + indSST;

deltaTx = deltaT; deltaTx(deltaTx > 0.1) = 0.1; deltaTx(deltaTx < -0.1) = -0.1; deltaTx(isnan(deltaTx)) = 0;
px.ptemp = px.ptemp + deltaTx;

fracJUNK = fracWV; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_1 = px.gas_1 .* (1 + fracJUNK);

fracJUNK = fracO3; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_3 = px.gas_3 .* (1 + fracJUNK);

cmip6_OLR_spectral_olr.ALL = compute_olr(h,px);
cmip6_OLR_spectral_olr.ALL_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.ALL,2)); title('ALL change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ALL_ecRad.clr); title('ALL change')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = p;
px.gas_2 = px.gas_2*(1+2.3/400);
px.gas_6 = px.gas_6*(1+5.1/1860);

indSST    = results(:,6)' + resultsunc(:,6)';
px.stemp = px.stemp + indSST;

deltaTx = deltaT + sign(deltaT).*deltaTunc; deltaTx(deltaTx > 0.1) = 0.1; deltaTx(deltaTx < -0.1) = -0.1; deltaTx(isnan(deltaTx)) = 0;
px.ptemp = px.ptemp + deltaTx;

fracJUNK = fracWV + sign(fracWV).*fracWVunc; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_1 = px.gas_1 .* (1 + fracJUNK);

fracJUNK = fracO3 + sign(fracO3).*fracO3unc; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_3 = px.gas_3 .* (1 + fracJUNK);

cmip6_OLR_spectral_olr.ALLunc = compute_olr(h,px);
cmip6_OLR_spectral_olr.ALLunc_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
figure(1); plot(h.vchan,nanmean(cmip6_OLR_spectral_olr.olr0 - cmip6_OLR_spectral_olr.ALLunc,2)); title('ALLunc change')
figure(2); plot(px.rlat,cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ALLunc_ecRad.clr); title('ALLunc change')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

olr_delta_tracegas      = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.tracegas_ecRad.clr;
olr_delta_wv            = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.wv_ecRad.clr;
olr_delta_o3            = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.o3_ecRad.clr;
olr_delta_skt           = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.skt_ecRad.clr;
olr_delta_ptemp         = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ptemp_ecRad.clr;
olr_delta_ALL           = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ALL_ecRad.clr;
olr_delta_ALLunc        = cmip6_OLR_spectral_olr.olr0_ecRad.clr - cmip6_OLR_spectral_olr.ALLunc_ecRad.clr;

plot(rlat,nanmean(reshape(olr_delta_tracegas,72,64),1))

iNumYears = 19;

plot(rlat,nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,'yx-',rlat,nanmean(reshape(olr_delta_ALLunc,72,64),1)*iNumYears,'gx-')
olr_unc = abs((olr_delta_ALL-olr_delta_ALLunc));
plot(rlat,nanmean(reshape(olr_delta_ALL,72,64),1),'b',rlat,nanmean(reshape(olr_unc,72,64),1)/sqrt(72),'c'); plotaxis2;

plot(rlat,nanmean(reshape(olr_delta_skt,72,64),1)*iNumYears,'g',rlat,nanmean(reshape(olr_delta_ptemp,72,64),1)*iNumYears,'r',rlat,nanmean(reshape(olr_delta_wv,72,64),1)*iNumYears,'c',...
     rlat,nanmean(reshape(olr_delta_tracegas,72,64),1)*iNumYears,'m',rlat,nanmean(reshape(olr_delta_o3,72,64),1)*iNumYears,'k',rlat,nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,'yx-','linewidth',2)
  plotaxis2; hl = legend('SKT','T(z)','WV(z)','CO2/CH4','O3','ALL','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now load in CERES
ceres_fnameS = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF-TOA_Ed4.1_Subset_200209-202108.nc';  %% what I brought
ceresS = load_ceres_data(ceres_fnameS,+1);

ceres_fnameR = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF_Ed4.1_Subset_200209-202108.nc';      %% what Ryan suggests
ceresR = load_ceres_data(ceres_fnameR,-1);

ceres = ceresR;

addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

iCnt = 0;
for yyx = 2002 : 2021
  mmS = 1; mmE = 12;
  if yyx == 2002
    mmS = 09;
  elseif yyx == 2021
    mmE = 08;
  end
  for ii = mmS : mmE
    iCnt = iCnt + 1;
    all.yy(iCnt) = yyx;
    all.mm(iCnt) = ii;
    all.dd(iCnt) = 15;
  end
end
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

for ii = 1 : 180
  data = ceres.lwdata(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_lw(ii) = B(2);  
    trend_ceres_lw_err(ii) = stats.se(2);
  else
    trend_ceres_lw(ii) = NaN;
    trend_ceres_lw_err(ii) = NaN;
  end

  data = ceres.lwdata_clr(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_lw_clr(ii) = B(2);  
    trend_ceres_lw_clr_err(ii) = stats.se(2);
  else
    trend_ceres_lw_clr(ii) = NaN;
    trend_ceres_lw_clr_err(ii) = NaN;
  end
end

for ii = 1 : 180
  data = ceres.swdata(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_sw(ii) = B(2);  
    trend_ceres_sw_err(ii) = stats.se(2);
  else
    trend_ceres_sw(ii) = NaN;
    trend_ceres_sw_err(ii) = NaN;
  end

  data = ceres.swdata_clr(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_sw_clr(ii) = B(2);  
    trend_ceres_sw_clr_err(ii) = stats.se(2);
  else
    trend_ceres_sw_clr(ii) = NaN;
    trend_ceres_sw_clr_err(ii) = NaN;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS
iA    = 1;
iNorD = 1;
iAorOrL = 0;
airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL);

figure(1); clf
iNumYears = 19;
plot(rlat,nanmean(reshape(olr_delta_skt,72,64),1)*iNumYears,'g',rlat,nanmean(reshape(olr_delta_ptemp,72,64),1)*iNumYears,'r',rlat,nanmean(reshape(olr_delta_wv,72,64),1)*iNumYears,'c',...
     rlat,nanmean(reshape(olr_delta_tracegas,72,64),1)*iNumYears,'m',rlat,nanmean(reshape(olr_delta_o3,72,64),1)*iNumYears,'k',rlat,nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,'yx-','linewidth',2)
  plotaxis2; hl = legend('SKT','T(z)','WV(z)','CO2/CH4','O3','ALL','location','best','fontsize',8);
  title('CMIP6 FORCING = ORIG-NEW'); ylabel('W/m2'); xlabel('Latitude');

figure(2); clf
iNumYears = 19;
plot(rlat,-nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,'yo-',ceres.lat,trend_ceres_lw*iNumYears,'c--',ceres.lat,trend_ceres_lw_clr*iNumYears,'b',...
     rlat,nanmean(airsChoice.thestats64x72_other.olrrate,1)*iNumYears,'m--',rlat,nanmean(airsChoice.thestats64x72_other.clrolrrate,1)*iNumYears,'r','linewidth',2)
  plotaxis2; hl = legend('CMIP6','CERES cld','CERES clr','AIRS L3 cld','AIRS L3 clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');

figure(3); clf; plot(ceres.lat,trend_ceres_lw*iNumYears,'c--',ceres.lat,trend_ceres_lw_clr*iNumYears,'b',ceres.lat,trend_ceres_sw*iNumYears,'m--',ceres.lat,trend_ceres_sw_clr*iNumYears,'r','linewidth',2)
  plotaxis2; hl = legend('CERES LW cld','CERES LW clr','CERES SW cld','CERES SW clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');
figure(4); clf; plot(ceres.lat,(trend_ceres_lw+trend_ceres_sw)*iNumYears,'b--',ceres.lat,(trend_ceres_lw_clr+trend_ceres_sw_clr)*iNumYears,'r','linewidth',2)
  plotaxis2; hl = legend('CERES LW+SW cld','CERES LW+SW clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');

olr_unc64     = nanmean(reshape(olr_unc,72,64),1)/sqrt(72);
airs_unc64    = nanmean(airsChoice.thestats64x72_other.olrratestd,1)/sqrt(72);
airsclr_unc64 = nanmean(airsChoice.thestats64x72_other.clrolrratestd,1)/sqrt(72);
figure(5); clf;
  hold on; errorbar(rlat,-nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,olr_unc64*iNumYears,'color','y','linewidth',2);
  hold on; errorbar(ceres.lat,trend_ceres_lw*iNumYears,trend_ceres_lw_err*iNumYears,'color','c','linewidth',2);
  hold on; errorbar(ceres.lat,trend_ceres_lw_clr*iNumYears,trend_ceres_lw_clr_err*iNumYears,'color','b','linewidth',2);
  hold on; errorbar(rlat,nanmean(airsChoice.thestats64x72_other.olrrate,1)*iNumYears,airs_unc64*iNumYears,'color','m','linewidth',2);
  hold on; errorbar(rlat,nanmean(airsChoice.thestats64x72_other.clrolrrate,1)*iNumYears,airsclr_unc64*iNumYears,'color','r','linewidth',2);
  hold off  
  plotaxis2; hl = legend('CMIP6','CERES cld','CERES clr','AIRS L3 cld','AIRS L3 clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');

figure(6); clf;
  hold on; errorbar(rlat,-nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,olr_unc64*iNumYears,'color','r','linewidth',2,'marker','o');
  hold on; errorbar(ceres.lat,trend_ceres_lw*iNumYears,trend_ceres_lw_err*iNumYears,'color','c','linewidth',2);
  hold on; errorbar(ceres.lat,trend_ceres_lw_clr*iNumYears,trend_ceres_lw_clr_err*iNumYears,'color','b','linewidth',2);
  hold off  
  plotaxis2; hl = legend('CMIP6','CERES cld','CERES clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');
%% figure(6); aslprint('fluxOLR_witherrorbars.pdf')

figure(7); clf;
  hold on; plot(rlat,-nanmean(reshape(olr_delta_ALL,72,64),1)*iNumYears,'r','linewidth',2,'marker','o');
  hold on; plot(ceres.lat,trend_ceres_lw*iNumYears,'c','linewidth',2);
  hold on; plot(ceres.lat,trend_ceres_lw_clr*iNumYears,'b','linewidth',2);
  hold off  
  plotaxis2; hl = legend('CMIP6','CERES cld','CERES clr','location','best','fontsize',8);
  title('\delta Flux = Final-Orig'); ylabel('W/m2'); xlabel('Latitude');
%% figure(7); aslprint('fluxOLR.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the data

ceres_trends.lat                    = ceres.lat;
ceres_trends.trend_ceres_lw         = trend_ceres_lw;
ceres_trends.trend_ceres_lw_clr     = trend_ceres_lw_clr;
ceres_trends.trend_ceres_sw         = trend_ceres_sw;
ceres_trends.trend_ceres_sw_clr     = trend_ceres_sw_clr;
ceres_trends.trend_ceres_lw_err     = trend_ceres_lw_err;
ceres_trends.trend_ceres_lw_clr_err = trend_ceres_lw_clr_err;
ceres_trends.trend_ceres_sw_err     = trend_ceres_sw_err;
ceres_trends.trend_ceres_sw_clr_err = trend_ceres_sw_clr_err;

airsL3_trends.rlat = rlat;
airsL3_trends.airsL3         = nanmean(airsChoice.thestats64x72_other.olrrate,1);
airsL3_trends.airsL3_clr     = nanmean(airsChoice.thestats64x72_other.clrolrrate,1);
airsL3_trends.airsL3_unc     = airs_unc64;
airsL3_trends.airsL3_clr_unc = airsclr_unc64;

cmip6_trends.rlat         = rlat;
cmip6_trends.cmip6_all     = -nanmean(reshape(olr_delta_ALL,72,64),1);
cmip6_trends.cmip6_all_err = olr_unc64;
cmip6_trends.cmip6_skt     = nanmean(reshape(olr_delta_skt,72,64),1);
cmip6_trends.cmip6_ptemp   = nanmean(reshape(olr_delta_ptemp,72,64),1);
cmip6_trends.cmip6_wv      = nanmean(reshape(olr_delta_wv,72,64),1);
cmip6_trends.cmip6_o3      = nanmean(reshape(olr_delta_o3,72,64),1);
cmip6_trends.cmip6_co2_ch4 = nanmean(reshape(olr_delta_tracegas,72,64),1);

outflux.ceres_trends  = ceres_trends;
outflux.airsL3_trends = airsL3_trends;
outflux.cmip6_trends   = cmip6_trends;

fprintf(1,'CERES LW area weighted 19 year flux   cld = %8.6f clr = %8.6f W/m2 \n',[sum(cos(ceres.lat*pi/180).*trend_ceres_lw'*iNumYears)   sum(cos(ceres.lat*pi/180).*trend_ceres_lw_clr'*iNumYears)]/sum(cos(ceres.lat*pi/180)))
fprintf(1,'CERES SW area weighted 19 year flux   cld = %8.6f clr = %8.6f W/m2 \n',[sum(cos(ceres.lat*pi/180).*trend_ceres_sw'*iNumYears)   sum(cos(ceres.lat*pi/180).*trend_ceres_sw_clr'*iNumYears)]/sum(cos(ceres.lat*pi/180)))
fprintf(1,'AIRS L3  area weighted 19 year flux   cld = %8.6f clr = %8.6f W/m2 \n',[sum(cos(rlat*pi/180).*airsL3_trends.airsL3'*iNumYears)  sum(cos(rlat*pi/180).*airsL3_trends.airsL3_clr'*iNumYears)]/sum(cos(rlat*pi/180)))
fprintf(1,'CMIP6    area weighted 19 year flux   cld = %8.6f clr = %8.6f W/m2 \n',[NaN  sum(cos(rlat*pi/180).*cmip6_trends.cmip6_all'*iNumYears)]/sum(cos(rlat*pi/180)))

%{
saveOLRname = [savename(1:end-4) '_olr_cmip6_vs_ceres.mat'];
saver = ['save ' saveOLRname ' outflux'];
eval(saver);
%}