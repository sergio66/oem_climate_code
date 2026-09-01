addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');
imagesc(reshape(p.stemp,72,64)'); colorbar; colormap jet
for ii = 2 : 64  
  %% start with 2 since no ocean at rlat==1
  ind = (1:72) + (ii-1)*72;
  ocean = find(p.landfrac(ind) == 0);
  stemp(ii-1) = nanmean(p.stemp(ind));
  stempOcean(ii-1) = nanmean(p.stemp(ind(ocean)));
  rlat(ii-1) = nanmean(p.rlat(ind));
  pOcean.rlat(ii-1) = nanmean(p.rlat(ind(ocean)));
  pOcean.rlon(ii-1) = nanmean(p.rlon(ind(ocean)));
  pOcean.plat(ii-1) = nanmean(p.plat(ind(ocean)));
  pOcean.plon(ii-1) = nanmean(p.plon(ind(ocean)));
  pOcean.scanang(ii-1) = nanmean(p.scanang(ind(ocean)));
  pOcean.satzen(ii-1) = nanmean(p.satzen(ind(ocean)));
  pOcean.solzen(ii-1) = nanmean(p.solzen(ind(ocean)));
  pOcean.nemis(ii-1) = nanmean(p.nemis(ind(ocean)));
  pOcean.efreq(:,ii-1) = nanmean(p.efreq(:,ind(ocean)),2);
  pOcean.emis(:,ii-1) = nanmean(p.emis(:,ind(ocean)),2);
  pOcean.rho(:,ii-1) = nanmean(p.rho(:,ind(ocean)),2);
  pOcean.stemp(ii-1) = nanmean(p.stemp(ind(ocean)));
  pOcean.nlevs(ii-1) = floor(nanmean(p.nlevs(ind(ocean))));
  pOcean.rtime(ii-1) = nanmean(p.rtime(ind(ocean)));
  pOcean.landfrac(ii-1) = nanmean(p.landfrac(ind(ocean)));
  pOcean.spres(ii-1) = nanmean(p.spres(ind(ocean)));
  pOcean.salti(ii-1) = nanmean(p.salti(ind(ocean)));
  pOcean.zobs(ii-1) = nanmean(p.zobs(ind(ocean)));
  pOcean.plevs(:,ii-1) = nanmean(p.plevs(:,ind(ocean)),2);
  pOcean.palts(:,ii-1) = nanmean(p.palts(:,ind(ocean)),2);
  pOcean.ptemp(:,ii-1) = nanmean(p.ptemp(:,ind(ocean)),2);
  pOcean.gas_1(:,ii-1) = nanmean(p.gas_1(:,ind(ocean)),2);
  pOcean.gas_2(:,ii-1) = nanmean(p.gas_2(:,ind(ocean)),2);
  pOcean.gas_3(:,ii-1) = nanmean(p.gas_3(:,ind(ocean)),2);
  pOcean.gas_4(:,ii-1) = nanmean(p.gas_4(:,ind(ocean)),2);
  pOcean.gas_5(:,ii-1) = nanmean(p.gas_5(:,ind(ocean)),2);
  pOcean.gas_6(:,ii-1) = nanmean(p.gas_6(:,ind(ocean)),2);
  pOcean.gas_9(:,ii-1) = nanmean(p.gas_9(:,ind(ocean)),2);
  pOcean.gas_12(:,ii-1) = nanmean(p.gas_12(:,ind(ocean)),2);
  pOcean.rcalc(:,ii-1) = nanmean(p.rcalc(:,ind(ocean)),2);
end
plot(rlat,stemp,rlat,stempOcean)

pOcean.spres = 1013 * ones(size(pOcean.spres));
pOcean.nlevs = 98   * ones(size(pOcean.spres));

cdRRTMback = ['cd ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iLambda_UseGlobalSST')
  iLambda_UseGlobalSST = +1;   %% new way, use global avg SST, much safer (less likely to be 0)
  iLambda_UseGlobalSST = -1;   %% old way till March 2022, using computed SST per tile instead of glabal average; dangerous because if dSST = 0 oopie when dividing
end

addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/

px = pOcean;
meanOLR.olr0 = compute_olr(h,px);
meanOLR.olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);
plot(pOcean.rlat,meanOLR.olr0_ecRad.clr); title('ECRAD OLR0'); ylabel('W/m2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/PLOTTER
rlat41 = equal_area_spherical_bands(20);
rlat40 = meanvaluebin(rlat41);

for ii = 1 : 40
  pOcean40.rlat(ii) = rlat40(ii);
  pOcean40.rlon(ii) = interp1(pOcean.rlat,pOcean.rlon,pOcean40.rlat(ii),[],'extrap');
  pOcean40.plat(ii) = interp1(pOcean.rlat,pOcean.plat,pOcean40.rlat(ii),[],'extrap');
  pOcean40.plon(ii) = interp1(pOcean.rlat,pOcean.plon,pOcean40.rlat(ii),[],'extrap');
  pOcean40.scanang(ii) = interp1(pOcean.rlat,pOcean.scanang,pOcean40.rlat(ii),[],'extrap');
  pOcean40.satzen(ii) = interp1(pOcean.rlat,pOcean.satzen,pOcean40.rlat(ii),[],'extrap');
  pOcean40.solzen(ii) = interp1(pOcean.rlat,pOcean.solzen,pOcean40.rlat(ii),[],'extrap');

  pOcean40.stemp(ii) = interp1(pOcean.rlat,pOcean.stemp,pOcean40.rlat(ii),[],'extrap');
  pOcean40.nlevs(ii) = interp1(pOcean.rlat,pOcean.nlevs,pOcean40.rlat(ii),[],'extrap');
  pOcean40.rtime(ii) = interp1(pOcean.rlat,pOcean.rtime,pOcean40.rlat(ii),[],'extrap');
  pOcean40.landfrac(ii) = interp1(pOcean.rlat,pOcean.landfrac,pOcean40.rlat(ii),[],'extrap');
  pOcean40.spres(ii) = interp1(pOcean.rlat,pOcean.spres,pOcean40.rlat(ii),[],'extrap');
  pOcean40.salti(ii) = interp1(pOcean.rlat,pOcean.salti,pOcean40.rlat(ii),[],'extrap');
  pOcean40.zobs(ii) = interp1(pOcean.rlat,pOcean.zobs,pOcean40.rlat(ii),[],'extrap');
end

for ii = 1 : 40
  pOcean40.nemis(ii) = interp1(pOcean.rlat,pOcean.nemis,pOcean40.rlat(ii),[],'extrap');
  for jj = 1 : 19
    pOcean40.efreq(jj,ii) = interp1(pOcean.rlat,pOcean.efreq(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.emis(jj,ii) = interp1(pOcean.rlat,pOcean.emis(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.rho(jj,ii) = interp1(pOcean.rlat,pOcean.rho(jj,:),pOcean40.rlat(ii),[],'extrap');
  end

  for jj = 1 : 101
    pOcean40.palts(jj,ii) = interp1(pOcean.rlat,pOcean.palts(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.plevs(jj,ii) = interp1(pOcean.rlat,pOcean.plevs(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.ptemp(jj,ii) = interp1(pOcean.rlat,pOcean.ptemp(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_1(jj,ii) = interp1(pOcean.rlat,pOcean.gas_1(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_2(jj,ii) = interp1(pOcean.rlat,pOcean.gas_2(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_3(jj,ii) = interp1(pOcean.rlat,pOcean.gas_3(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_4(jj,ii) = interp1(pOcean.rlat,pOcean.gas_4(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_5(jj,ii) = interp1(pOcean.rlat,pOcean.gas_5(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_6(jj,ii) = interp1(pOcean.rlat,pOcean.gas_6(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_9(jj,ii) = interp1(pOcean.rlat,pOcean.gas_9(jj,:),pOcean40.rlat(ii),[],'extrap');
    pOcean40.gas_12(jj,ii) = interp1(pOcean.rlat,pOcean.gas_12(jj,:),pOcean40.rlat(ii),[],'extrap');
  end
end

for ii = 1 : 40
  for jj = 1 : 2645
    pOcean40.rcalc(jj,ii) = interp1(pOcean.rlat,pOcean.rcalc(jj,:),pOcean40.rlat(ii),[],'extrap');
  end
end
ohoh = find((pOcean40.nlevs-floor(pOcean40.nlevs)) ~= 0)

pOcean40.nlevs = round(pOcean40.nlevs);
px = pOcean40;
meanOLR40.olr0 = compute_olr(h,px);
meanOLR40.olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);
plot(pOcean.rlat,meanOLR.olr0_ecRad.clr,pOcean40.rlat,meanOLR40.olr0_ecRad.clr); title('ECRAD OLR0'); ylabel('W/m2')

pOcean40.plays = plevs2plays(pOcean40.plevs);
pcolor(pOcean40.rlat,log10(nanmean(pOcean40.plays(1:97,:),2)),pOcean40.ptemp(1:97,:)); colormap jet; colorbar
pcolor(pOcean40.rlat,nanmean(pOcean40.plays(1:97,:),2),pOcean40.ptemp(1:97,:)); colormap jet; colorbar; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); shading interp
pcolor(pOcean40.rlat,nanmean(pOcean40.plays(1:97,:),2),log10(pOcean40.gas_1(1:97,:))); colormap jet; colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see eg driver_olr_fluxchanges_UMBC.m
anom40_ST = load('../AIRS_new_clear_scan_August2019_AMT2020PAPER/anomaly_0dayavg_cal_results_LW_fewerN2ochans.mat');
anom40_T  = load('../AIRS_new_clear_scan_August2019_AMT2020PAPER/era_ptempanom.mat');  %% REALLY NEED TO GET A BETTER ONE

addpath ../../oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
warning off
for ii = 1 : 40
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  junk = squeeze(anom40_T.era_ptempanom(ii,:,:));
  datime = (1:16:365*16);
  for jj = 1 : 20
    data = junk(jj,:);
    [B, stats, err] = Math_tsfit_lin_robust(datime,data,4);
    trendT(jj,ii) = B(2);
  end
end
fprintf(1,'\n');
warning on
pcolor(trendT); shading interp; colorbar
p101 = load('/home/sergio/MATLABCODE/airslevels.dat');
plays100 = plevs2plays(flipud(p101));
for jj = 1 : 20
  ind = (1:5) + (jj-1)*5;
  plays20(jj) = mean(plays100(ind));
end
addpath /home/sergio/MATLABCODE/COLORMAP/
pcolor(pOcean40.rlat,plays20,trendT); shading interp; colormap(usa2); caxis([-1 +1]*0.15); 
pcolor(pOcean40.rlat,plays20,trendT); shading interp; colormap(usa2); caxis([-1 +1]*0.05); xlim([-1 +1]*60);
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp; colorbar
load /home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat;  colormap(llsmap5);

airsL3_spectral_olr.planck = compute_olr(h,px);
airsL3_spectral_olr.planck_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% for ii = 1 : 40
%%   junk = squeeze(anom40_T.era_ptempanom(ii,:,:));
%%   for jj = 1 : 20
%%     ind = (1:5) + (jj-1)*5;
%%     deltaT(ind,:) = ones(5,1)*junk(jj,:);
%%   end
%%   deltaTx = deltaT; 
%%   %deltaTx(deltaTx > 0.1) = 0.1; deltaTx(deltaTx < -0.1) = -0.1; deltaTx(isnan(deltaTx)) = 0;
%%   [hxx,pxx] = replicate_rtp_headprof(h,pOcean40,ii,365);
%%   pxx.ptemp(1:100,:) = pxx.ptemp(1:100,:) + deltaTx;
%%   junk = superdriver_run_ecRad_rtp_loop_over_profiles(hxx,pxx,-1);               
%%   olr40(ii,:) = junk.clr;
%% end
%% 
%% pcolor(olr40); shading interp; colorbar
%% for ii = 1 : 40
%%   [B, stats, err] = Math_tsfit_lin_robust(datime,olr40(ii,:),4);
%%   olrtrend(ii) = B(2);
%% end
%% plot(pOcean40.rlat,smoothn(-olrtrend,5)); plotaxis2;
%% 
%% 
%% olr0   = meanOLR40.olr0_ecRad.clr' * ones(1,365);
%% dstemp = anom40_ST.stemp;
%% 
%% dxstemp = dstemp; dxstemp(abs(dxstemp) < 0.001) = 0.001;
%% %dxstemp = nanmean(dstemp(:)) * ones(size(dstemp));
%% 
%% %% see eg compute_feedbacks_umbc_rrtm.m
%% planck_feedback = (olr40 - olr0)./dxstemp;
%% pcolor(planck_feedback); shading interp; colorbar; caxis([-3 +3])
%% plot(pOcean40.rlat,-nanmean(planck_feedback,2)); ylim([-1 +1]*10); plotaxis2;
%% plot(pOcean40.rlat,-smoothn(nanmean(planck_feedback,2),10)); ylim([-1 +1]*10); plotaxis2;
%% 
%% figure(7); plot(dstemp,olr40-olr0,'.')
%% addpath /home/sergio/MATLABCODE/SHOWSTATS/
%% [n,nx,ny,nmean,nstd] = myhist2d(dxstemp(:),olr40(:)-olr0(:),-1:0.1:+1,-1:0.1:+1);
%% errorbar(-1:0.1:+1,nmean,nstd); xlabel('\Delta ST (K)'); ylabel('\Delta OLR W/m2')
%% 
%% addpath /home/sergio/MATLABCODE/NANROUTINES
%% for ii = 1 : 40
%%   a = dxstemp(ii,:);
%%   b = olr40(ii,:)-olr0(ii,:);
%%   oo = find(isfinite(a) & isfinite(b));
%%   junk = polyfit(a(oo),b(oo),1);
%%   slope(ii) = junk(1);
%% end
%% plot(pOcean40.rlat,slope,pOcean40.rlat,smoothn(nanmean(planck_feedback,2),10))
%% 
%% oo = find(isfinite(dstemp(:)) & isfinite(olr40(:)-olr0(:)));
%% a = dstemp(:);
%% b = olr40(:)-olr0(:);
%% figure(7); plot(a(oo),b(oo),'.'); nanpolyfit(a(oo),b(oo),1)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see eg compute_feedbacks_airsL3_ecRad.m
%%  px.stemp = px.stemp + globalSST;
%% px.ptemp = px.ptemp + ones(101,4608)*globalSST;
px = pOcean40;
px.ptemp = px.ptemp + 1;
px.stemp = px.stemp + 1;
junk = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
olr0   = meanOLR40.olr0_ecRad.clr;
olr40 = junk.clr;
%% see eg compute_feedbacks_umbc_rrtm.m
planck_feedback = (olr40 - olr0)./1;
plot(pOcean40.rlat,-smoothn(planck_feedback,10)); ylim([-1 +1]*10); plotaxis2;


globalSST = nanmean(anom40_ST.stemp(:));
px = pOcean40;
px.ptemp = px.ptemp + globalSST;
px.stemp = px.stemp + globalSST;
junk = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
olr0   = meanOLR40.olr0_ecRad.clr;
olr40 = junk.clr;
%% see eg compute_feedbacks_umbc_rrtm.m
planck_feedback = (olr40 - olr0)./globalSST;
plot(pOcean40.rlat,-smoothn(planck_feedback,10)); ylim([-1 +1]*10); plotaxis2;
