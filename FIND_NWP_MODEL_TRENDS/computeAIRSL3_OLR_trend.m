yy0 = 2002; mm0 = 8;
yy = yy0; mm = mm0;
for JOBB = 1 : 240
  yy = yy; mm = mm + 1;
  if mod(JOBB-5,12) == 0
    yy = yy + 1;
    mm = 1;
  end
  saveYY(JOBB) = yy; saveMM(JOBB) = mm;
  fprintf(1,'JOBB = %3i %4i/%2i \n',JOBB,yy,mm)
end

if ~exist('xolr') & ~exist('OLR_ecRad/AIRSL3/all_airsl3_olr.mat')
  disp('need to load 12*20 = 240 files "+" = 100, "." = 10')
  for ii = 1 : 240
    if mod(ii,100) == 0
      fprintf(1,'+')
    elseif mod(ii,10) == 0
      fprintf(1,'.')
    end
  
    fnameOUT = ['OLR_ecRad/AIRSL3/airsL3_olr_' num2str(ii,'%03d') '.mat'];
  
    if exist(fnameOUT)
      loader = ['load ' fnameOUT];
      eval(loader);
      xolr(ii,:) = olr.clr;
      xolr_airsL3(ii,:) = olr_AIRSL3;
      xrtime(ii,:) = utc2taiSergio(saveYY(ii),saveMM(ii),15,1.5);
      xstemp(ii,:) = stemp;
      xmmw(ii,:) = mmw;
  
    else
      xolr(ii,:) = NaN;
      xolr_airsL3(ii,:) = NaN;
      xrtime(ii,:) = utc2taiSergio(saveYY(ii),saveMM(ii),15,1.5);
      xstemp(ii,:) = NaN;
      xmmw(ii,:) = NaN;
    end
    
    clear fnameIN fnameOUT
  end
  save OLR_ecRad/AIRSL3/all_airsl3_olr.mat xmmw xrtime xstemp xolr*
else
  load OLR_ecRad/AIRSL3/all_airsl3_olr.mat
end

figure(1); pcolor(xstemp); colormap(jet); shading interp; colorbar; title('stemp');
figure(2); pcolor(xolr);   colormap(jet); shading interp; colorbar; title('olr')
figure(3); pcolor(xmmw);   colormap(jet); shading interp; colorbar; title('mmw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stopY = 2022;
iNumYears = 20;

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

if ~exist('trend_olr')
  disp('making trends')
  iCnt = 0;
  for yyx = 2002 : stopY
    mmS = 1; mmE = 12;
    if yyx == 2002
      mmS = 09;
    elseif yyx == stopY
      mmE = 08;
    end
    for ii = mmS : mmE
      iCnt = iCnt + 1;
      thetime.yy(iCnt) = yyx;
      thetime.mm(iCnt) = ii;
      thetime.dd(iCnt) = 15;
    end
  end
  dayOFtime = change2days(thetime.yy,thetime.mm,thetime.dd,2002);
  
  for ii = 1 : 4608
    data = xolr(:,ii);
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
      trend_olr(ii) = B(2);
      trend_olr_err(ii) = stats.se(2);
    else
      trend_olr(ii) = NaN;
      trend_olr_err(ii) = NaN;
    end
  
    data = xolr_airsL3(:,ii);
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
      trend_olr_airsL3(ii) = B(2);
      trend_olr_airsL3_err(ii) = stats.se(2);
    else
      trend_olr_airsL3(ii) = NaN;
      trend_olr_airsL3_err(ii) = NaN;
    end
  
    data = xstemp(:,ii);
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
      trend_stemp(ii) = B(2);
      trend_stemp_err(ii) = stats.se(2);
    else
      trend_stemp(ii) = NaN;
      trend_stemp_err(ii) = NaN;
    end

    data = xmmw(:,ii);
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
      trend_mmw(ii) = B(2);
      trend_mmw_err(ii) = stats.se(2);
    else
      trend_mmw(ii) = NaN;
      trend_mmw_err(ii) = NaN;
    end 
  end
  
  save OLR_ecRad/AIRSL3/all_airsl3_olr.mat xrtime xstemp xmmw xolr trend*
else
  load OLR_ecRad/AIRSL3/all_airsl3_olr.mat
end

do_XX_YY_from_X_Y

figure(2); clf
plot(rlat,nanmean(reshape(trend_stemp,72,64),1),'b',rlat,nanmean(reshape(trend_olr,72,64),1),'r',...
                                                    rlat,nanmean(reshape(trend_olr_airsL3,72,64),1),'k',...
                                                    rlat,nanmean(reshape(trend_mmw,72,64),1),'c',...
                                                    'linewidth',2)
%plot(rlat,nanmean(reshape(trend_stemp,64,72),2),'b',rlat,nanmean(reshape(trend_olr,64,72),2),'r',...
%                                                    rlat,nanmean(reshape(trend_olr_airsL3,64,72),2),'k','linewidth',2)
plotaxis2; hl = legend('stemp','OLR using ecRad+AIRS L3 profiles','OLR from AIRS L3 directly','mmw','location','best','fontsize',10);
xlim([-90 +90])
xlabel('Latitude');

figure(3); clf
ctmp = coast;

[h0,ha,p0,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.ip.rtp');
pcolor(reshape(p0.rlon,72,64),reshape(p0.rlat,72,64),reshape(trend_stemp,72,64)); colormap(usa2); caxis([-1 +1]*0.15); title('SKT trend K/yr'); shading interp;
%pcolor(reshape(p0.rlon,72,64)',reshape(p0.rlat,72,64)',reshape(trend_stemp,64,72)); colormap(usa2); caxis([-1 +1]*0.15); title('SKT trend K/yr'); shading interp;
hold on; plot(ctmp(:,2), ctmp(:,1), 'k','linewidth',2); grid on; hold off

figure(4); clf
pcolor(reshape(p0.rlon,72,64),reshape(p0.rlat,72,64),reshape(trend_mmw,72,64)); colormap(usa2); caxis([-1 +1]*0.15); title('MMW trend K/yr'); shading interp;
%pcolor(reshape(p0.rlon,72,64)',reshape(p0.rlat,72,64)',reshape(trend_mmw,64,72)); colormap(usa2); caxis([-1 +1]*0.15); title('MMW trend K/yr'); shading interp;
hold on; plot(ctmp(:,2), ctmp(:,1), 'k','linewidth',2); grid on; hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see ../AIRS_gridded_STM_May2021_trendsonlyCLR/compare_OLR_trend.m

junk = load(['../AIRS_gridded_STM_May2021_trendsonlyCLR/ceres_trends_' num2str(iNumYears,'%02d') 'year_T.mat']);
  ceres_trend = junk.ceres_trend;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;

%% axax = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat');
miaowX0 = airsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
miaowX1 = airsL3_spectral_olr.perts9999.atm_skt_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
miaowX2 = airsL3_spectral_olr.perts9999.atm_only_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
miaow = miaowX1;
miaow = miaowX0;

figure(5); clf; 
plot(ceres_trend.trend_lat,ceres_trend.trend_lw,'y:',ceres_trend.trend_lat,ceres_trend.trend_lw_clr,...
                 meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'kx-',meanvaluebin(latB2),nanmean(reshape(trend_olr,72,64),1),'gx-','linewidth',2);
plotaxis2; legend('CERES allsky','CERES clrsky','Sergio AIRSL3 profile trends -> OLR trends','Sergio ecRad AIRSL3 timeseries','location','best','fontsize',10);
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr');

figure(5); plot(ceres_trend.trend_lat,ceres_trend.trend_lw,'y:',ceres_trend.trend_lat,ceres_trend.trend_lw_clr,...
                 meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'kx-',meanvaluebin(latB2),nanmean(reshape(trend_olr,72,64),1),'gx-',...
                 meanvaluebin(latB2),nanmean(reshape(trend_olr_airsL3,72,64),1),'b-','linewidth',2);
hold on;
shadedErrorBar(ceres_trend.trend_lat,ceres_trend.trend_lw,ceres_trend.trend_lw_err,'y',0.05);
shadedErrorBar(ceres_trend.trend_lat,ceres_trend.trend_lw_clr,ceres_trend.trend_lw_clr_err,'r',0.05);
shadedErrorBar(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),ones(1,64)*0.025,'k',0.05);
shadedErrorBar(meanvaluebin(latB2),nanmean(reshape(trend_olr,72,64),1),nanmean(reshape(trend_olr_err,72,64),1),'g',0.05);
shadedErrorBar(meanvaluebin(latB2),nanmean(reshape(trend_olr_airsL3,72,64),1),nanmean(reshape(trend_olr_airsL3_err,72,64),1),'b',0.05);
hold off
plotaxis2; 
legend('CERES allsky','CERES clrsky','Sergio AIRSL3 profile trends -> OLR trends','Sergio ecRad AIRSL3 timeseries','Sergio direct from AIRS L3','location','best','fontsize',10);
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr');

