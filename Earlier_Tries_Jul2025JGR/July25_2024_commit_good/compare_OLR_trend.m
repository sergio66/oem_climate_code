addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools

[~,~,pjunk,~] = rtpread('summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp');

load llsmap5
load latB64.mat

iNumYears = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCorT = +1;  %% what I started with
iCorT = -1;  %% Ryan suggestion

iCorT

junkC = load(['ceres_trends_' num2str(iNumYears,'%02d') 'year_C.mat']); %% original of mine
junkT = load(['ceres_trends_' num2str(iNumYears,'%02d') 'year_T.mat']); %% ryan suggestion
if iCorT < 0
  junk = junkT;
else
  junk = junkC;
end
ceres_trend = junk.ceres_trend;
z11A = junkT.ceres_trend.trend_toa_lw_all_4608;
z12A = junkC.ceres_trend.trend_toa_lw_all_4608;
  plotoptions2.str11 = 'CERES C';   plotoptions2.str12 = 'CERES T';
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'ALLSKY dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.5; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_1x2tiledlayout(z11A,z12A,1,plotoptions2);
  diffALL = z11A-z12A;
z11C = junkT.ceres_trend.trend_toa_lw_clr_4608;
z12C = junkC.ceres_trend.trend_toa_lw_clr_4608;
  plotoptions2.str11 = 'CERES C';   plotoptions2.str12 = 'CERES T';
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'CLRSKY dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.5; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_1x2tiledlayout(z11C,z12C,2,plotoptions2);
  diffCLR = z11C-z12C;
z11 = diffALL;
z12 = diffCLR;
  plotoptions2.str11 = 'CERES ALL';   plotoptions2.str12 = 'CERES CLR';
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'T-C dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.5; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_1x2tiledlayout(z11,z12,3,plotoptions2);
figure(4); 
  plot(meanvaluebin(latB2),nanmean(reshape(z11A,72,64),1),'r',meanvaluebin(latB2),nanmean(reshape(z12A,72,64),1),'m',...
       junkT.ceres_trend.trend_lat,junkT.ceres_trend.trend_lw,'r--',junkC.ceres_trend.trend_lat,junkC.ceres_trend.trend_lw,'m--',...
       meanvaluebin(latB2),nanmean(reshape(z11C,72,64),1),'b',meanvaluebin(latB2),nanmean(reshape(z12C,72,64),1),'c',...
       junkT.ceres_trend.trend_lat,junkT.ceres_trend.trend_lw_clr,'b--',junkC.ceres_trend.trend_lat,junkC.ceres_trend.trend_lw_clr,'c--',...
       'linewidth',2);
  plotaxis2; 
  hl = legend('ALLSKY 72x64 T','ALLSKY 72x64 C','ALLSKY 180 T','ALLSKY 180 C','CLRSKY 72x64 T','CLRCKY 72x64 C','CLRSKY 180 T','CLRSKY 180 C',...
                         'location','best','fontsize',10);
  xlim([-90 +90])
  xlabel('Latitude / [deg]'); ylabel('OLR trend/ [W/m2/yr]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is from taking 12 months x 20 years ptofiles, and running through ecRad
%% see ../FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_OLR.m
%% see ../FIND_NWP_MODEL_TRENDS/computeERA5_OLR_trend.m
junk = load('../FIND_NWP_MODEL_TRENDS/OLR_ecRad/ERA5/all_era5_olr.mat');
  ecRad_ERA5_20years = junk;
  junk4608.era5_ecRad = junk.trend_olr;

%% this is from downloading ERA5 TOA OLR
junk = load('era5_monthly_olrtrends_directcomputation_20years.mat');
  directERA5 = junk;

plot(meanvaluebin(latB2),nanmean(reshape(directERA5.trend_olr_clr,72,64),1),'rx-',...
     meanvaluebin(latB2),nanmean(reshape(ecRad_ERA5_20years.trend_olr_ERA5,72,64),1),'m',meanvaluebin(latB2),nanmean(reshape(ecRad_ERA5_20years.trend_olr,72,64),1),'b','linewidth',2)
     plotaxis2; title('ERA5 OLR trends'); hl = legend('direct from their OLR timeseries','direct from their OLR timeseries (again)','ERA5 profiles --> ecRad --> trends','location','best');
plot(meanvaluebin(latB2),nanmean(reshape(directERA5.trend_olr_clr,72,64),1),'rx-',meanvaluebin(latB2),nanmean(reshape(ecRad_ERA5_20years.trend_olr,72,64),1),'b','linewidth',2)
  plotaxis2; title('ERA5 OLR trends'); hl = legend('direct from their OLR timeseries','ERA5 profiles --> ecRad --> trends','location','best');

%% this is from computing ERA5 T(z),WV(z),SKT trends, then computing OLR trends
junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'era5_spectral_olr');
  era5_spectral_olr = junk.era5_spectral_olr;

iFig = 5;
figure(iFig); clf
miaow = era5_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - era5_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'k','linewidth',2);
  junk64.era5 = nanmean(reshape(miaow,72,64),1);
  junk4608.era5 = miaow;
miaow = directERA5.trend_olr_clr; 
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'bx-','linewidth',2); 
  junk64.era5_direct = nanmean(reshape(miaow,72,64),1);
  junk4608.era5_direct = miaow;
miaow = ecRad_ERA5_20years.trend_olr;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'g','linewidth',2); 
  junk4608.era5_ecRad = miaow;
hold off; plotaxis2; 
hl = legend('from geophysical trends','from direct OLR in Copernicus','from ecRad runs of monthly profiles','location','best');
title('ERA5 flux trends')

figure(iFig+1); clf
z11 = junk4608.era5_direct;
z12 = junk4608.era5;
z22 = junk4608.era5_ecRad;
  plotoptions2.str11 = 'ERA5 direct OLR'; plotoptions2.str12 = 'ERA5 geophys trends --> OLR';   
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'ERA5 dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.5; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_1x2tiledlayout(z11,z12,iFig+1,plotoptions2);

  plotoptions2.str11 = 'ERA5 direct OLR';   plotoptions2.str12 = 'diff :  direct OLR - from geotrends';
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'ERA5 dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.25; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_1x2tiledlayout(z11,z11-z12,iFig+1,plotoptions2);

  plotoptions2.str11 = 'ERA5 direct OLR';   plotoptions2.str12 = 'ERA5 geotrends --> olr trends';          plotoptions2.str13 = 'diff';
  plotoptions2.str21 = 'ERA5 direct OLR';   plotoptions2.str22 = 'ERA5 monthly --> ecRad --> olr trends';  plotoptions2.str23 = 'diff';
  plotoptions2.ystr = 'Latitude'; plotoptions2.xstr = 'Longitude';
  plotoptions2.maintitle = 'ERA5 dOLR/dt W/m2/yr';
  plotoptions2.cx = [-1 +1]*0.25; plotoptions2.cmap = llsmap5; plotoptions2.yReverseDir = -1; plotoptions2.yLinearOrLog = +1;
  aslmap_2x3tiledlayout(z11,z12,z11-z12,z11,z22,z11-z22,iFig+1,plotoptions2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;

junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iFig+2); clf
plot(ceres_trend.trend_lat,ceres_trend.trend_lw_clr,ceres_trend.trend_lat,ceres_trend.trend_lw,'linewidth',6);
  junk64.ceres = (interp1(ceres_trend.trend_lat,ceres_trend.trend_lw_clr,meanvaluebin(latB2),[],'extrap')');
  if iCorT > 0
    junk4608.ceres    = junkC.ceres_trend.trend_toa_lw_clr_4608;
    junk4608.ceresall = junkC.ceres_trend.trend_toa_lw_all_4608;
  elseif iCorT < 0
    junk4608.ceres    = junkT.ceres_trend.trend_toa_lw_clr_4608;
    junk4608.ceresall = junkT.ceres_trend.trend_toa_lw_all_4608;
  end
miaow = era5_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - era5_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
  junk64.era5 = nanmean(reshape(miaow,72,64),1);
  junk4608.era5 = miaow;
miaow = merra2_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - merra2_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'s-','linewidth',2);
  junk64.merra2 = nanmean(reshape(miaow,72,64),1);
  junk4608.merra2 = miaow;
miaow = umbc_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'x-','linewidth',4);
  junk64.umbc = nanmean(reshape(miaow,72,64),1);
  junk4608.umbc = miaow;
miaow = airsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
  junk64.airsL3 = nanmean(reshape(miaow,72,64),1);
  junk4608.airsL3 = miaow;
miaow = climcapsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - climcapsL3_spectral_olr.olr0_ecRad.clr;;
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'linewidth',2);
  junk64.climcaps = nanmean(reshape(miaow,72,64),1);
  junk4608.climcaps = miaow;
miaow = directERA5.trend_olr_clr; 
  hold on; plot(meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'x-','linewidth',2); 
  junk64.era5_direct = nanmean(reshape(miaow,72,64),1);
  junk4608.era5_direct = miaow;
hold off; plotaxis2; 
hl = legend('CERES clrsky','CERES allsky','ERA5','MERRA2','THIS WORK','AIRSL3','CLIMAPS L3','ERA5 clrsky directly from files','location','best','fontsize',10); 
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

figure(iFig+3); clf;
plot(ceres_trend.trend_lat,ceres_trend.trend_lw_clr,meanvaluebin(latB2),nanmean(reshape(junk4608.ceres,72,64),1),'linewidth',2); 
plotaxis2; hl = legend('from 180 latbins','from 72x64','location','best'); title('CERES OLR clr trends')
plot(ceres_trend.trend_lat,ceres_trend.trend_lw_clr,ceres_trend.trend_lat,smooth(ceres_trend.trend_lw_clr,3),meanvaluebin(latB2),nanmean(reshape(junk4608.ceres,72,64),1),'linewidth',2);
plotaxis2; hl = legend('from 180 latbins','smmoth(180 latbins,3)','from 72x64','location','best'); title('CERES OLR clr trends')

figure(iFig+3); clf; 
plot(meanvaluebin(latB2),[junk64.ceres; junk64.era5; junk64.merra2; junk64.umbc; junk64.airsL3; junk64.climcaps; junk64.era5_direct],'linewidth',2);
plotaxis2; hl = legend('CERES clrsky','ERA5','MERRA2','THIS WORK','AIRSL3','CLIMAPS L3','ERA5 clrsky directly from files','location','best','fontsize',10); 
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

figure(iFig+4); clf;
plot(meanvaluebin(latB2),junk64.era5_direct - [junk64.ceres; junk64.era5; junk64.merra2; junk64.umbc; junk64.airsL3; junk64.climcaps; junk64.era5_direct],'linewidth',2);
plotaxis2; hl = legend('CERES clrsky','ERA5','MERRA2','THIS WORK','AIRSL3','CLIMAPS L3','ERA5 clrsky directly from files','location','best','fontsize',10); 
xlabel('Latitude'); title('\delta Flux Trend : ERA5direct - X'); ylabel('Flux/yr W/m2/yr'); 

disp('\n Comparing to ERA_direct ...')
thediff = junk64.era5_direct - [junk64.ceres; junk64.era5; junk64.merra2; junk64.umbc; junk64.airsL3; junk64.climcaps; junk64.era5_direct];
thediff = sum(thediff.*thediff,2);
fprintf(1,'chisqr with CERES    trends = %8.6f \n',thediff(1));
fprintf(1,'chisqr with ERA5     trends = %8.6f \n',thediff(2));
fprintf(1,'chisqr with MERRA2   trends = %8.6f \n',thediff(3));
fprintf(1,'chisqr with UMBC     trends = %8.6f \n',thediff(4));
fprintf(1,'chisqr with AIRSL3   trends = %8.6f \n',thediff(5));
fprintf(1,'chisqr with CLIMCAPS trends = %8.6f \n',thediff(6));

disp('\n Comparing to CERES ...')
thediff = junk64.ceres - [junk64.ceres; junk64.era5; junk64.merra2; junk64.umbc; junk64.airsL3; junk64.climcaps; junk64.era5_direct];
thediff = sum(thediff.*thediff,2);
fprintf(1,'chisqr with ERA5       trends = %8.6f \n',thediff(2));
fprintf(1,'chisqr with MERRA2     trends = %8.6f \n',thediff(3));
fprintf(1,'chisqr with UMBC       trends = %8.6f \n',thediff(4));
fprintf(1,'chisqr with AIRSL3     trends = %8.6f \n',thediff(5));
fprintf(1,'chisqr with CLIMCAPS   trends = %8.6f \n',thediff(6));
fprintf(1,'chisqr with ERA5direct trends = %8.6f \n',thediff(7));

z11    = junk4608.ceres;
z11all = junk4608.ceresall;
z12 = junk4608.era5;
z21 = junk4608.merra2;
z22 = junk4608.umbc;
z31 = junk4608.airsL3; 
z32 = junk4608.climcaps;
  plotoptions6.str11 = 'CERES';   plotoptions6.str12 = 'ERA5';
  plotoptions6.str21 = 'MERRA2';  plotoptions6.str22 = 'UMBC';
  plotoptions6.str31 = 'AIRS L3'; plotoptions6.str32 = 'CLIMCAPS';
  plotoptions6.ystr = 'Latitude'; plotoptions6.xstr = 'Longitude';
  plotoptions6.maintitle = 'dOLR/dt W/m2/yr';
  plotoptions6.cx = [-1 +1]*0.25; plotoptions6.cmap = llsmap5; plotoptions6.yReverseDir = -1; plotoptions6.yLinearOrLog = +1;
  aslmap_3x2tiledlayout(z11,z12,z21,z22,z31,z32,iFig+5,plotoptions6);
z12x = junk4608.era5_direct;
zall = [z11; z12; z12x; z21; z22; z31; z32];
corrcoef(double(zall'))
corr(double(zall'))

figure(iFig+6); clf
varnames = {'CERES','ERA5','ERA5 direct','MERRA2','UMBC','AIRS L3','CLIMCAPS'};
[R,Pvalue] = corrplot(double(zall'),Varnames = varnames)
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
arrayfun(@(ax) xlim(ax, [-1 +1]*0.5), haxes);

[meanOLR(1) stdOLR(1)] = weighted_mean_stddev(z11,cos(pjunk.rlat*pi/180));
[meanOLR(2) stdOLR(2)] = weighted_mean_stddev(z11all,cos(pjunk.rlat*pi/180));
[meanOLR(3) stdOLR(3)] = weighted_mean_stddev(z12,cos(pjunk.rlat*pi/180));
[meanOLR(4) stdOLR(4)] = weighted_mean_stddev(z21,cos(pjunk.rlat*pi/180));
[meanOLR(5) stdOLR(5)] = weighted_mean_stddev(z22,cos(pjunk.rlat*pi/180));
[meanOLR(6) stdOLR(6)] = weighted_mean_stddev(z31,cos(pjunk.rlat*pi/180));
[meanOLR(7) stdOLR(7)] = weighted_mean_stddev(z32,cos(pjunk.rlat*pi/180));
fprintf(1,'cosine averaged OLR trends : CERESclr = %8.4f   CERESall = %8.6f    ERA5 = %8.4f   MERRA2 = %8.4f    UMBC = %8.4f   airsL3 = %8.4f   climcapsL3 = %8.4f W/m2/yr \n',meanOLR);
fprintf(1,'                    stddev : CERESclr = %8.4f   CERESall = %8.6f    ERA5 = %8.4f   MERRA2 = %8.4f    UMBC = %8.4f   airsL3 = %8.4f   climcapsL3 = %8.4f W/m2/yr \n',stdOLR);

figure(iFig+7); clf; 
plot(meanvaluebin(latB2),nanmean(reshape(z11,72,64),1),meanvaluebin(latB2),nanmean(reshape(z11all,72,64),1),meanvaluebin(latB2),nanmean(reshape(z12,72,64),1),meanvaluebin(latB2),nanmean(reshape(z21,72,64),1),...
     meanvaluebin(latB2),nanmean(reshape(z22,72,64),1),meanvaluebin(latB2),nanmean(reshape(z31,72,64),1),meanvaluebin(latB2),nanmean(reshape(z32,72,64),1),'linewidth',2);
plotaxis2; hl = legend('CERES clrsky','CERES allsky','ERA5','MERRA2','THIS WORK','AIRSL3','CLIMAPS L3','location','best','fontsize',10); 
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iFig+8); clf
axax = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat');
miaow = airsL3_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - airsL3_spectral_olr.olr0_ecRad.clr;;
plot(ceres_trend.trend_lat,ceres_trend.trend_lw,ceres_trend.trend_lat,ceres_trend.trend_lw_clr,...
     meanvaluebin(latB2),nanmean(reshape(miaow,72,64),1),'kx-',...
     meanvaluebin(latB2),nanmean(reshape(axax.thestats64x72_other.clrolrrate,72,64),1),...
     meanvaluebin(latB2),nanmean(reshape(axax.thestats64x72_other.olrrate,72,64),1),'linewidth',2);
plotaxis2; legend('CERES allsky','CERES clrsky','Sergio AIRS L3','Joel AIRS L3 clr','Joel AIRS L3 allsky','location','best','fontsize',10);
xlabel('Latitude'); title('Flux Trend'); ylabel('Flux/yr W/m2/yr'); 

%hmm I did not do CLIMCAPS OLR
%/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat
%axax = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat');
%axax = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
