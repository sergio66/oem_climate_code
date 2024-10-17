addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE_Git/SHOWSTATS
addpath /umbc/xfs2/strow/asl/matlib/maps/

if ~exist('YY')
  do_XX_YY_from_X_Y
end

if ~exist('llsmap5')
  load llsmap5
end

if ~exist('era5_tzrate')
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_stemprate');
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_ozrate');
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_tzrate');
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_wvrate');
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_rhrate');
  load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat','era5_plays');
end

if ~exist('iNavgLay')
  iNavgLay = floor(100/iNumLay);
end

[mmm,nnn] = size(resultsT);
for ixx = 1 : nnn
  junk = (1:iNavgLay) + (ixx-1)*iNavgLay;
  era5_ozrateNNN(:,ixx) = nanmean(era5_ozrate(:,junk),2);
  era5_rhrateNNN(:,ixx) = nanmean(era5_rhrate(:,junk),2);
  era5_wvrateNNN(:,ixx) = nanmean(era5_wvrate(:,junk),2);
  era5_tzrateNNN(:,ixx) = nanmean(era5_tzrate(:,junk),2);
end

disp(' ')
disp('SKT')
figure(60); clf;
figure(61); clf;
figure(62); clf;
aslmap(60,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1),   [-90 +90],[-180 +180]); title('dST/dt AIRS\_RT');     colorbar; caxis([-1 +1]*0.101); colormap(llsmap5)
aslmap(61,rlat65,rlon73,smoothn((reshape(era5_stemprate',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt ERA5');         colorbar; caxis([-1 +1]*0.101); colormap(llsmap5)
figure(62); [Y,I] = sort(era5_stemprate); plot(era5_stemprate(I),results(I,6),'b.',era5_stemprate,era5_stemprate,'k'); 
xlabel('ERA5 dSK/dt'); ylabel('UMBC dSK/dt'); axis([-1 +1 -1 +1]*0.2)

nanlinearcorrelation(era5_stemprate,results(:,6)');
wahaE = era5_stemprate; wahaU = results(:,6);
latlat = ones(72,1) * rlat'; latlat = cos(pi*latlat/180);
%fprintf(1,'ALL SKT : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'          UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'ALL SKT : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',m,s)

latlat = ones(72,1) * rlat'; 
oo = find(abs(latlat) < 60);
latlat = cos(pi*latlat/180);
wahaE = wahaE(oo); wahaU = wahaU(oo); latlat = latlat(oo);
%fprintf(1,'TRP/MIDLAT SKT : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'          UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'TRP/MIDLAT SKT : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',m,s)

figure(62); clf
clear plotoptions
  plotoptions.cx = [-1 +1]*0.101;
  plotoptions.cmap = llsmap5;
  plotoptions.maintitle = ' ';
aslmap_2x1tiledlayout(era5_stemprate',results(:,6)',62,plotoptions);
%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('WVfrac')

latlat = ones(100,1) * rlat';
latlat = cos(latlat * pi/180);
figure(63);  clf; wahaU = squeeze(nanmean(reshape(fracWV,101,72,64),2)); wahaU = wahaU(1:100,:); pcolor(rlat,era5_plays,wahaU); 
  title('d WVfrac/dt AIRS\_RT'); colorbar; caxis([-1 +1]*0.0101); colormap(llsmap5); set(gca,'ydir','reverse'); ylim([100 1000]); shading interp;
figure(64);  clf; wahaE = squeeze(nanmean(reshape(era5_wvrate,100,72,64),2)); wahaE = wahaE(1:100,:); pcolor(rlat,era5_plays,wahaE); 
  title('d WVfrac/dt ERA5'); colorbar; caxis([-1 +1]*0.0101); colormap(llsmap5); set(gca,'ydir','reverse'); ylim([100 1000]); shading interp;

figure(65); clf; plot(wahaE(:),wahaU(:),'b.',wahaE(:),wahaE(:),'k');
  xlabel('ERA5 dfracWV/dt'); ylabel('UMBC dfracWV/dt'); axis([-1 +1 -1 +1]*0.02)
nanlinearcorrelation(wahaE(:),wahaU(:));
%fprintf(1,'ALL WV : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'ALL WV : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',m,s)
[~,~,~,nmeanWV1,nstdWV1] = myhist2d(wahaE(:),wahaU(:),linspace(-1*0.02,+1*0.02,250),linspace(-1*0.02,+1*0.02,250)); 
    xlabel('ERA5 dfracWV/dt'); ylabel('UMBC dfracWV/dt'); axis([-1 +1 -1 +1]*0.02/2); title('ALL'); disp('ret to continue'); pause

junkH = find(abs(era5_plays) > 300 & abs(era5_plays) < 800);
junkL = find(abs(rlat) < 60);
wahaE = wahaE(junkH,junkL); wahaE = wahaE(:);
wahaU = wahaU(junkH,junkL); wahaU = wahaU(:);
latlat = latlat(junkH,junkL); latlat = latlat(:);
nanlinearcorrelation(wahaE(:),wahaU(:));
%fprintf(1,'TRP/MIDlat WV : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'              : UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'TRP/MIDlat WV : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'              : UMBC trend = %8.6f +/- %8.6f /yr \n',m,s);
[~,~,~,nmeanWV2,nstdWV2] = myhist2d(wahaE(:),wahaU(:),linspace(-1*0.02,+1*0.02,250),linspace(-1*0.02,+1*0.02,250));
    xlabel('ERA5 dfracWV/dt'); ylabel('UMBC dfracWV/dt'); axis([-1 +1 -1 +1]*0.02/2); title('+/- 60 lat, 300-800 mb'); disp('ret to continue'); pause

figure(65); clf
wahaU = squeeze(nanmean(reshape(fracWV,101,72,64),2)); wahaU = wahaU(1:100,:);
wahaE = squeeze(nanmean(reshape(era5_wvrate,100,72,64),2)); wahaE = wahaE(1:100,:);
clear plotoptions
  plotoptions.cx = [-1 +1]*0.0101;
  plotoptions.cmap = llsmap5;
  %plotoptions.maintitle = ' ';
  %plotoptions.maintitle = 'dfracWV(z,lat)/dt'; 
  plotoptions.yLimits = [100 1000]; plotoptions.xLimits = [-85 +85];
  plotoptions.yReverseDir = +1;
  plotoptions.yLinearOrLog = +1;
  plotoptions.ystr = 'Pressure [mb]';
  plotoptions.xstr = 'Latitude [deg]';
profile_plots_2x1tiledlayout(rlat,era5_plays,wahaE,wahaU,65,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Tz')

latlat = ones(100,1) * rlat';
latlat = cos(latlat * pi/180);
figure(66);  clf; wahaU = squeeze(nanmean(reshape(deltaT,101,72,64),2)); wahaU = wahaU(1:100,:); pcolor(rlat,era5_plays,wahaU); 
  title('dT/dt AIRS\_RT'); colorbar; caxis([-1 +1]*0.101); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp;
figure(67);  clf; wahaE = squeeze(nanmean(reshape(era5_tzrate,100,72,64),2)); wahaE = wahaE(1:100,:); pcolor(rlat,era5_plays,wahaE); 
  title('dT/dt ERA5'); colorbar; caxis([-1 +1]*0.101); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp;

figure(68); clf; plot(wahaE(:),wahaU(:),'b.',wahaE(:),wahaE(:),'k');
  xlabel('ERA5 dT/dt'); ylabel('UMBC dT/dt'); axis([-1 +1 -1 +1]*0.2)
nanlinearcorrelation(wahaE(:),wahaU(:));
%fprintf(1,'ALL Tz : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'ALL Tz : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'         UMBC trend = %8.6f +/- %8.6f /yr \n',m,s)
[~,~,~,nmeanT1,nstdT1] = myhist2d(wahaE(:),wahaU(:),linspace(-1*0.2,+1*0.2,100),linspace(-1*0.2,+1*0.2,250)); 
    xlabel('ERA5 dT/dt'); ylabel('UMBC dT/dt'); axis([-1 +1 -1 +1]*0.2/2); title('ALL'); disp('ret to continue'); pause

junkH = find(abs(era5_plays) > 300 & abs(era5_plays) < 800);
junkH = find(abs(era5_plays) > 050 & abs(era5_plays) < 900);
junkL = find(abs(rlat) < 60);
wahaE = wahaE(junkH,junkL); wahaE = wahaE(:);
wahaU = wahaU(junkH,junkL); wahaU = wahaU(:);
nanlinearcorrelation(wahaE(:),wahaU(:));
latlat = latlat(junkH,junkL); latlat = latlat(:);
%fprintf(1,'TRP/MIDlat Tz : ERA5 trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaE(:)),nanstd(wahaE(:)))
%fprintf(1,'              : UMBC trend = %8.6f +/- %8.6f /yr \n',nanmean(wahaU(:)),nanstd(wahaU(:)))
[m,s,m0,s0] = weighted_mean_stddev(wahaE(:),latlat(:));
fprintf(1,'TRP/MIDLAT Tz : ERA5 trend = %8.6f +/- %8.6f /yr \n',m,s)
[m,s,m0,s0] = weighted_mean_stddev(wahaU(:),latlat(:));
fprintf(1,'                UMBC trend = %8.6f +/- %8.6f /yr \n',m,s)
[~,~,~,nmeanT2,nstdT2] = myhist2d(wahaE(:),wahaU(:),linspace(-1*0.2,+1*0.2,250),linspace(-1*0.2,+1*0.2,250)); 
    xlabel('ERA5 dT/dt'); ylabel('UMBC dT/dt'); axis([-1 +1 -1 +1]*0.2/2); title('+/- 60 lat, 050-900 mb'); disp('ret to continue'); pause

figure(68); clf
wahaU = squeeze(nanmean(reshape(deltaT,101,72,64),2)); wahaU = wahaU(1:100,:);
wahaE = squeeze(nanmean(reshape(era5_tzrate,100,72,64),2)); wahaE = wahaE(1:100,:);
clear plotoptions
  plotoptions.cx = [-1 +1]*0.101;
  plotoptions.cmap = llsmap5;
  %plotoptions.maintitle = ' ';
  %plotoptions.maintitle = 'dT(z,lat)/dt'; 
  plotoptions.yLimits = [10 1000]; plotoptions.xLimits = [-85 +85];
  plotoptions.yReverseDir = +1;
  plotoptions.yLinearOrLog = -1;
  plotoptions.ystr = 'Pressure [mb]';
  plotoptions.xstr = 'Latitude [deg]';
profile_plots_2x1tiledlayout(rlat,era5_plays,wahaE,wahaU,68,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S/N
show_S_N

figure(69); clf
  sn = squeeze(nanmean(reshape(resultsT'./resultsTunc',49,72,64),2)); pcolor(rlat,pavg,abs(sn)); colormap jet; colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('retrieval S/N T')
figure(70); clf
  sn = squeeze(nanmean(reshape(resultsWV'./resultsWVunc',49,72,64),2)); pcolor(rlat,pavg,abs(sn)); colormap jet; colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); title('retrieval S/N WV')

figure(69); clf;
  wah = abs(resultsT'./resultsTunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); title('S/N : T/\sigma T'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp; 
  set(gca,'yscale','log'); ylim([1 1000])
figure(70); clf;
  wah = abs(resultsWV'./resultsWVunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); title('S/N : WV/\sigma WV'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp

figure(69); clf;
  wah = abs(resultsTunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); set(gca,'ydir','reverse'); ylim([1 1000]); colormap jet; shading interp; % title('\sigma T'); ; plotaxis2;
  set(gca,'yscale','log'); ylim([10 1000])
  ylim([10 1000]); caxis([0 0.02]);  xlabel('Latitude [deg]'); ylabel('Pressure [mb]'); 
  % colorbar('vertical'); text(90,08,'K/yr'); set(gca,'fontsize',12)
  colorbar('horizontal'); text(-90,2000,'K/yr'); set(gca,'fontsize',12)

figure(70); clf;
  wah = abs(resultsWVunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); set(gca,'ydir','reverse'); ylim([100 1000]); colormap jet; shading interp; % title('\sigma WV'); plotaxis2; 
  set(gca,'yscale','linear'); ylim([100 1000])
  ylim([100 1000]); caxis([0 0.003]); xlabel('Latitude [deg]'); ylabel('Pressure [mb]'); 
  % colorbar('vertical'); text(90,50,'1/yr'); set(gca,'fontsize',12)
  colorbar('horizontal'); text(-100,1250,'1/yr'); set(gca,'fontsize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%

figure(71); clf;
wahaU = squeeze(nanmean(reshape(deltaT,101,72,64),2)); wahaU = wahaU(1:100,:); 
wahaE = squeeze(nanmean(reshape(era5_tzrate,100,72,64),2)); wahaE = wahaE(1:100,:); 
pcolor(rlat,era5_plays,abs(wahaE-wahaU)); 
  title('dT/dt ERA5-UMBC'); colorbar; caxis([0 +1]*0.02); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp;

figure(72); clf;
wah = abs(resultsTunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
for ii = 1 : 64
  wahx(:,ii) = interp1(log(pavg),wah(:,ii),log(era5_plays),[],'extrap');
end
pcolor(rlat,era5_plays,abs(wahaE-wahaU)./wahx); 
  title('dT/dt (ERA5-UMBC)/uncT'); colorbar; caxis([0 +1]*5); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp;

%%%%%

figure(73); clf
wahaU = squeeze(nanmean(reshape(fracWV,101,72,64),2)); wahaU = wahaU(1:100,:);
wahaE = squeeze(nanmean(reshape(era5_wvrate,100,72,64),2)); wahaE = wahaE(1:100,:); 
pcolor(rlat,era5_plays,abs(wahaE-wahaU)); 
  title('dfracWV/dt ERA5-UMBC'); colorbar; caxis([0 +1]*0.003); colormap(jet); set(gca,'ydir','reverse'); ylim([100 1000]); shading interp;

figure(74); clf;
wah = abs(resultsWVunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
for ii = 1 : 64
  wahx(:,ii) = interp1(log(pavg),wah(:,ii),log(era5_plays),[],'extrap');
end
pcolor(rlat,era5_plays,abs(wahaE-wahaU)./wahx); 
  title('dfracWV/dt (ERA5-UMBC)/uncWV'); colorbar; caxis([0 +1]*5); colormap(jet); set(gca,'ydir','reverse'); ylim([100 1000]); shading interp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i400 = find(era5_plays >= 400,1);
figure(75); clf; 
aslmap(75,rlat65,rlon73,smoothn((reshape(fracWV(i400,:)',72,64)') ,1),     [-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb'); 
figure(76); clf; 
aslmap(76,rlat65,rlon73,smoothn((reshape(era5_wvrate(i400,:)',72,64)') ,1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5); % title('dfracWV/dt ERA5 400 mb');    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iOffset = 80;
statistical_significance_sergio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(62); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/compare_era5_vs_calc_retrieval_dsst_dt');
figure(65); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/compare_era5_vs_calc_retrieval_dwvfrac_dt');
figure(68); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/compare_era5_vs_calc_retrieval_dtz_dt');

%% Xyou need to make sure you have eg  load /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS.mat   %% this is for Fig 14 of trends paper (observation)
%% Xor from ../FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%% X  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MLS.mat';    iNumYears = 20;     %% use CarbonTracker CO2 trends, MLS a priori
%% X  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';        iNumYears = 20;     %% use CarbonTracker CO2 trends
%% X  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat'; iNumYears = 20;     %% use CarbonTracker CO2 trends
%% X                          vs   eg  load /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat   %% this is for Fig  7 of trends paper (simulation)
%% X 
%% Xsee FIND_NWP_MODEL_TRENDS/driver_compare_trends_Day_vs_Night.m --> FIND_NWP_MODEL_TRENDS/get_umbc_day_night_name.m

%% << driver_trendspaper_figure15_uncertainty_ztest.m >> 
%% 
%% figure(69); set(gca,'fontsize',12); text(-85,20,'(A)','fontsize',12);  sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/unc_dtz_dt');
%% figure(70); set(gca,'fontsize',12); text(-85,150,'(C)','fontsize',12); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/unc_dwvfrac_dt');
%% figure(75); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/umbc_dwvfrac_dt_400mb');
%% figure(76); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/era5_dwvfrac_dt_400mb');
%%
%% figure(86); set(gca,'fontsize',12); text(-85,20,'(B)','fontsize',12); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/zvalue_dtz_dt');
%% figure(88); set(gca,'fontsize',12); % sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/zvalue_drh_dt');
%% figure(90); set(gca,'fontsize',12); text(-85,150,'(D)','fontsize',12); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/zvalue_dwvfrac_dt');

%} 
