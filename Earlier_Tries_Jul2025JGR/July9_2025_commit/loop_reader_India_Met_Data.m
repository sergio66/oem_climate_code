iConsecutive = 0;

iTile = input('Enter which tile in India to look at : ');
if length(iTile) == 0
  iTile = 2788;
end

figure(1); clf

iMaxOrMin = input('Enter (+1/default) for max temp (-1) for min temp (0) for rainfall: ');
if length(iMaxOrMin) == 0
  iMaxOrMin = +1;
end
for yy = 2002 : 2022
  if iMaxOrMin > 0
    fin = ['/asl/s1/sergio/NSD_India_Met_Data/Maxtemp_MaxT_' num2str(yy) '.GRD'];
    iMax = 31;
    jMax = 31;
    imd_lat = linspace(07.5,37.5,31);
    imd_lon = linspace(67.5,97.5,31);
  elseif iMaxOrMin < 0
    fin = ['/asl/s1/sergio/NSD_India_Met_Data/Mintemp_MinT_' num2str(yy) '.GRD'];
    iMax = 31;
    jMax = 31;
    imd_lat = linspace(07.5,37.5,31);
    imd_lon = linspace(67.5,97.5,31);
  elseif iMaxOrMin == 0
    fin = ['/asl/s1/sergio/NSD_India_Met_Data/Rainfall_ind' num2str(yy) '_rfp25.grd'];
    iMax = 135;
    jMax = 129;
    imd_lat = linspace(06.5,038.5,129);
    imd_lon = linspace(66.5,100.0,135);
  end

  fprintf(1,'%4i : %s \n',yy,fin);
  T = reader_India_Met_Data(fin,iMaxOrMin);
  if abs(iMaxOrMin) == 1
    T(T > 99) = NaN;
    T(T < 00) = NaN;
  else
    T(T < 0) = NaN;
  end
  if mod(yy,4) == 0
    daysINmonth = [31 29 31 30 31 30 31 31 30 31 30 31];
  else
    daysINmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
  end
  iCnt = 0;
  days = 0;
  for mm = 1 : 12
    iCnt = iCnt + 1;
    days = (1:daysINmonth(mm)) + days(end);
    T12(iCnt,:,:) = squeeze(nanmean(T(days,:,:),1));

    iConsecutive = iConsecutive + 1;
    Tconsecutive(iConsecutive,:,:) = squeeze(nanmean(T(days,:,:),1));;

    yysave(iConsecutive) = yy;
    mmsave(iConsecutive) = mm;
    daysSince2002IMD(iConsecutive) = change2days(yy,mm,15,2002);
  end
  Tall(yy-2002+1,:,:,:) = T12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
whos T*
for mm = 1 : 12
  junk = squeeze(Tall(:,mm,:,:));
  junk = squeeze(nanmean(junk,1));
  pcolor(junk); colormap jet; caxis([20 40]); colorbar; 
  if abs(iMaxOrMin) == 1
    caxis([20 40]); 
  end
  title(['mm = ' num2str(mm,'%02d')]);; shading interp
  pause(0.25);
end

figure(1); clf;
[YIMD,XIMD] = meshgrid(imd_lat,imd_lon);
moo = load('coast.mat');

pcolor(XIMD,YIMD,squeeze(nanmean(Tconsecutive,1))'+273.15); colormap jet; caxis([20 40]+273.15); colorbar; title('Mean IMD Data read in, 20 years'); shading flat
pcolor(XIMD,YIMD,squeeze(nanmean(Tconsecutive,1))');        colormap jet; caxis([20 40]);        colorbar; title('Mean IMD Data read in, 20 years'); shading flat
hold on; plot(moo.long,moo.lat,'k','linewidth',2); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tS = find(yysave == 2002 & mmsave == 9);
tE = find(yysave == 2022 & mmsave == 8);
timeSE = daysSince2002IMD(tS:tE);
warning off
for jj = 1 : jMax
  for ii = 1 : iMax
    data = squeeze(Tconsecutive(tS:tE,jj,ii));
    boo = find(isfinite(data));
    if length(boo) > 10      
      [B,err] = Math_tsfit_lin_robust(timeSE(boo),data(boo),4);
      trend(jj,ii) = B(2);
      trend_err(jj,ii) = err.se(2);
    else
      trend(jj,ii) = NaN;
      trend_err(jj,ii) = NaN;
    end
  end
end
warning on
load llsmap5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
pcolor(XIMD,YIMD,trend'); title('Trend K/yr'); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; shading interp; shading flat
figure(2); ax = axis;
hold on; plot(moo.long,moo.lat,'k','linewidth',2); hold off
axis(ax); 
hold on; plot(imd_lon,ones(1,iMax)*15,'k',imd_lon,20*nanmin(trend)+15,'b',imd_lon,ones(1,iMax)*25,'k',imd_lon,20*nanmean(trend)+25,'g',imd_lon,ones(1,iMax)*35,'k',imd_lon,20*nanmax(trend)+35,'r','linewidth',2); hold off
if iMaxOrMin > 0
  title('MaxT Trend K/yr');
elseif iMaxOrMin < 0
  title('MinT Trend K/yr');
elseif iMaxOrMin == 0
  title('Rainfall Trend K/yr');
end

for ii = 1 : 2
  figure(ii); hold on; plot(moo.long,moo.lat,'k','linewidth',2); hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(iMaxOrMin) == 1
  do_XX_YY_from_X_Y

  obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat','b_asc','b_desc');
  if iMaxOrMin > 0
    era5spectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat','trend');
    obsspectra.x = (reshape(obsspectra.b_asc,4608,2645))';
  else
    era5spectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat','trend');
    obsspectra.x = (reshape(obsspectra.b_desc,4608,2645))';
  end
  figure(3); clf; plot_72x64_tiles([],obsspectra.x(1520,:)); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('AIRS BT1231 trends')
  figure(4); clf; plot_72x64_tiles([],era5spectra.trend(1520,:)); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('ERA5 BT1231 trends')
  figure(3); clf; scatter_coast(X(:),Y(:),100,obsspectra.x(1520,:)); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('AIRS BT1231 trends')
  figure(4); clf; scatter_coast(X(:),Y(:),100,era5spectra.trend(1520,:)); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('ERA5 BT1231 trends')
  %figure(3); clf; simplemap(Y(:),X(:),obsspectra.x(1520,:)',5); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('AIRS BT1231 trends')
  %figure(4); clf; simplemap(Y(:),X(:),era5spectra.trend(1520,:)',5); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('ERA5 BT1231 trends')

  F1 = griddedInterpolant(X,Y,reshape(obsspectra.x(1520,:),72,64));
  F2 = griddedInterpolant(X,Y,reshape(era5spectra.trend(1520,:),72,64));
  junk1 = F1(XIMD,YIMD);
  junk2 = F2(XIMD,YIMD);
  figure(5); clf; simplemap(YIMD,XIMD,junk1); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('AIRS BT1231 trends')
  figure(6); clf; simplemap(YIMD,XIMD,junk2); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('ERA5 BT1231 trends')

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% see ../FIND_NWP_MODEL_TRENDS/prep_colWV_T_WV_trends_Day_vs_Night.m
  airsL3_day_file   = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
  airsL3_night_file = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

  climcapsL3_day_file   = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
  climcapsL3_night_file = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

  if iMaxOrMin > 0
    era5trend = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc_surf.mat','trend_stemp');
    obsspectra = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day_removeemisstrends.mat','results');
    obsspectra.x = obsspectra.results(:,6);
    airsL3 = load(airsL3_day_file,'thestats64x72');
    climL3 = load(climcapsL3_day_file,'thestats64x72');
  else
    era5trend = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_surf.mat','trend_stemp');
    obsspectra = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_removeemisstrends.mat','results');
    obsspectra.x = obsspectra.results(:,6);
    airsL3 = load(airsL3_night_file,'thestats64x72');
    climL3 = load(climcapsL3_night_file,'thestats64x72');
  end

  F1 = griddedInterpolant(X,Y,reshape(obsspectra.x,72,64));
  F2 = griddedInterpolant(X,Y,reshape(era5trend.trend_stemp,72,64));
  F3A = griddedInterpolant(X,Y,airsL3.thestats64x72.stemprate);
  F3C = griddedInterpolant(X,Y,climL3.thestats64x72.stemprate);
  junk1  = F1(XIMD,YIMD);
  junk2  = F2(XIMD,YIMD);
  junk3A = F3A(XIMD,YIMD);
  junk3C = F3C(XIMD,YIMD);
  figure(7); clf; simplemap(YIMD,XIMD,junk1);  axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('UMBC     STEMP trends')
  figure(8); clf; simplemap(YIMD,XIMD,junk2);  axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('ERA5     STEMP trends')
  figure(9); clf; simplemap(YIMD,XIMD,junk3A); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('AIRS L3  STEMP trends')
  figure(10);clf; simplemap(YIMD,XIMD,junk3C); axis(ax); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; title('CLIMCAPS STEMP trends')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  for ii = 3 : 10
    figure(ii); hold on; plot(moo.long,moo.lat,'k','linewidth',2); hold off;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(iMaxOrMin) > 0
  % more /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/ReadmeQuick
  % (39-1)*72 + 52
  % iTile = 2788;
      indY = floor(iTile/72);
      indX = iTile - indY*72;
      indY = indY + 1;
  dir0 = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';
  trend_data = load([dir0 '/LatBin' num2str(indY,'%02i') '/LonBin' num2str(indX,'%02i') '/iQAX_3_summarystats_LatBin' num2str(indY,'%02i') '_LonBin' num2str(indX,'%02i') '_timesetps_001_457_V1.mat']);
  figure(11); vertices = plot_72x64_tiles(iTile);

  daysSince2002 = change2days(trend_data.year_asc,trend_data.month_asc,trend_data.day_asc,2002);

  for ii = 3 : 10
    figure(ii); rectangle('position',[vertices(1) vertices(3) vertices(2)-vertices(1) vertices(4)-vertices(3)],'EdgeColor','k','linewidth',2);
  end

  [mooX,mooY] = find(XIMD >= vertices(1) & XIMD <= vertices(2) & YIMD >= vertices(3) & YIMD <= vertices(4));
  [moo] = find(XIMD >= vertices(1) & XIMD <= vertices(2) & YIMD >= vertices(3) & YIMD <= vertices(4));
  verticesIMD = [min(XIMD(moo)) max(XIMD(moo)) min(YIMD(moo)) max(YIMD(moo))];
  for ii = 1 : 2
    figure(ii); rectangle('position',[verticesIMD(1) verticesIMD(3) verticesIMD(2)-verticesIMD(1) verticesIMD(4)-verticesIMD(3)],'EdgeColor','k','linewidth',2);
  end

  figure(9); clf;   
  nSmoothYears = 4;
  nSmoothYears = 1;
  nSmoothYears = 2;
  wonkIMD = squeeze(Tconsecutive(:,mooX,mooY)); wonkIMD = nanmean(squeeze(nanmean(wonkIMD,2)),2);
  figure(9); plot(2002+daysSince2002/365,trend_data.quantile1231_asc(:,1)-273.15,2002+daysSince2002/365,smooth(trend_data.quantile1231_asc(:,1),23*2)-273.15)
  figure(9); plot(2002+daysSince2002/365,smooth(trend_data.quantile1231_asc(:,3),23*nSmoothYears)-273.15,'rx-',2002+daysSince2002/365,smooth(trend_data.quantile1231_desc(:,3),23*nSmoothYears)-273.15,'bo-',...
                  2002+daysSince2002/365,smooth(trend_data.meanBT_asc(:,1520),23*nSmoothYears)-273.15,'m',2002+daysSince2002/365,smooth(trend_data.meanBT_desc(:,1520),23*nSmoothYears)-273.15,'c',...
                  2002+daysSince2002IMD/365,smooth(wonkIMD,12*nSmoothYears),'k','linewidth',2); title('quantile1231 Q03 + MeanBT1231')
  hl = legend('Q03 ASC','Q03 DESC','mean ASC','mean DESC','IMD','location','best','fontsize',8); ylabel('T [deg C]'); xlabel('Time')
  
  if iMaxOrMin > 0
    figure(9); plot(2002+daysSince2002/365,smooth(trend_data.quantile1231_asc(:,3),23*nSmoothYears)-273.15,'rx-',2002+daysSince2002IMD/365,smooth(wonkIMD,12*nSmoothYears),'k','linewidth',2); title('quantile1231 Q03 + MeanBT1231')
    hl = legend('Q03 ASC','IMD','location','best','fontsize',8); ylabel('T [deg C]'); xlabel('Time')
  else
    figure(9); plot(2002+daysSince2002/365,smooth(trend_data.quantile1231_desc(:,3),23*nSmoothYears)-273.15,'bx-',2002+daysSince2002IMD/365,smooth(wonkIMD,12*nSmoothYears),'k','linewidth',2); title('quantile1231 Q03 + MeanBT1231')
    hl = legend('Q03 DESC','IMD','location','best','fontsize',8); ylabel('T [deg C]'); xlabel('Time')
  end

  figure(10); clf; plot(2002+daysSince2002/365,trend_data.meanBT_asc(:,1520),2002+daysSince2002/365,trend_data.meanBT_desc(:,1520),'linewidth',2); hl = legend('asc','desc'); title('BT1231 mean')
  figure(10); clf; plot(2002+daysSince2002/365,trend_data.quantile1231_asc,2002+daysSince2002/365,trend_data.meanBT_asc(:,1520),'kx-','linewidth',2); 
    hl = legend('1','2','3','4','5','mean','location','best','fontsize',10);; title('BT1231 Q01-05 and mean')
  figure(10); clf; plot(2002+daysSince2002/365,trend_data.quantile1231_asc - 273.15,2002+daysSince2002/365,trend_data.meanBT_asc(:,1520) - 273.15,'kx-','linewidth',2); 
    hl = legend('1','2','3','4','5','mean','location','best','fontsize',10);; title('BT1231 Q01-05 and mean')
  disp('Mean 1231 asc quantile (in deg C'); printarray(nanmean(trend_data.quantile1231_asc,1)-273.15)

  jett = jet(256); jett(1,:) = 1;
  figure(11); pcolor(2002+daysSince2002/365,trend_data.dbt,trend_data.hist1231_desc'); shading flat; colorbar; title('BT1231 DESC')
  figure(12); pcolor(2002+daysSince2002/365,trend_data.dbt,trend_data.hist1231_asc'); shading flat; colorbar; title('BT1231 DESC')
  figure(11); pcolor(2002+daysSince2002/365,trend_data.dbt-273.15,log10(trend_data.hist1231_desc')); shading flat; colorbar; title('BT1231 DESC'); colormap(jett); caxis([-3.5 -1.5]); xlabel('Time'); ylabel('BT1231 deg C')
    hold on; plot(2002+daysSince2002/365,trend_data.quantile1231_desc(:,3)-273.15,'kx-','linewidth',2); hold off
  figure(12); pcolor(2002+daysSince2002/365,trend_data.dbt-273.15,log10(trend_data.hist1231_asc'));  shading flat; colorbar; title('BT1231 ASC'); colormap(jett); caxis([-3.5 -1.5]); xlabel('Time'); ylabel('BT1231 deg C')
    hold on; plot(2002+daysSince2002/365,trend_data.quantile1231_asc(:,3)-273.15,'kx-','linewidth',2); hold off
end
