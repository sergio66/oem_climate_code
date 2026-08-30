function wind = wind_speed_changes();

if ~exist('YY')

  %% run this as stand alone to check anonaly code
  %% run this as stand alone to check anonaly code
  %% run this as stand alone to check anonaly code

  do_XX_YY_from_X_Y
  addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
  addpath /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE_Git/COLORMAP
  addpath /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE_Git/PLOTTER
  addpath /home/sergio/MATLABCODE/TIME
  addpath /asl/matlib/aslutil
end

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/driver_wspeed_trends.m
iJunkModel = input('Enter (+5 [DEFAULT]) ERA5 monthly (+2) MERRA2 : ');
if length(iJunkModel) == 0
  iJunkModel = 5;
end

if iJunkModel == 5
  wind = load('/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/wspeed_2002_09_2024_06.mat');
  junkstr = ' ERA5 ';
else
  wind = load('/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA/merra2_wspeed_2002_09_2024_06.mat');
  junkstr = ' MERRA2 ';
end

[num_months_wind,~] = size(wind.wspeedsave);
junkdoy = (1:num_months_wind)*365.25/12;
junkdoy = change2days(wind.yysave,wind.mmsave,ones(size(wind.yysave))*15,2002);

junkcos = ones(num_months_wind,1) * cos(YY*pi/180);
wind.yymm = 2002.75 + ((1:num_months_wind)-1)/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ocean = find(wind.emissivitychange.landfrac == 0);
stemp = wind.emissivitychange.p0.stemp;
plot(wind.emissivitychange.pnew.efreq(:,ocean),wind.emissivitychange.pnew.emis(:,ocean) - wind.emissivitychange.p0.emis(:,ocean));
plot(wind.emissivitychange.pnew.efreq(:,ocean),nanmean(wind.emissivitychange.pnew.emis(:,ocean) - wind.emissivitychange.p0.emis(:,ocean),2));
plot(wind.emissivitychange.pnew.efreq(:,ocean),nanmean(abs(wind.emissivitychange.pnew.emis(:,ocean) - wind.emissivitychange.p0.emis(:,ocean)),2));

junk900 = find(wind.emissivitychange.pnew.efreq(:,ocean(1)) >= 900,1);
junkrad900_0 = zeros(1,4608);
junkrad900_F = zeros(1,4608);
junkrad900_0(ocean) = wind.emissivitychange.p0.emis(junk900,ocean).*ttorad(900,stemp(ocean));
junkrad900_F(ocean) = wind.emissivitychange.pnew.emis(junk900,ocean).*ttorad(900,stemp(ocean));
scatter_coast(reshape(XX,72,64),reshape(YY,72,64),100,reshape(junkrad900_F - junkrad900_0,72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*0.0005)
  title('Change in radiance at 900 cm-1 \newline due to Wspeed --> Emiss change')

SB = 5.67e-8; %% stefan boltzmann
junkflux900_0 = zeros(1,4608);
junkflux900_F = zeros(1,4608);
junkflux900_0(ocean) = wind.emissivitychange.p0.emis(junk900,ocean).*(stemp(ocean).^4 * SB);
junkflux900_F(ocean) = wind.emissivitychange.pnew.emis(junk900,ocean).*(stemp(ocean).^4 * SB);
scatter_coast(reshape(XX,72,64),reshape(YY,72,64),100,reshape(junkflux900_F - junkflux900_0,72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*0.001)
  title('Change in Planck OLR Flux \newline due to Wspeed --> Emiss change')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jett = jet(128); jett(1,:) = 1;
scatter_coast(XX,YY,100,nanmean(wind.wspeedsave,1)); title([junkstr 'Wind Speed']);        caxis([0 1]*10); colormap(jett)
scatter_coast(XX,YY,100,wind.trend.wspeed);          title([junkstr 'Wind Speed Trends']); caxis([-1 +1]*0.1); colormap(usa2)

disp('computing anomalies for 4608 tiles ... +=1000  .=100')
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+');
  elseif mod(ii,100) == 0
    fprintf(1,'.');
  end

  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.stempsave(:,ii),4);
  wind.stempanom(:,ii) = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.stempsave(:,ii),-1);

  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.wspeedsave(:,ii),4);
  wind.wspeedanom(:,ii) = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.wspeedsave(:,ii),-1);

  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.tccsave(:,ii),4);
  wind.tccanom(:,ii) = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.tccsave(:,ii),-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcolor(reshape(XX,72,64),reshape(YY,72,64),reshape(nanmean(wind.wspeedanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
scatter_coast(reshape(XX,72,64),reshape(YY,72,64),100,reshape(nanmean(wind.wspeedanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1)
title([junkstr 'Mean WSpeed Anomaly \newline over 263 months m/s'])

%moo = wind.wspeedsave(num_months_wind,:) .* cos(YY*pi/180); moo = sum(moo)/sum(cos(YY*pi/180));
wind.wspeedglobalavg = sum(wind.wspeedsave.*junkcos,2)./sum(junkcos,2);
  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.wspeedglobalavg,4);
  wind.wspeedglobalavganom = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.wspeedglobalavg,-1);
plot(wind.yymm,wind.wspeedglobalavg); plotaxis2;                                         title([junkstr 'Global Avg Windspeed m/s']); 
%plot(wind.yymm,wind.wspeedglobalavg - nanmean(wind.wspeedglobalavg)); plotaxis2;        title([junkstr 'Global Avg Windspeed Anomaly m/s']); 
plot(wind.yymm,wind.wspeedglobalavg,wind.yymm,B(1)+wind.wspeedglobalavganom); plotaxis2; title([junkstr 'Global Avg Windspeed Anomaly m/s']); 
plot(wind.yymm,wind.wspeedglobalavg' - (B(1)+wind.wspeedglobalavganom)); plotaxis2;      title([junkstr 'Global Avg Windspeed Signal-Anomaly m/s']); 
plot(wind.yymm,wind.wspeedglobalavganom); plotaxis2;                                     title([junkstr 'Global Avg Windspeed Anomaly m/s']); 

wind.wspeedglobalavg_anom = sum(wind.wspeedanom.*junkcos,2)./sum(junkcos,2);
plot(wind.yymm,wind.wspeedglobalavg_anom); plotaxis2; title([junkstr 'Global Avg Windspeed Anomaly m/s']); 

%%%%%%%%%%%%%%%%%%%%%%%%%

jett = jet(128); jett(1,:) = 1;
scatter_coast(XX,YY,100,nanmean(wind.stempsave,1)); title([junkstr 'Stemp']); caxis([200 320]); colormap(jett)
scatter_coast(XX,YY,100,wind.trend.stemp);          title([junkstr 'STEMP Trends']); caxis([-1 +1]*0.1); colormap(usa2)

%for ii = 1 : 4608
%  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.stempsave(:,ii),4);
%  wind.stempanom(:,ii) = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.stempsave(:,ii),-1);
%end
pcolor(reshape(XX,72,64),reshape(YY,72,64),reshape(nanmean(wind.stempanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*2)
scatter_coast(reshape(XX,72,64),reshape(YY,72,64),100,reshape(nanmean(wind.stempanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*2)
title([junkstr 'Mean Stemp Anomaly \newline over 263 months K'])

%moo = wind.stempsave(num_months_wind,:) .* cos(YY*pi/180); moo = sum(moo)/sum(cos(YY*pi/180));
wind.stempglobalavg      = sum(wind.stempsave.*junkcos,2)./sum(junkcos,2);
  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.stempglobalavg,4);
  wind.stempglobalavganom = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.stempglobalavg,-1);
plot(wind.yymm,wind.stempglobalavg); plotaxis2;                                   title([junkstr 'Global Avg STEMP Anomaly K']); 
%plot(wind.yymm,wind.stempglobalavg - nanmean(wind.stempglobalavg)); plotaxis2;   title([junkstr 'Global Avg STEMP Anomaly K']); 
plot(wind.yymm,wind.stempglobalavg' - (B(1)+wind.stempglobalavganom)); plotaxis2; title([junkstr 'Global Avg STEMP Signal-Anomaly K']); 
plot(wind.yymm,wind.stempglobalavganom); plotaxis2;                               title([junkstr 'Global Avg ERA5 STEMP Anomaly K']); 

%moo = wind.stempanom(num_months_wind,:) .* cos(YY*pi/180); moo = sum(moo)/sum(cos(YY*pi/180));
wind.stempglobalavg_anom = sum(wind.stempanom.*junkcos,2)./sum(junkcos,2);
plot(wind.yymm,wind.stempglobalavg_anom); plotaxis2; title([junkstr 'Global Avg STEMP Anomaly K']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jett = jet(128); jett(1,:) = 1;
scatter_coast(XX,YY,100,nanmean(wind.tccsave,1)); title([junkstr 'Tcc']); caxis([200 320]); colormap(jett)
scatter_coast(XX,YY,100,wind.trend.tcc);          title([junkstr 'TCC Trends']); caxis([-1 +1]*0.01); colormap(usa2)

%for ii = 1 : 4608
%  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.tccsave(:,ii),4);
%  wind.tccanom(:,ii) = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.tccsave(:,ii),-1);
%end
pcolor(reshape(XX,72,64),reshape(YY,72,64),reshape(nanmean(wind.tccanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*2)
scatter_coast(reshape(XX,72,64),reshape(YY,72,64),100,reshape(nanmean(wind.tccanom,1),72,64)); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*0.1)
title([junkstr 'Mean Tcc Anomaly \newline over 263 months K'])

%moo = wind.tccsave(num_months_wind,:) .* cos(YY*pi/180); moo = sum(moo)/sum(cos(YY*pi/180));
wind.tccglobalavg      = sum(wind.tccsave.*junkcos,2)./sum(junkcos,2);
  [B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.tccglobalavg,4);
  wind.tccglobalavganom = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.tccglobalavg,-1);
plot(wind.yymm,wind.tccglobalavg); plotaxis2;                                 title([junkstr 'Global Avg TCC Anomaly K']); 
%plot(wind.yymm,wind.tccglobalavg - nanmean(wind.tccglobalavg)); plotaxis2;   title([junkstr 'Global Avg TCC Anomaly']); 
plot(wind.yymm,wind.tccglobalavg' - (B(1)+wind.tccglobalavganom)); plotaxis2; title([junkstr 'Global Avg TCC Signal-Anomaly']); 
plot(wind.yymm,wind.tccglobalavganom); plotaxis2;                             title([junkstr 'Global Avg ERA5 TCC Anomaly K']); 

%moo = wind.tccanom(num_months_wind,:) .* cos(YY*pi/180); moo = sum(moo)/sum(cos(YY*pi/180));
wind.tccglobalavg_anom = sum(wind.tccanom.*junkcos,2)./sum(junkcos,2);
plot(wind.yymm,wind.tccglobalavg_anom); plotaxis2; title([junkstr 'Global Avg TCC Anomaly']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; scatter_coast(XX,YY,100,wind.trend.wspeed); title([junkstr 'Wind Speed Trends']); caxis([-1 +1]*0.1); colormap(usa2)
figure(2); clf; scatter_coast(XX,YY,100,wind.trend.stemp);  title([junkstr 'STEMP Trends']); caxis([-1 +1]*0.1); colormap(usa2)
figure(3); clf; scatter_coast(XX,YY,100,wind.trend.tcc);    title([junkstr 'TCC Trends']); caxis([-1 +1]*0.01); colormap(usa2)

addpath /home/sergio/MATLABCODE/COLORMAP/COLORBREWER/cbrewer/cbrewer
blues = flipud(cbrewer('seq', 'Blues', 256));

figure(4); clf; plot(wind.yymm,wind.wspeedglobalavg_anom); plotaxis2; title([junkstr 'Global Avg Windspeed Anomaly m/s']); 
figure(5); clf; plot(wind.yymm,wind.stempglobalavg_anom);  plotaxis2; title([junkstr 'Global Avg STEMP Anomaly K']); 
figure(6); clf; plot(wind.yymm,wind.tccglobalavg_anom);    plotaxis2; title([junkstr 'Global Avg TCC Anomaly']); 

figure(7); clf; junk = wind.wspeedsave; junk = squeeze(nanmean(junk,1)); scatter_coast(XX,YY,100,junk); title([junkstr 'Wind Speed Avg']); caxis([0 10]);    colormap(jet)
figure(8); clf; junk = wind.stempsave;  junk = squeeze(nanmean(junk,1)); scatter_coast(XX,YY,100,junk); title([junkstr 'STEMP Avg']);      caxis([200 300]); colormap(jet)
figure(9); clf; junk = wind.tccsave;    junk = squeeze(nanmean(junk,1)); scatter_coast(XX,YY,100,junk); title([junkstr 'TCC Avg']);        caxis([0 1]);     colormap(blues)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modiscloud = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/MODIS_L3_AEROSOL_TRENDS/modis_L3_cloud_trends_2020_0007_2024_0006.mat');

%figure(10); clf; junk = modis.avg_od_liq;  junk = squeeze(nanmean(junk,1)); scatter_coast(XX,YY,100,junk); title('MODIS L3 Wind Speed Avg'); caxis([0 10]);    colormap(jet)
%figure(11); clf; junk = modis.avg_cldtop;  junk = squeeze(nanmean(junk,1)); scatter_coast(XX,YY,100,junk); title('MODIS L3 STEMP Avg');      caxis([200 300]); colormap(jet)
 figure(12); clf; junk = modiscloud.avg_cldfrac(:); scatter_coast(XX,YY,100,junk); title('MODIS L3 TCC Avg');        caxis([0 1]);     colormap(blues)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
[B, stats] = Math_tsfit_lin_robust(junkdoy,wind.stempglobalavg,4);
junkanom = compute_anomaly(1:num_months_wind,junkdoy,B,[],wind.stempglobalavg,-1);

[B, stats] = Math_tsfit_lin_robust(junkdoy-junkdoy(1),wind.stempglobalavg,4);
junkanom = compute_anomaly(1:num_months_wind,junkdoy-junkdoy(1),B,[],wind.stempglobalavg,-1);

[B2, stats2, junkanom2] = compute_anomaly_wrapper(1:length(junkdoy),junkdoy,wind.stempglobalavg,4,[],-1,+1);
plot(wind.yymm,junkanom,'b.',wind.yymm,junkanom2,'r'); plotaxis2;

%see /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/compute_anomaly.m
%[y g] = Math_timeseries_2( dtime(k)-dtime(k(1)), B );   %% y = B(1) + B(2)*dtime + sum_j=1^4 B(j)cos(n*dtime*2*pi) + C(j)sin(n*dtime*2*pi)  = fitted signal IS WRONG
%[y g] = Math_timeseries_2( dtime(k), B );               %% y = B(1) + B(2)*dtime + sum_j=1^4 B(j)cos(n*dtime*2*pi) + C(j)sin(n*dtime*2*pi)  = fitted signal IS CORRECT
%r_anom(k) = (radiance(k) - y') + g(:,2)*B(2);           %% anomaly = (raw signal - fitted signal) + B(2)*dtime
junky = Math_timeseries_2(junkdoy - junkdoy(1), B );
%junky = Math_timeseries_2(junkdoy, B );

sergio = B(1) + B(2)*((junkdoy-junkdoy(1))/365) + B(3)*cos((junkdoy-junkdoy(1))*2*pi/365) + B(4)*sin((junkdoy-junkdoy(1))*2*pi/365); 
sergio2 = zeros(size(sergio));
for ii = 1 : 12
  boo = (1:12:length(sergio2)) + (ii-1);
  boo = boo(boo <= length(sergio));
  sergio2(boo) = wind.stempglobalavg(boo) - mean(wind.stempglobalavg(boo));
end
P = polyfit(wind.yymm,wind.stempglobalavg,1); junkY = polyval(P,wind.yymm); plot(wind.yymm,wind.stempglobalavg,wind.yymm,junkY)

plot(wind.yymm,wind.stempglobalavg,'b.-',wind.yymm,sergio,'g',wind.yymm,junky,'r','linewidth',2); xlim([2020 2025])
plot(wind.yymm,wind.stempglobalavg' - sergio,'b',wind.yymm,sergio2,'cx-',wind.yymm,wind.stempglobalavg' - junky,'g',wind.yymm,wind.stempglobalavganom,'r','linewidth',2); 
  plotaxis2; hl = legend('Obs - sergio_{simple}','Obs - sergio_{simple2}','Obs - math_{timeseries2}','Actual computed anomaly','location','best','fontsize',8); xlim([2020 2025])
%}
