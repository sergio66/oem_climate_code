disp('see Readme_SKTtrends_tilecenter_tilemean.m')

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

yy = 2002;
mm = 08;

iaNumYears = [20];
for ii = 1 : iaNumYears*12

  mm = mm + 1;
  if mm > 12
    yy = yy + 1;
    mm = 1;
  end
  ysave(ii) = yy;
  msave(ii) = mm;

  doy(ii) = change2days(yy,mm,01,1800);
end

iGISS   = -1;
iL3     = -1;     %% does AIRS v7 and CLIMCAPS
iMERRA2 = -1;
iERA5   = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iGISS > 0
  disp('doing GISS')
  a = load('L3_SKT_TIMESERIES_2002_09_to_2022_08/giss_skt_mean.mat');
  for jj = 1 : 64
    for ii = 1 : 72

      boo = a.cntr_stempIJ_giss(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),1);
        trend.giss.cntr(ii,jj) = B(2);
        trend.giss.cntr_err(ii,jj) = stats.se(2);
      else
        trend.giss.cntr(ii,jj) = nan;
        trend.giss.cntr_err(ii,jj) = nan;
      end

      boo = a.mean_stempIJ_giss(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),1);
        trend.giss.mean(ii,jj) = B(2);
        trend.giss.mean_err(ii,jj) = stats.se(2);
      else
        trend.giss.mean(ii,jj) = nan;
        trend.giss.mean_err(ii,jj) = nan;
      end
    end
  end

  figure(1); clf; scatter_coast(a.mean_rlonIJ_giss,a.mean_rlatIJ_giss,100,trend.giss.cntr); title('GISS CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); clf; scatter_coast(a.mean_rlonIJ_giss,a.mean_rlatIJ_giss,100,trend.giss.mean); title('GISS MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); clf; scatter_coast(a.mean_rlonIJ_giss,a.mean_rlatIJ_giss,100,trend.giss.cntr - trend.giss.mean); title('GISS MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)
  pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear a

if iL3  > 0
  disp('doing L3')
  for tt = 1 : 240
    ab = load(['L3_SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(ysave(tt)) '_' num2str(msave(tt),'%02d') '.mat']);
    a.rlonIJ = ab.rlonIJ;
    a.rlatIJ = ab.rlatIJ;
    a.cntrA_stempIJ_clim(tt,:,:) = ab.cntrA_stempIJ_clim;
    a.meanA_stempIJ_clim(tt,:,:) = ab.meanA_stempIJ_clim;
    a.cntrD_stempIJ_clim(tt,:,:) = ab.cntrD_stempIJ_clim;
    a.meanD_stempIJ_clim(tt,:,:) = ab.meanD_stempIJ_clim;
    a.cntrA_stempIJ_v7(tt,:,:) = ab.cntrA_stempIJ_v7;
    a.meanA_stempIJ_v7(tt,:,:) = ab.meanA_stempIJ_v7;
    a.cntrD_stempIJ_v7(tt,:,:) = ab.cntrD_stempIJ_v7;
    a.meanD_stempIJ_v7(tt,:,:) = ab.meanD_stempIJ_v7;
  end

  disp('nannnnnnnning')
  trend.climcaps.meanA = nan(72,64);
  trend.climcaps.cntrA = nan(72,64);
  trend.climcaps.meanD = nan(72,64);
  trend.climcaps.cntrD = nan(72,64);
  trend.v7.meanA = nan(72,64);
  trend.v7.cntrA = nan(72,64);
  trend.v7.meanD = nan(72,64);
  trend.v7.cntrD = nan(72,64);

  for jj = 1 : 64
    for ii = 1 : 72

      boo = a.cntrA_stempIJ_clim(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.climcaps.cntrA(ii,jj) = B(2);
        trend.climcaps.cntrA_err(ii,jj) = stats.se(2);
      else
        trend.climcaps.cntrA(ii,jj) = nan;
        trend.climcaps.cntrA_err(ii,jj) = nan;
      end

      boo = a.cntrD_stempIJ_clim(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.climcaps.cntrD(ii,jj) = B(2);
        trend.climcaps.cntrD_err(ii,jj) = stats.se(2);
      else
        trend.climcaps.cntrD(ii,jj) = nan;
        trend.climcaps.cntrD_err(ii,jj) = nan;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%

      boo = a.meanA_stempIJ_clim(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.climcaps.meanA(ii,jj) = B(2);
        trend.climcaps.meanA_err(ii,jj) = stats.se(2);
      else
        trend.climcaps.meanA(ii,jj) = nan;
        trend.climcaps.meanA_err(ii,jj) = nan;
      end

      boo = a.meanD_stempIJ_clim(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.climcaps.meanD(ii,jj) = B(2);
        trend.climcaps.meanD_err(ii,jj) = stats.se(2);
      else
        trend.climcaps.meanD(ii,jj) = nan;
        trend.climcaps.meanD_err(ii,jj) = nan;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%

      boo = a.cntrA_stempIJ_v7(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        if B(2) < -10
          plot(doy,boo)
          keyboard_nowindow
        end
        trend.v7.cntrA(ii,jj) = B(2);
        trend.v7.cntrA_err(ii,jj) = stats.se(2);
      else
        trend.v7.cntrA(ii,jj) = nan;
        trend.v7.cntrA_err(ii,jj) = nan;
      end

      boo = a.cntrD_stempIJ_v7(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        if B(2) < -10
          plot(doy,boo)
          keyboard_nowindow
        end
        trend.v7.cntrD(ii,jj) = B(2);
        trend.v7.cntrD_err(ii,jj) = stats.se(2);
      else
        trend.v7.cntrD(ii,jj) = nan;
        trend.v7.cntrD_err(ii,jj) = nan;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%

      boo = a.meanA_stempIJ_v7(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        if B(2) < -10
          plot(doy,boo)
          keyboard_nowindow
        end
        trend.v7.meanA(ii,jj) = B(2);
        trend.v7.meanA_err(ii,jj) = stats.se(2);
      else
        trend.v7.meanA(ii,jj) = nan;
        trend.v7.meanA_err(ii,jj) = nan;
      end

      boo = a.meanD_stempIJ_v7(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        if B(2) < -10
          plot(doy,boo)
          keyboard_nowindow
        end
%        if ii == 43 & jj == 14
%          error('lksg;lsk;lkgsg')
%        end
        trend.v7.meanD(ii,jj) = B(2);
        trend.v7.meanD_err(ii,jj) = stats.se(2);
      else
        trend.v7.meanD(ii,jj) = nan;
        trend.v7.meanD_err(ii,jj) = nan;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%

    end
  end

  figure(1); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.climcaps.cntrA); title('CLIMCAPS CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.climcaps.meanA); title('CLIMCAPS MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.climcaps.cntrA - trend.climcaps.meanA); title('CLIMCAPS MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)
  pause(0.1)

  figure(1); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.v7.cntrA); title('V7 CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.v7.meanA); title('V7 MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.v7.cntrA - trend.v7.meanA); title('V7 MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)
  pause(0.1)

  trend.rlonIJ = a.rlonIJ;
  trend.rlatIJ = a.rlatIJ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear a

if iMERRA2  > 0
  disp('doing MERRA2')
  for tt = 1 : 240
    ab = load(['/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA/SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(ysave(tt)) '_' num2str(msave(tt),'%02d') '.mat']);
    a.rlonIJ = ab.rlonIJ;
    a.rlatIJ = ab.rlatIJ;
    a.cntr_stempIJ(tt,:,:) = ab.cntr_stempIJ;
    a.mean_stempIJ(tt,:,:) = ab.mean_stempIJ;
  end

  for jj = 1 : 64
    for ii = 1 : 72

      boo = a.mean_stempIJ(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.merra2.mean(ii,jj) = B(2);
        trend.merra2.mean_err(ii,jj) = stats.se(2);
      else
        trend.merra2.mean(ii,jj) = nan;
        trend.merra2.mean_err(ii,jj) = nan;
      end

      boo = a.cntr_stempIJ(:,ii,jj);
      good = find(isfinite(boo));
      if length(good) > 20
        [B stats] = Math_tsfit_lin_robust(doy(good),double(boo(good)),4);
        trend.merra2.cntr(ii,jj) = B(2);
        trend.merra2.cntr_err(ii,jj) = stats.se(2);
      else
        trend.merra2.cntr(ii,jj) = nan;
        trend.merra2.cntr_err(ii,jj) = nan;
      end

    end
  end

  figure(1); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.merra2.cntr); title('MERRA2 CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.merra2.mean); title('MERRA2 MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.merra2.cntr - trend.merra2.mean); title('MERRA2 MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)

  trend.rlonIJ = a.rlonIJ;
  trend.rlatIJ = a.rlatIJ;

  pause(0.1)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iERA5 > 0
  clear a
  a = load('~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/skt_timeseries_20years_desc.mat');
  trend.ERA5.meanD     = a.trendERA5.mean;
  trend.ERA5.meanD_err = a.trendERA5.mean_err;;
  trend.ERA5.cntrD     = a.trendERA5.cntr;
  trend.ERA5.cntrD_err = a.trendERA5.cntr_err;;

  clear a
  a = load('~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/skt_timeseries_20years_asc.mat');
  trend.ERA5.meanA     = a.trendERA5.mean;
  trend.ERA5.meanA_err = a.trendERA5.mean_err;;
  trend.ERA5.cntrA     = a.trendERA5.cntr;
  trend.ERA5.cntrA_err = a.trendERA5.cntr_err;;

  trend.rlonIJ = a.rlonIJ;
  trend.rlatIJ = a.rlatIJ;

  figure(1); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.ERA5.cntrA); title('ERA5 CNTR'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.ERA5.meanA); title('ERA5 MEAN'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); clf; scatter_coast(a.rlonIJ,a.rlatIJ,100,trend.ERA5.cntrA - trend.ERA5.meanA); title('ERA5 MEAN-CNTR'); caxis([-1 +1]*0.15/10); colormap(llsmap5)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
comment = 'see driver_SKTtrends_tilecenter_tilemean.m and Readme_SKTtrends_tilecenter_tilemean';
save skt_trends_tilecenter_tilemean.mat trend comment
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('trend')
  load skt_trends_tilecenter_tilemean.mat
end
%% new to test pushing smoothn

clear newz*
%% see compare_SKT_trends_Day_vs_Night.m
newz6 = load('L3_SKT_TIMESERIES_2002_09_to_2022_08/compare_SKT_trends_Day_vs_Night_fig22.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using ORIG from Jan 2025 submission
disp('using ORIG/CNTR data')
newz11 = newz6.newz11;  %% umbc
newz12 = newz6.newz12;  %% v7
newz13 = newz6.newz13;  %% climcaps
newz21 = newz6.newz21;  %% giss
newz22 = newz6.newz22;  %% era5
newz23 = newz6.newz23;  %% merra2
make_table3('ORIG',newz11,newz12,newz13,newz21,newz22,newz23);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% using CNTR
disp('using NEW/CNTR data')
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.cntrA + trend.v7.cntrD)*0.5;
  newz13 = (trend.climcaps.cntrA + trend.climcaps.cntrD)*0.5;
newz21 = trend.giss.cntr;                                   newz22 = (trend.ERA5.cntrD + trend.ERA5.cntrA)*0.5;                             newz23 = trend.merra2.cntr;
make_table3('CNTR',newz11,newz12,newz13,newz21,newz22,newz23);

%% using MEAN
disp('using NEW/MEAN data')
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.meanA + trend.v7.meanD)*0.5;
  newz13 = (trend.climcaps.meanA + trend.climcaps.meanD)*0.5;
newz21 = trend.giss.mean;                                   newz22 = (trend.ERA5.meanD + trend.ERA5.meanA)*0.5;                             newz23 = trend.merra2.mean;
make_table3('MEAN',newz11,newz12,newz13,newz21,newz22,newz23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.meanA + trend.v7.meanD)*0.5;
  newz13 = (trend.climcaps.meanA + trend.climcaps.meanD)*0.5;
newz21 = trend.giss.mean;                                   newz22 = (trend.ERA5.meanD + trend.ERA5.meanA)*0.5;                             newz23 = trend.merra2.mean;

iFig = 22;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = 1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

iFig = 23;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = -1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.cntrA + trend.v7.cntrD)*0.5;
  newz13 = (trend.climcaps.cntrA + trend.climcaps.cntrD)*0.5;
newz21 = trend.giss.cntr;                                   newz22 = (trend.ERA5.cntrD + trend.ERA5.cntrA)*0.5;                             newz23 = trend.merra2.cntr;

iFig = 24;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = 1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);


iFig = 25;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = -1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFig = 20; 
figure(iFig); clf; 

xnewz11 = newz6.newz11 - newz11;
xnewz12 = newz6.newz12 - newz12;
xnewz13 = newz6.newz13 - newz13;
xnewz21 = newz6.newz21 - newz21;
xnewz22 = newz6.newz22 - newz22;
xnewz23 = newz6.newz23 - newz23;
xplotoptions = plotoptions;
xplotoptions.cx = [-1 +1]*0.151/10;
xplotoptions.barstr = 'DELTA OLD-NEW dSKT/dt [K/yr]';
figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(xnewz11,xnewz12,xnewz13,xnewz21,xnewz22,xnewz23,iFig,xplotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_XX_YY_from_X_Y
figure(12); sizefig; ; clf; aslmap(12,rlat65,rlon73,smoothn(reshape(newz6.newz11,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT SMOOTHN');
  caxis([-1 +1]*0.151); colormap(llsmap5);
window = 7;
window = 4;
window = 5;
figure(13); sizefig; ; clf; aslmap(13,rlat65,rlon73,smoothdata(reshape(newz6.newz11,72,64)',"movmean",window),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT SMOOTH MEAN');
  caxis([-1 +1]*0.151); colormap(llsmap5);
figure(14); sizefig; ; clf; aslmap(14,rlat65,rlon73,reshape(newz6.newz11,72,64)',[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT NO SMOOTH');
  caxis([-1 +1]*0.151); colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(22); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_smoothnn');
figure(23); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_nosmooth');
figure(24); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_smoothnn');
figure(25); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_nosmooth');

figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_smoothnn.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_mean_smoothnn');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_nosmooth.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_mean_nosmooth');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_smoothnn.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_cntr_smoothnn');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_nosmooth.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_cntr_nosmooth');
%}
