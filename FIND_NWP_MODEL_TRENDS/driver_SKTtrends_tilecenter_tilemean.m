disp('see Readme_SKTtrends_tilecenter_tilemean.m')

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/NANROUTINES
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
else
  disp('skipping processing SKT time series from GISS, will read in saved mat file')
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

else
  disp('skipping processing SKT time series from AIRSv7 and CLIMCAPS, will read in saved mat file')
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

else
  disp('skipping processing SKT time series from MERRA2, will read in saved mat file')
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

else
  disp('skipping processing SKT time series from ERA5, will read in saved mat file')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
comment = 'see driver_SKTtrends_tilecenter_tilemean.m and Readme_SKTtrends_tilecenter_tilemean';
save skt_trends_tilecenter_tilemean.mat trend comment
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this makes figs 181,182,183 which are fig8 of the accepted paper
process_SKTtrends_tilecenter_tilemean

