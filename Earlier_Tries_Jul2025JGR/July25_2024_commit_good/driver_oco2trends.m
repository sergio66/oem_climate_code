addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

iCnt = 0;
for yy = 2015:2021
  for mm = 1 : 12
    iCnt = iCnt + 1;
    dirname = '//asl/s1/sergio/OCO2_L3/';
    fname = ['oco2_GEOS_L3CO2_month_' num2str(yy) num2str(mm,'%02d') '_B10206Ar.nc4'];
    daysince2002(iCnt) = change2days(yy,mm,15,2002);
    if exist([dirname '/' fname])
      iaFound(iCnt) = 1;
    else
      iaFound(iCnt) = 0;
    end
  end
end

sum(iaFound)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iCnt = 0;
for yy = 2015:2021
  for mm = 1 : 12
    iCnt = iCnt + 1;
    dirname = '//asl/s1/sergio/OCO2_L3/';
    fname = ['oco2_GEOS_L3CO2_month_' num2str(yy) num2str(mm,'%02d') '_B10206Ar.nc4'];

    a = read_netcdf_lls([dirname '/' fname]);

    clear co2_64x72
    for jj = 1 : length(rlat65)-1
      for ii = 1 : length(rlon73)-1
        booY = find(a.lat >= rlat65(jj) & a.lat < rlat65(jj+1));
        booX = find(a.lon >= rlon73(ii) & a.lon < rlon73(ii+1));
        XCO2 = a.XCO2(booX,booY);
        co2_64x72(ii,jj) = nanmean(XCO2(:));
      end
    end
    timeseries_co2_64x72(iCnt,:,:) = co2_64x72;
    
    figure(1); pcolor(a.XCO2'); shading interp; colorbar; colormap jet; title([num2str(yy) '/' num2str(mm,'%02d')])
    figure(2); pcolor(co2_64x72'); shading interp; colorbar; colormap jet; title([num2str(yy) '/' num2str(mm,'%02d')])
    pause(0.1);
  end
end

save oco2_timeseries.mat timeseries_co2_64x72 rlat65 rlon73
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 72
  fprintf(1,'ii = %2i of 72 \n',ii);
  for jj = 1 : 64
    data = squeeze(timeseries_co2_64x72(:,ii,jj));
    [B, stats] = Math_tsfit_lin_robust(daysince2002,data*1e6,4);
    co2_trend(ii,jj)     = B(2);
    co2_trend_err(ii,jj) = stats.se(2);
  end
end
save oco2_timeseries.mat timeseries_co2_64x72 rlat65 rlon73 co2_trend co2_trend_err
figure(1); clf; pcolor(co2_trend'); title('OCO2 trends 2015/01 to 2021/12'); shading interp; colorbar
figure(2); clf; pcolor(co2_trend_err');

addpath /umbc/xfs2/strow/asl/matlib/maps
aslmap(1,rlat65,rlon73,smoothn((reshape(co2_trend,72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  
title('OCO2 trends 2015/01 to 2021/12 d/dt CO2');
