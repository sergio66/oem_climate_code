addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

if ~exist('olr_clr')  
  for ii = 1 : 240
    if mod(ii,100) == 0
      fprintf(1,'+')
    elseif mod(ii,10) == 0
      fprintf(1,'.')
    end
    fname = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    if exist(fname)
      iaFound(ii) = 1;
      boo = load(fname);
      olr_clr(ii,:) = boo.pnew_op.olr_clr;
      olr_cld(ii,:) = boo.pnew_op.olr;
      [yy(ii),mm(ii),dd(ii),hh] = tai2utcSergio(mean(boo.pnew_op.rtime));
      daysSince2002(ii) = change2days(yy(ii),mm(ii),dd(ii),2002);
    else
      iaFound(ii) = 0;
    end
  end
end

fprintf(1,'\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('have read data, doing trends for 4608 tiles')

for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  good = find(iaFound == 1);
  thetime = daysSince2002(good);

  thedata = olr_cld(good,ii);
  [B, stats] = Math_tsfit_lin_robust(thetime,thedata,4);
  trend_olr_cld(ii) = B(2);
  trend_olr_cld_err(ii) = stats.se(2);

  thedata = olr_clr(good,ii);
  [B, stats] = Math_tsfit_lin_robust(thetime,thedata,4);
  trend_olr_clr(ii) = B(2);
  trend_olr_clr_err(ii) = stats.se(2);
end

figure(1); clf; pcolor(reshape(trend_olr_clr,72,64)'); caxis([-1 +1]); title('clrsky'); shading interp;
figure(2); clf; pcolor(reshape(trend_olr_cld,72,64)'); caxis([-1 +1]); title('allsky'); shading interp;
% save era5_monthly_olrtrends_directcomputation_20years.mat trend_olr_clr*  trend_olr_cld* 

boo = reshape(trend_olr_clr,72,64); zonal_olr_clr = nanmean(boo,1);
boo = reshape(trend_olr_cld,72,64); zonal_olr_cld = nanmean(boo,1);

load latB2.txt
rlat = meanvaluebin(latB2);
plot(rlat,zonal_olr_clr,rlat,zonal_olr_cld); plotaxis2; hl = legend('clr','cld');

load ceres_trends_20year.mat
plot(rlat,zonal_olr_clr,rlat,zonal_olr_cld,ceres_trend.trend_lat,ceres_trend.trend_lw_clr,ceres_trend.trend_lat,ceres_trend.trend_lw,'linewidth',2); 
plotaxis2; hl = legend('ERA5 clr','ERA5 cld','CERES clr','CERES cld','location','best','fontsize',10);

disp('now also run compare_OLR_trend.m')
