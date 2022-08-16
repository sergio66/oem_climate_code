disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

%% esee eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/fit_robust_one_lat.m for lag 1
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing surface64 trends .=10')

warning off
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'.')
  end

  data = reshape(all.stemp,228,72,64); data = nanmean(data,2); data = data(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend64_stemp(ii) = B(2);  trend64_stemp_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_stemp_lag(ii) = l(1);
   else
    trend64_stemp_lag(ii) = NaN;
   end
    
  data = reshape(all.TwSurf,228,72,64); data = nanmean(data,2); data = data(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend64_TwSurf(ii) = B(2); trend64_TwSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_TwSurf_lag(ii) = l(1);
   else
    trend64_TwSurf_lag(ii) = NaN;
   end

  data = reshape(all.RHSurf,228,72,64); data = nanmean(data,2); data = data(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend64_RHSurf(ii) = B(2); trend64_RHSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_RHSurf_lag(ii) = l(1);
   else
    trend64_RHSurf_lag(ii) = NaN;
   end

  data = reshape(all.mmw,228,72,64); data = nanmean(data,2); data = data(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend64_mmw(ii) = B(2);    trend64_mmw_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_mmw_lag(ii) = l(1);
   else
    trend64_mmw_lag(ii) = NaN;
   end

end
fprintf(1,'\n')
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('doing surface trends +=1000,x=100,.=10')
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+ \n')
  elseif mod(ii,100) == 0
    fprintf(1,'x')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  data = all.stemp(:,ii);   [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_stemp(ii) = B(2);  trend_stemp_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_stemp_lag(ii) = l(1);
   else
    trend_stemp_lag(ii) = NaN;
   end
    
  data = all.TwSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_TwSurf_lag(ii) = l(1);
   else
    trend_TwSurf_lag(ii) = NaN;
   end

  data = all.RHSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_RHSurf_lag(ii) = l(1);
   else
    trend_RHSurf_lag(ii) = NaN;
   end

  data = all.mmw(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_mmw_lag(ii) = l(1);
   else
    trend_mmw_lag(ii) = NaN;
   end

end
fprintf(1,'\n')
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wah = reshape(trend_stemp,72,64); whos wah;     plot(1:64,nanmean(wah,1),1:64,trend64_stemp); title('Stemp Trend');
wah = reshape(trend_stemp_lag,72,64); whos wah; plot(1:64,nanmean(wah,1),1:64,trend64_stemp_lag); title('Stemp Lag');

wah = reshape(trend_mmw,72,64); whos wah;     plot(1:64,nanmean(wah,1),1:64,trend64_mmw); title('Mmw Trend');
wah = reshape(trend_mmw_lag,72,64); whos wah; plot(1:64,nanmean(wah,1),1:64,trend64_mmw_lag); title('Mmw Lag');


