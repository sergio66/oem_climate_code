disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

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
  data = all.TwSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
  data = all.RHSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
  data = all.mmw(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);
end
fprintf(1,'\n')
warning on
