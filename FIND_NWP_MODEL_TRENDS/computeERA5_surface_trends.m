disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

if ~exist('iAllorSeasonal')
  iAllorSeasonal = +1;
end

fprintf(1,'computeERA5_surface_trends.m : iAllorSeasonal = %2i \n',iAllorSeasonal)

disp('doing surface trends +=1000,x=100,.=10')
warning off
if iAllorSeasonal == +1
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
else
  if iAllorSeasonal == -1
    thetimeSeason = find(all.mm == 12 | all.mm == 01 | all.mm == 02);
  elseif iAllorSeasonal == -2
    thetimeSeason = find(all.mm == 03 | all.mm == 04 | all.mm == 05);
  elseif iAllorSeasonal == -3
    thetimeSeason = find(all.mm == 06 | all.mm == 07 | all.mm == 08);
  elseif iAllorSeasonal == -4
    thetimeSeason = find(all.mm == 09 | all.mm == 10 | all.mm == 11);
  end

  for ii = 1 : 4608
    if mod(ii,1000) == 0
      fprintf(1,'+ \n')
    elseif mod(ii,100) == 0
      fprintf(1,'x')
    elseif mod(ii,10) == 0
      fprintf(1,'.')
    end
    data = all.stemp(thetimeSeason,ii);   [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_stemp(ii) = B(2);  trend_stemp_err(ii) = stats.se(2);
    data = all.TwSurf(thetimeSeason,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
    data = all.RHSurf(thetimeSeason,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
    data = all.mmw(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);
  end

end

fprintf(1,'\n')
warning on
