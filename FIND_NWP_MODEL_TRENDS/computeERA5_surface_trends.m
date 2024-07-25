disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

if ~exist('iAllorSeasonal')
  iAllorSeasonal = +1;
end

fprintf(1,'computeERA5_surface_trends.m : iAllorSeasonal = %2i \n',iAllorSeasonal)

if iOLR > 0
  disp('   .... iOLR > 0 so doing d2m/t2m/olr/ilr etc trends ...')
end

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

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = all.d2m(:,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_d2m(ii) = B(2);      trend_d2m_err(ii) = stats.se(2);
        data = all.t2m(:,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_t2m(ii) = B(2);      trend_t2m_err(ii) = stats.se(2);
        data = all.RH2m(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RH2m(ii) = B(2);     trend_RH2m_err(ii) = stats.se(2);
        data = all.e2a(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_e2m(ii) = B(2);      trend_e2m_err(ii) = stats.se(2);
        data = all.e2a(:,ii)/nanmean(all.e2a(:,ii));     
                                  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_frac_e2m(ii) = B(2); trend_frac_e2m_err(ii) = stats.se(2);
        data = all.ecs(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ecs(ii) = B(2);      trend_ecs_err(ii) = stats.se(2);
        data = all.Rld(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_Rld(ii) = B(2);  trend_ilr_Rld_err(ii) = stats.se(2);
      end        
      data = all.olr(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr(ii) = B(2);      trend_olr_err(ii) = stats.se(2);
      data = all.olr_clr(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr_clr(ii) = B(2);  trend_olr_clr_err(ii) = stats.se(2);
      data = all.ilr(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr(ii) = B(2);      trend_ilr_err(ii) = stats.se(2);
      data = all.ilr_clr(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_clr(ii) = B(2);  trend_ilr_clr_err(ii) = stats.se(2);
      data = all.ilr_adj(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_adj(ii) = B(2);  trend_ilr_adj_err(ii) = stats.se(2);
    end

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

    if iOLR > 0
      data = all.d2m(thetimeSeason,ii);       [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_d2m(ii) = B(2);      trend_d2m_err(ii) = stats.se(2);
      data = all.t2m(thetimeSeason,ii);       [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_t2m(ii) = B(2);      trend_t2m_err(ii) = stats.se(2);
      data = all.RH2m(thetimeSeason,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RH2m(ii) = B(2);     trend_RH2m_err(ii) = stats.se(2);

      data = all.olr(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr(ii) = B(2);      trend_olr_err(ii) = stats.se(2);
      data = all.olr_clr(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr_clr(ii) = B(2);  trend_olr_clr_err(ii) = stats.se(2);
      data = all.ilr(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr(ii) = B(2);      trend_ilr_err(ii) = stats.se(2);
      data = all.ilr_clr(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_clr(ii) = B(2);  trend_ilr_clr_err(ii) = stats.se(2);
      data = all.ilr_adj(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_adj(ii) = B(2);  trend_ilr_adj_err(ii) = stats.se(2);
      data = all.e2a(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_e2m(ii) = B(2);      trend_e2m_err(ii) = stats.se(2);
      data = all.e2a(thetimeSeason,ii)/nanmean(all.e2a(thetimeSeason,ii));     
                                [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_frac_e2m(ii) = B(2); trend_frac_e2m_err(ii) = stats.se(2);
      data = all.ecs(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ecs(ii) = B(2);      trend_ecs_err(ii) = stats.se(2);
      data = all.Rld(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_Rld(ii) = B(2);  trend_ilr_Rld_err(ii) = stats.se(2);
    end

  end

end

fprintf(1,'\n')
warning on
