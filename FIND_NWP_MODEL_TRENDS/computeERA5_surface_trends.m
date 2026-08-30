disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

if ~exist('iAllorSeasonal')
  iAllorSeasonal = +1;
end

fprintf(1,'computeERA5_surface_trends.m : iAllorSeasonal = %2i \n',iAllorSeasonal)

if ~exist('iOLR')
  iOLR = -1;
end  
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
    data = pall.stemp(:,ii);   [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_stemp(ii) = B(2);  trend_stemp_err(ii) = stats.se(2);
    data = pall.TwSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
    data = pall.RHSurf(:,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
    data = pall.mmw(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = pall.d2m(:,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_d2m(ii) = B(2);      trend_d2m_err(ii) = stats.se(2);
        data = pall.t2m(:,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_t2m(ii) = B(2);      trend_t2m_err(ii) = stats.se(2);
        data = pall.RH2m(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RH2m(ii) = B(2);     trend_RH2m_err(ii) = stats.se(2);
        data = pall.e2a(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_e2m(ii) = B(2);      trend_e2m_err(ii) = stats.se(2);
        data = pall.e2a(:,ii)/nanmean(pall.e2a(:,ii));     
                                  [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_frac_e2m(ii) = B(2); trend_frac_e2m_err(ii) = stats.se(2);
        data = pall.ecs(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ecs(ii) = B(2);      trend_ecs_err(ii) = stats.se(2);
        data = pall.Rld(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_Rld(ii) = B(2);  trend_ilr_Rld_err(ii) = stats.se(2);
      end        
      data = pall.olr(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr(ii) = B(2);      trend_olr_err(ii) = stats.se(2);
      data = pall.olr_clr(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr_clr(ii) = B(2);  trend_olr_clr_err(ii) = stats.se(2);
      data = pall.ilr(:,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr(ii) = B(2);      trend_ilr_err(ii) = stats.se(2);
      data = pall.ilr_clr(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_clr(ii) = B(2);  trend_ilr_clr_err(ii) = stats.se(2);
      data = pall.ilr_adj(:,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_adj(ii) = B(2);  trend_ilr_adj_err(ii) = stats.se(2);
    end

  end

else
  if iAllorSeasonal == -1
    thetimeSeason = find(pall.mm == 12 | pall.mm == 01 | pall.mm == 02);
  elseif iAllorSeasonal == -2
    thetimeSeason = find(pall.mm == 03 | pall.mm == 04 | pall.mm == 05);
  elseif iAllorSeasonal == -3
    thetimeSeason = find(pall.mm == 06 | pall.mm == 07 | pall.mm == 08);
  elseif iAllorSeasonal == -4
    thetimeSeason = find(pall.mm == 09 | pall.mm == 10 | pall.mm == 11);
  end

  for ii = 1 : 4608
    if mod(ii,1000) == 0
      fprintf(1,'+ \n')
    elseif mod(ii,100) == 0
      fprintf(1,'x')
    elseif mod(ii,10) == 0
      fprintf(1,'.')
    end
    data = pall.stemp(thetimeSeason,ii);   [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_stemp(ii) = B(2);  trend_stemp_err(ii) = stats.se(2);
    data = pall.TwSurf(thetimeSeason,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
    data = pall.RHSurf(thetimeSeason,ii);  [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
    data = pall.mmw(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime(thetimeSeason),data,0); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);

    if iOLR > 0
      data = pall.d2m(thetimeSeason,ii);       [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_d2m(ii) = B(2);      trend_d2m_err(ii) = stats.se(2);
      data = pall.t2m(thetimeSeason,ii);       [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_t2m(ii) = B(2);      trend_t2m_err(ii) = stats.se(2);
      data = pall.RH2m(thetimeSeason,ii);      [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_RH2m(ii) = B(2);     trend_RH2m_err(ii) = stats.se(2);

      data = pall.olr(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr(ii) = B(2);      trend_olr_err(ii) = stats.se(2);
      data = pall.olr_clr(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_olr_clr(ii) = B(2);  trend_olr_clr_err(ii) = stats.se(2);
      data = pall.ilr(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr(ii) = B(2);      trend_ilr_err(ii) = stats.se(2);
      data = pall.ilr_clr(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_clr(ii) = B(2);  trend_ilr_clr_err(ii) = stats.se(2);
      data = pall.ilr_adj(thetimeSeason,ii); [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_adj(ii) = B(2);  trend_ilr_adj_err(ii) = stats.se(2);
      data = pall.e2a(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_e2m(ii) = B(2);      trend_e2m_err(ii) = stats.se(2);
      data = pall.e2a(thetimeSeason,ii)/nanmean(pall.e2a(thetimeSeason,ii));     
                                [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_frac_e2m(ii) = B(2); trend_frac_e2m_err(ii) = stats.se(2);
      data = pall.ecs(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ecs(ii) = B(2);      trend_ecs_err(ii) = stats.se(2);
      data = pall.Rld(thetimeSeason,ii);     [B, stats] = Math_tsfit_lin_robust(dayOFtime,data,4); trend_ilr_Rld(ii) = B(2);  trend_ilr_Rld_err(ii) = stats.se(2);
    end

  end

end

fprintf(1,'\n')
warning on
