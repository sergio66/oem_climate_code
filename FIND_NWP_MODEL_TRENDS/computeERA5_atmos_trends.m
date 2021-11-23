disp(' ')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

disp('doing atmos trends +=1000,x=100,.=10')
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+ \n')
  elseif mod(ii,100) == 0
    fprintf(1,'x')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  if isfield(all,'nwp_plevs')
    disp('doing N levs IP ptemp,rh trends')
    [mmmm,nnnn,oooo] = size(all.nwp_plevs);
    for ll = 1 : nnnn
      data = squeeze(all.nwp_ptemp(:,ll,ii));  
      boo = find(isfinite(data));
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend_nwp_ptemp(ll,ii) = B(2);  trend_nwp_ptemp_err(ll,ii) = stats.se(2);
      else
        trend_nwp_ptemp(ll,ii) = NaN; 
        trend_nwp_ptemp_err(ll,ii) = NaN;
      end

      data = squeeze(all.nwp_rh(:,ll,ii));  
      boo = find(isfinite(data));
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend_nwp_rh(ll,ii) = B(2);  trend_nwp_rh_err(ll,ii) = stats.se(2);
      else
        trend_nwp_rh(ll,ii) = NaN; 
        trend_nwp_rh_err(ll,ii) = NaN;
      end
    end
  end

  for ll = 1 : 100    
    disp('doing 100 layers OP ptemp,rh,gas_1,gas_3 trends')
    data = squeeze(all.ptemp(:,ll,ii));  
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend_ptemp(ll,ii) = B(2);  trend_ptemp_err(ll,ii) = stats.se(2);
    else
      trend_ptemp(ll,ii) = NaN; 
      trend_ptemp_err(ll,ii) = NaN;
    end

    data = squeeze(all.gas_1(:,ll,ii));  data = data/mean(data); 
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend_gas_1(ll,ii) = B(2);  trend_gas_1_err(ll,ii) = stats.se(2);
    else
      trend_gas_1(ll,ii) = NaN; 
      trend_gas_1_err(ll,ii) = NaN;
    end

    data = squeeze(all.gas_3(:,ll,ii));  data = data/mean(data); 
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend_gas_3(ll,ii) = B(2);  trend_gas_3_err(ll,ii) = stats.se(2);
    else
      trend_gas_3(ll,ii) = NaN; 
      trend_gas_3_err(ll,ii) = NaN;
    end

    data = squeeze(all.RH(:,ll,ii));                             
    boo = find(isfinite(data));
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend_RH(ll,ii) = B(2);     trend_RH_err(ll,ii) = stats.se(2);
    else
      trend_RH(ll,ii) = NaN; 
      trend_RH_err(ll,ii) = NaN;
    end
  end
end
fprintf(1,'\n')
warning on
