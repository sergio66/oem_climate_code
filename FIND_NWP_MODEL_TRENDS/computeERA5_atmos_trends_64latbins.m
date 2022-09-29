disp(' ')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

disp('doing atmos trends .=10')
warning off
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'.')
  end

  if isfield(all,'nwp_plevs')
    disp('doing N levs IP ptemp,rh trends, frac, ppmv and gg trends')
    [mmmm,nnnn,oooo] = size(all.nwp_plevs);
    for ll = 1 : nnnn
      data = squeeze(all.nwp_ptemp(:,ll,:));  
        data = reshape(data,mmmm,72,64); 
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend64_nwp_ptemp(ll,ii) = B(2);  trend64_nwp_ptemp_err(ll,ii) = stats.se(2);
      else
        trend64_nwp_ptemp(ll,ii) = NaN; 
        trend64_nwp_ptemp_err(ll,ii) = NaN;
      end
      k = remove_nan(data);
      if length(k) > 50
        l = xcorr(data(k),1,'coeff');
        trend64_nwp_ptemp_lag(ll,ii) = l(1);
      else
        trend64_nwp_ptemp_lag(ll,ii) = NaN;
      end

      data = squeeze(all.nwp_rh(:,ll,:));  
        data = reshape(data,mmmm,72,64);
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend64_nwp_rh(ll,ii) = B(2);  trend64_nwp_rh_err(ll,ii) = stats.se(2);
      else
        trend64_nwp_rh(ll,ii) = NaN; 
        trend64_nwp_rh_err(ll,ii) = NaN;
      end
      k = remove_nan(data);
      if length(k) > 50
        l = xcorr(data(k),1,'coeff');
        trend64_nwp_rh_lag(ll,ii) = l(1);
      else
        trend64_nwp_rh_lag(ll,ii) = NaN;
      end

      data = squeeze(all.nwp_gas_1(:,ll,:));
        data = reshape(data,mmmm,72,64);
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend64_nwp_gg(ll,ii) = B(2);  trend64_nwp_gg_err(ll,ii) = stats.se(2);
      else
        trend64_nwp_gg(ll,ii) = NaN; 
        trend64_nwp_gg_err(ll,ii) = NaN;
      end
      k = remove_nan(data);
      if length(k) > 50
        l = xcorr(data(k),1,'coeff');
        trend64_nwp_gg_lag(ll,ii) = l(1);
      else
        trend64_nwp_gg_lag(ll,ii) = NaN;
      end

      data  = squeeze(all.nwp_gas_1(:,ll,:));
        data = reshape(data,mmmm,72,64);
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      dataP = squeeze(all.nwp_plevs(:,ll,:));
        dataP = reshape(dataP,mmmm,72,64);
        dataP = squeeze(nanmean(dataP,2));
        dataP = dataP(:,ii);
      dataT = squeeze(all.nwp_ptemp(:,ll,:));
        dataT = reshape(dataT,mmmm,72,64);
        dataT = squeeze(nanmean(dataT,2));
        dataT = dataT(:,ii);
      data = toppmv(dataP,dataT,data,18,21);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend64_nwp_ppmv(ll,ii) = B(2);  trend64_nwp_ppmv_err(ll,ii) = stats.se(2);
      else
        trend64_nwp_ppmv(ll,ii) = NaN; 
        trend64_nwp_ppmv_err(ll,ii) = NaN;
      end
      k = remove_nan(data);
      if length(k) > 50
        l = xcorr(data(k),1,'coeff');
        trend64_nwp_ppmv_lag(ll,ii) = l(1);
      else
        trend64_nwp_ppmv_lag(ll,ii) = NaN;
      end

      data = squeeze(all.nwp_gas_1(:,ll,:));
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
      data = data/nanmean(data);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); 
        trend64_nwp_frac(ll,ii) = B(2);  trend64_nwp_frac_err(ll,ii) = stats.se(2);
      else
        trend64_nwp_frac(ll,ii) = NaN; 
        trend64_nwp_frac_err(ll,ii) = NaN;
      end
      k = remove_nan(data);
      if length(k) > 50
        l = xcorr(data(k),1,'coeff');
        trend64_nwp_gas_1_lag(ll,ii) = l(1);
      else
        trend64_nwp_gas_1_lag(ll,ii) = NaN;
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ll = 1 : 100    
    disp('doing 100 layers OP ptemp,rh,gas_1,gas_3 trends')
    data = squeeze(all.ptemp(:,ll,:));  
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    boo = find(isfinite(data)  & iaIndexUse' == 1);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend64_ptemp(ll,ii) = B(2);  trend64_ptemp_err(ll,ii) = stats.se(2);
    else
      trend64_ptemp(ll,ii) = NaN; 
      trend64_ptemp_err(ll,ii) = NaN;
    end
    k = remove_nan(data);
    if length(k) > 50
      l = xcorr(data(k),1,'coeff');
      trend64_ptemp_lag(ll,ii) = l(1);
    else
      trend64_ptemp_lag(ll,ii) = NaN;
    end

    data = squeeze(all.gas_1(:,ll,:));  
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    data = data/mean(data); 
    boo = find(isfinite(data)  & iaIndexUse' == 1);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend64_gas_1(ll,ii) = B(2);  trend64_gas_1_err(ll,ii) = stats.se(2);
    else
      trend64_gas_1(ll,ii) = NaN; 
      trend64_gas_1_err(ll,ii) = NaN;
    end
    k = remove_nan(data);
    if length(k) > 50
      l = xcorr(data(k),1,'coeff');
      trend64_gas_1_lag(ll,ii) = l(1);
    else
      trend64_gas_1_lag(ll,ii) = NaN;
    end

    data = squeeze(all.gas_3(:,ll,:));  
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    data = data/mean(data); 
    boo = find(isfinite(data) & iaIndexUse' == 1);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend64_gas_3(ll,ii) = B(2);  trend64_gas_3_err(ll,ii) = stats.se(2);
    else
      trend64_gas_3(ll,ii) = NaN; 
      trend64_gas_3_err(ll,ii) = NaN;
    end
    k = remove_nan(data);
    if length(k) > 50
      l = xcorr(data(k),1,'coeff');
      trend64_gas_3_lag(ll,ii) = l(1);
    else
      trend64_gas_3_lag(ll,ii) = NaN;
    end

    data = squeeze(all.RH(:,ll,:));                             
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    boo = find(isfinite(data) & iaIndexUse' == 1);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4); trend64_RH(ll,ii) = B(2);     trend64_RH_err(ll,ii) = stats.se(2);
    else
      trend64_RH(ll,ii) = NaN; 
      trend64_RH_err(ll,ii) = NaN;
    end
    k = remove_nan(data);
    if length(k) > 50
      l = xcorr(data(k),1,'coeff');
      trend64_RH_lag(ll,ii) = l(1);
    else
      trend64_RH_lag(ll,ii) = NaN;
    end
  end
end
fprintf(1,'\n')
warning on
