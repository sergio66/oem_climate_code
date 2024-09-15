disp(' ')
disp('computeERA5_atmos_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_atmos_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_atmos_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp(' ')

% usetime = 1 : length(dayOFtime);

disp('doing atmos trends and anomalies .=10')
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
        k = find(isfinite(data(usetime))); 
        [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_nwp_ptemp(ll,ii,:) = junkanom;
      else
        anom64_nwp_ptemp(ll,ii,:) = NaN;
      end

      data = squeeze(all.nwp_rh(:,ll,:));  
        data = reshape(data,mmmm,72,64);
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        k = find(isfinite(data(usetime))); 
        [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_nwp_rh(ll,ii,:) = junkanom;
      else
        anom64_nwp_rh(ll,ii,:) = NaN; 
      end

      data = squeeze(all.nwp_gas_1(:,ll,:));
        data = reshape(data,mmmm,72,64);
        data = squeeze(nanmean(data.*maskLF,2));
        data = data(:,ii);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        k = find(isfinite(data(usetime))); 
        [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_nwp_gg(ll,ii,:) = junkanom;
      else
        anom64_nwp_gg(ll,ii,:) = NaN; 
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
        k = find(isfinite(data(usetime))); 
        [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_nwp_ppmv(ll,ii,:) = junkanom;
      else
        anom64_nwp_ppmv(ll,ii,:) = NaN; 
      end

      data = squeeze(all.nwp_gas_1(:,ll,:));
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
      data = data/nanmean(data);
      boo = find(isfinite(data) & iaIndexUse' == 1);
      if length(boo) > 20
        k = find(isfinite(data(usetime))); 
        [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_nwp_frac(ll,ii,:) = junkanom;
      else
        anom64_nwp_frac(ll,ii,:) = NaN; 
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
      k = find(isfinite(data(usetime))); 
      [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_ptemp(ll,ii,:) = junkanom;
    else
      anom64_ptemp(ll,ii,:) = NaN; 
    end

    data = squeeze(all.gas_1(:,ll,:));  
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    data = data/mean(data); 
    boo = find(isfinite(data)  & iaIndexUse' == 1);
    if length(boo) > 20
      k = find(isfinite(data(usetime))); 
      [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_gas_1(ll,ii,:) = junkanom;
    else
      anom64_gas_1(ll,ii,:) = NaN; 
    end

    data = squeeze(all.gas_3(:,ll,:));  
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    data = data/mean(data); 
    boo = find(isfinite(data) & iaIndexUse' == 1);
    if length(boo) > 20
      k = find(isfinite(data(usetime))); 
      [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_gas_3(ll,ii,:) = junkanom;
    else
      anom64_gas_3(ll,ii,:) = NaN; 
    end

    data = squeeze(all.RH(:,ll,:));                             
      data = reshape(data,mmmm,72,64);
      data = squeeze(nanmean(data.*maskLF,2));
      data = data(:,ii);
    boo = find(isfinite(data) & iaIndexUse' == 1);
    if length(boo) > 20
      k = find(isfinite(data(usetime))); 
      [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_RH(ll,ii,:) = junkanom;
    else
      anom64_RH(ll,ii,:) = NaN; 
    end
  end
end
fprintf(1,'\n')
warning on
