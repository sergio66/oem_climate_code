error('too long and painful, forget this')

disp(' ')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_atmos_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

if ~exist('iAllorSeasonal')
  iAllorSeasonal = +1;
end

if iAllorSeasonal == 1
  thetimeSeason = 1 : length(dayOFtime);
  iNumCycles = 4;
elseif iAllorSeasonal == -1
  thetimeSeason = find(all.mm == 12 | all.mm == 01 | all.mm == 02);
  iNumCycles = 0;
elseif iAllorSeasonal == -2
  thetimeSeason = find(all.mm == 03 | all.mm == 04 | all.mm == 05);
  iNumCycles = 0;
elseif iAllorSeasonal == -3
  thetimeSeason = find(all.mm == 06 | all.mm == 07 | all.mm == 08);
  iNumCycles = 0;
elseif iAllorSeasonal == -4
  thetimeSeason = find(all.mm == 09 | all.mm == 10 | all.mm == 11);
  iNumCycles = 0;
end

fprintf(1,'doing 1..100 layers OP ptemp,rh,gas_1,gas_3 trends \n');
if isfield(all,'nwp_plevs')
  [mmmm,nnnn,oooo] = size(all.nwp_plevs);
  fprintf(1,' also doing N levs IP raw ERA5/MERRA2/AIRSL3/CMIP ptemp,rh trends, frac, ppmv and gg trends 1 .. %3i \n',nnnn)
end

disp(' ')
disp('<<< commencing atmos trends +=1000,x=100,.=10 >>> ')

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
    %disp('doing N levs IP ptemp,rh trends, frac, ppmv and gg trends')
    [mmmm,nnnn,oooo] = size(all.nwp_plevs);
    for ll = 1 : nnnn
      data = squeeze(all.nwp_ptemp(:,ll,ii));  
      boo = find(isfinite(data));
      boo = intersect(boo,thetimeSeason);
      if length(boo) > 20
        [B,stats,anom_nwp_ptemp(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
      else
        anom_nwp_ptemp(ll,ii,:) = NaN; 
      end

      data = squeeze(all.nwp_rh(:,ll,ii));  
      boo = find(isfinite(data));
      boo = intersect(boo,thetimeSeason);
      if length(boo) > 20
        [B,stats,anom_nwp_rh(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);  
      else
        anom_nwp_rh(ll,ii,:) = NaN; 
      end

      data = squeeze(all.nwp_gas_1(:,ll,ii));
      boo = find(isfinite(data));
      boo = intersect(boo,thetimeSeason);
      if length(boo) > 20
        [B,stats,anom_nwp_gg(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
      else
        anom_nwp_gg(ll,ii,:) = NaN; 
      end

      data  = squeeze(all.nwp_gas_1(:,ll,ii));
      dataP = squeeze(all.nwp_plevs(:,ll,ii));
      dataT = squeeze(all.nwp_ptemp(:,ll,ii));
      data = toppmv(dataP,dataT,data,18,21);
      boo = find(isfinite(data));
      boo = intersect(boo,thetimeSeason);
      if length(boo) > 20
        [B,stats,anom_nwp_ppmv(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
      else
        anom_nwp_ppmv(ll,ii,:) = NaN; 
      end

      data = squeeze(all.nwp_gas_1(:,ll,ii));
      data = data/nanmean(data);
      boo = find(isfinite(data));
      boo = intersect(boo,thetimeSeason);
      if length(boo) > 20
        [B,stats,anom_nwp_frac(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
      else
        anom_nwp_frac(ll,ii,:) = NaN; 
      end

    end
  end

  %fprintf(1,'doing 100 layers OP ptemp,rh,gas_1,gas_3 trends ii = %4i of 4608 \n',ii)
  for ll = 1 : 100    
    data = squeeze(all.ptemp(:,ll,ii));  
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B,stats,anom_ptemp(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
    else
      anom_ptemp(ll,ii,:) = NaN; 
    end

    data = squeeze(all.gas_1(:,ll,ii));  data = data/mean(data); 
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B,stats,anom_gas_1(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
    else
      anom_gas_1(ll,ii,:) = NaN; 
    end

    data = squeeze(all.gas_3(:,ll,ii));  data = data/mean(data); 
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B,stats,anom_gas_3(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
    else
      anom_gas_3(ll,ii,:) = NaN; 
    end

    data = squeeze(all.RH(:,ll,ii));                             
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B,stats,anom_RH(ll,ii,:)] = compute_anomaly_wrapper(boo,dayOFtime,data,iNumCycles,-1,-1);
    else
      anom_RH(ll,ii,:) = NaN; 
    end
  end
end
fprintf(1,'\n')
warning on
