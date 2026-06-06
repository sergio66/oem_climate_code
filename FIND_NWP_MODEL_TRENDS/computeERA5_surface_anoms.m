disp(' ')
disp('computeERA5_surface_anoms.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_surface_anoms.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_surface_anoms.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp(' ')

if ~exist('iAllorSeasonal')
  iAllorSeasonal = +1;
end

fprintf(1,'computeERA5_surface_anoms.m : iAllorSeasonal = %2i \n',iAllorSeasonal)

if iOLR > 0
  disp('   .... iOLR > 0 so doing d2m/t2m/olr/ilr etc anoms ...')
end

disp('doing surface anoms +=1000,x=100,.=10')
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
    data = pall.stemp(:,ii);   k = find(isfinite(data)); [B,err,anom_stemp(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.TwSurf(:,ii);  k = find(isfinite(data)); [B,err,anom_TwSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.RHSurf(:,ii);  k = find(isfinite(data)); [B,err,anom_RHSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.mmw(:,ii);     k = find(isfinite(data)); [B,err,anom_mmw(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = pall.d2m(:,ii);     k = find(isfinite(data)); [B,err,anom_d2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.t2m(:,ii);     k = find(isfinite(data)); [B,err,anom_t2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.RH2m(:,ii);    k = find(isfinite(data)); [B,err,anom_RH2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.e2a(:,ii);     k = find(isfinite(data)); [B,err,anom_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.e2a(:,ii)/nanmean(pall.e2a(:,ii));     
                                  k = find(isfinite(data)); [B,err,anom_frac_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = pall.ecs(:,ii);     k = find(isfinite(data)); [B,err,anom_ecs(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = pall.Rld(:,ii);     k = find(isfinite(data)); [B,err,anom_Rld(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

      end        
      data = pall.olr(:,ii);       k = find(isfinite(data)); [B,err,anom_olr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.olr_clr(:,ii);   k = find(isfinite(data)); [B,err,anom_olr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr(:,ii);       k = find(isfinite(data)); [B,err,anom_ilr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr_clr(:,ii);   k = find(isfinite(data)); [B,err,anom_ilr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr_adj(:,ii);   k = find(isfinite(data)); [B,err,anom_ilr_adj(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
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

    data = pall.stemp(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_stemp(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.TwSurf(thetimeSeason,ii);  k = find(isfinite(data)); [B,err,anom_TwSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.RHSurf(thetimeSeason,ii);  k = find(isfinite(data)); [B,err,anom_RHSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = pall.mmw(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_mmw(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = pall.d2m(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_d2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.t2m(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_t2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.RH2m(thetimeSeason,ii);    k = find(isfinite(data)); [B,err,anom_RH2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.e2a(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = pall.e2a(thetimeSeason,ii)/nanmean(pall.e2a(thetimeSeason,ii));     
                                  k = find(isfinite(data)); [B,err,anom_frac_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = pall.ecs(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_ecs(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = pall.Rld(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_Rld(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

      end        
      data = pall.olr(thetimeSeason,ii);       k = find(isfinite(data)); [B,err,anom_olr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.olr_clr(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_olr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr(thetimeSeason,ii);       k = find(isfinite(data)); [B,err,anom_ilr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr_clr(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_ilr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = pall.ilr_adj(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_ilr_adj(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    end
  end

end

fprintf(1,'\n')
warning on
