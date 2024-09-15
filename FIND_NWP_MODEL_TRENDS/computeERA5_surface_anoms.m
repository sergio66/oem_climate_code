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
    data = all.stemp(:,ii);   k = find(isfinite(data)); [B,err,anom_stemp(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.TwSurf(:,ii);  k = find(isfinite(data)); [B,err,anom_TwSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.RHSurf(:,ii);  k = find(isfinite(data)); [B,err,anom_RHSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.mmw(:,ii);     k = find(isfinite(data)); [B,err,anom_mmw(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = all.d2m(:,ii);     k = find(isfinite(data)); [B,err,anom_d2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.t2m(:,ii);     k = find(isfinite(data)); [B,err,anom_t2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.RH2m(:,ii);    k = find(isfinite(data)); [B,err,anom_RH2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.e2a(:,ii);     k = find(isfinite(data)); [B,err,anom_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.e2a(:,ii)/nanmean(all.e2a(:,ii));     
                                  k = find(isfinite(data)); [B,err,anom_frac_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = all.ecs(:,ii);     k = find(isfinite(data)); [B,err,anom_ecs(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = all.Rld(:,ii);     k = find(isfinite(data)); [B,err,anom_Rld(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

      end        
      data = all.olr(:,ii);       k = find(isfinite(data)); [B,err,anom_olr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.olr_clr(:,ii);   k = find(isfinite(data)); [B,err,anom_olr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr(:,ii);       k = find(isfinite(data)); [B,err,anom_ilr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr_clr(:,ii);   k = find(isfinite(data)); [B,err,anom_ilr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr_adj(:,ii);   k = find(isfinite(data)); [B,err,anom_ilr_adj(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
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

    data = all.stemp(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_stemp(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.TwSurf(thetimeSeason,ii);  k = find(isfinite(data)); [B,err,anom_TwSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.RHSurf(thetimeSeason,ii);  k = find(isfinite(data)); [B,err,anom_RHSurf(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    data = all.mmw(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_mmw(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

    if iOLR > 0
iYes2m = -1;
iYes2m = +1;
      if iYes2m > 0
        data = all.d2m(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_d2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.t2m(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_t2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.RH2m(thetimeSeason,ii);    k = find(isfinite(data)); [B,err,anom_RH2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.e2a(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1); 
        data = all.e2a(thetimeSeason,ii)/nanmean(all.e2a(thetimeSeason,ii));     
                                  k = find(isfinite(data)); [B,err,anom_frac_e2m(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = all.ecs(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_ecs(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
        data = all.Rld(thetimeSeason,ii);     k = find(isfinite(data)); [B,err,anom_Rld(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);

      end        
      data = all.olr(thetimeSeason,ii);       k = find(isfinite(data)); [B,err,anom_olr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.olr_clr(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_olr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr(thetimeSeason,ii);       k = find(isfinite(data)); [B,err,anom_ilr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr_clr(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_ilr_clr(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
      data = all.ilr_adj(thetimeSeason,ii);   k = find(isfinite(data)); [B,err,anom_ilr_adj(ii,:)] = compute_anomaly_wrapper(k,dayOFtime,data,4,-1,-1);
    end
  end

end

fprintf(1,'\n')
warning on
