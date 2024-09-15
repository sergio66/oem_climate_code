disp(' ')
disp('computeERA5_surface_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_surface_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp('computeERA5_surface_anoms_64latbins.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/compute_anomaly_wrapper.m')
disp(' ')

%% see eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/fit_robust_one_lat.m for lag 1
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /asl/matlib/h4tools

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing surface64 trends and anoms ... ')

addpath /asl/matlib/science/
[salti,p.landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);

xmaskLF = ones(72,64);
if iAorOorL == 0
  xmaskLF = ones(72,64);
elseif iAorOorL == 1
  xmaskLF = ones(1,4608) * NaN;
  xmaskLF(p.landfrac == 0) = 1;
  xmaskLF = reshape(xmaskLF,72,64);
elseif iAorOorL == -1
  xmaskLF = ones(1,4608) * NaN;
  xmaskLF(p.landfrac == 1) = 1;
  xmaskLF = reshape(xmaskLF,72,64);
end

[i228,~] = size(all.stemp);
for ii = 1 : i228
  maskLF(ii,:,:) = xmaskLF;
end

warning off
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'.')
  end

  boo = usetime;
  boo = find(iaIndexUse' == 1);

  data = reshape(all.stemp,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_stemp(ii,:) = junkanom;
    
  data = reshape(all.TwSurf,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_TwSurf(ii,:) = junkanom;

  data = reshape(all.RHSurf,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii); 
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_RHSurf(ii,:) = junkanom;

  data = reshape(all.mmw,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_mmw(ii,:) = junkanom;

  data = reshape(all.olr,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_olr(ii,:) = junkanom;

  data = reshape(all.olr_clr,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_olr_clr(ii,:) = junkanom;

  data = reshape(all.ilr,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_ilr(ii,:) = junkanom;

  data = reshape(all.ilr_clr,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_ilr_clr(ii,:) = junkanom;

  data = reshape(all.ilr_adj,i228,72,64); data = squeeze(nanmean(maskLF.*data,2)); data = data(:,ii);  
  k = find(isfinite(data(usetime))); [B,err,junkanom] = compute_anomaly_wrapper(k,dayOFtime(usetime),data(usetime),4,-1,-1); anom64_ilr_adj(ii,:) = junkanom;

end
fprintf(1,'\n')
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

