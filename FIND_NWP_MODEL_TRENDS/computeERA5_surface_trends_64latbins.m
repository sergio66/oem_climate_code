disp(' ')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp('computeERA5_surface_trends.m is calling /home/sergio/MATLABCODE/FIND_TRENDS/Math_tsfit_robust_filter.m with four arguments ==> only keeping POSITIVE numbers')
disp(' ')

%% see eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/fit_robust_one_lat.m for lag 1
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /asl/matlib/h4tools

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing surface64 trends ... ')

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

  data = reshape(all.stemp,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  boo = usetime;
  boo = find(iaIndexUse' == 1);

  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_stemp(ii) = B(2);  trend64_stemp_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_stemp_lag(ii) = l(1);
   else
    trend64_stemp_lag(ii) = NaN;
   end
    
  data = reshape(all.TwSurf,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_TwSurf(ii) = B(2); trend64_TwSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_TwSurf_lag(ii) = l(1);
   else
    trend64_TwSurf_lag(ii) = NaN;
   end

  data = reshape(all.RHSurf,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii); 
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_RHSurf(ii) = B(2); trend64_RHSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_RHSurf_lag(ii) = l(1);
   else
    trend64_RHSurf_lag(ii) = NaN;
   end

  data = reshape(all.mmw,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_mmw(ii) = B(2);    trend64_mmw_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_mmw_lag(ii) = l(1);
   else
    trend64_mmw_lag(ii) = NaN;
   end

  data = reshape(all.olr,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_olr(ii) = B(2);    trend64_olr_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_olr_lag(ii) = l(1);
   else
    trend64_olr_lag(ii) = NaN;
   end

  data = reshape(all.olr_clr,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_olr_clr(ii) = B(2);    trend64_olr_clr_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_olr_clr_lag(ii) = l(1);
   else
    trend64_olr_clr_lag(ii) = NaN;
   end

  data = reshape(all.ilr,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_ilr(ii) = B(2);    trend64_ilr_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_ilr_lag(ii) = l(1);
   else
    trend64_ilr_lag(ii) = NaN;
   end

  data = reshape(all.ilr_clr,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_ilr_clr(ii) = B(2);    trend64_ilr_clr_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_ilr_clr_lag(ii) = l(1);
   else
    trend64_ilr_clr_lag(ii) = NaN;
   end

  data = reshape(all.ilr_adj,i228,72,64); data = nanmean(maskLF.*data,2); data = data(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend64_ilr_adj(ii) = B(2);    trend64_ilr_adj_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend64_ilr_adj_lag(ii) = l(1);
   else
    trend64_ilr_adj_lag(ii) = NaN;
   end

end
fprintf(1,'\n')
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('doing 4608 surface trends to check ---- +=1000,x=100,.=10')
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+ \n')
  elseif mod(ii,100) == 0
    fprintf(1,'x')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  data = all.stemp(:,ii);   
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend_stemp(ii) = B(2);  trend_stemp_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_stemp_lag(ii) = l(1);
   else
    trend_stemp_lag(ii) = NaN;
   end
    
  data = all.TwSurf(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend_TwSurf(ii) = B(2); trend_TwSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_TwSurf_lag(ii) = l(1);
   else
    trend_TwSurf_lag(ii) = NaN;
   end

  data = all.RHSurf(:,ii);  
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend_RHSurf(ii) = B(2); trend_RHSurf_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_RHSurf_lag(ii) = l(1);
   else
    trend_RHSurf_lag(ii) = NaN;
   end

  data = all.mmw(:,ii);     
  [B, stats] = Math_tsfit_lin_robust(dayOFtime(usetime),data(usetime),4); trend_mmw(ii) = B(2);    trend_mmw_err(ii) = stats.se(2);
  k = remove_nan(data);
  if length(k) > 50
    l = xcorr(data(k),1,'coeff');
    trend_mmw_lag(ii) = l(1);
   else
    trend_mmw_lag(ii) = NaN;
   end

end
fprintf(1,'\n')
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); wah = reshape(trend_stemp,72,64); whos wah;     wah = wah.*xmaskLF; plot(1:64,nanmean(wah,1),1:64,trend64_stemp); title('Stemp Trend');
figure(2); wah = reshape(trend_stemp_lag,72,64); whos wah; wah = wah.*xmaskLF; plot(1:64,nanmean(wah,1),1:64,trend64_stemp_lag); title('Stemp Lag');

figure(3); wah = reshape(trend_mmw,72,64); whos wah;     wah = wah.*xmaskLF; plot(1:64,nanmean(wah,1),1:64,trend64_mmw); title('Mmw Trend');
figure(4); wah = reshape(trend_mmw_lag,72,64); whos wah; wah = wah.*xmaskLF; plot(1:64,nanmean(wah,1),1:64,trend64_mmw_lag); title('Mmw Lag');


