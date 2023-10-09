%% now load in CERES
ceres_fnameS = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF-TOA_Ed4.1_Subset_200209-202108.nc';  %% what I brought
ceresS = load_ceres_data(ceres_fnameS,+1);

ceres_fnameR = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF_Ed4.1_Subset_200209-202108.nc';      %% what Ryan suggests
ceresR = load_ceres_data(ceres_fnameR,-1);

ceres_fnameR = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF_Ed4.1_Subset_200209-202108.nc';      %% what Ryan suggests
ceres_fnameR = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF-TOA_Ed4.2_Subset_200209-202208.nc'; 
ceresR = load_ceres_data(ceres_fnameR,-1);

ceres = ceresR;

addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

bonk = findstr(ceres_fnameR,'.nc');
startD = ceres_fnameR(bonk-13:bonk-08); startY = str2num(startD(1:4));
stopD  = ceres_fnameR(bonk-06:bonk-01); stopY  = str2num(stopD(1:4));
iCnt = 0;

for yyx = 2002 : stopY
  mmS = 1; mmE = 12;
  if yyx == 2002
    mmS = 09;
  elseif yyx == stopY
    mmE = 08;
  end
  for ii = mmS : mmE
    iCnt = iCnt + 1;
    all.yy(iCnt) = yyx;
    all.mm(iCnt) = ii;
    all.dd(iCnt) = 15;
  end
end
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

for ii = 1 : 180
  data = ceres.lwdata(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_lw(ii) = B(2);  
    trend_ceres_lw_err(ii) = stats.se(2);
  else
    trend_ceres_lw(ii) = NaN;
    trend_ceres_lw_err(ii) = NaN;
  end

  data = ceres.lwdata_clr(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_lw_clr(ii) = B(2);  
    trend_ceres_lw_clr_err(ii) = stats.se(2);
  else
    trend_ceres_lw_clr(ii) = NaN;
    trend_ceres_lw_clr_err(ii) = NaN;
  end
end

for ii = 1 : 180
  data = ceres.swdata(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_sw(ii) = B(2);  
    trend_ceres_sw_err(ii) = stats.se(2);
  else
    trend_ceres_sw(ii) = NaN;
    trend_ceres_sw_err(ii) = NaN;
  end

  data = ceres.swdata_clr(ii,:);
  boo = find(isfinite(data));
  if length(boo) > 20
    [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),4);
    trend_ceres_sw_clr(ii) = B(2);  
    trend_ceres_sw_clr_err(ii) = stats.se(2);
  else
    trend_ceres_sw_clr(ii) = NaN;
    trend_ceres_sw_clr_err(ii) = NaN;
  end
end

trend_ceres_lat = ceres.lat;
plot(trend_ceres_lat,trend_ceres_lw,trend_ceres_lat,trend_ceres_lw_clr,'linewidth',2); plotaxis2; hl = legend('allsky','clrsky','location','best');
xlabel('Latitude'); ylabel('Trend'); title('W/m2/K/yr')

ceres_trend.fname            = ceres_fnameR;
ceres_trend.trend_lat        = trend_ceres_lat;
ceres_trend.trend_lw         = trend_ceres_lw;
ceres_trend.trend_lw_err     = trend_ceres_lw_err;
ceres_trend.trend_lw_clr     = trend_ceres_lw_clr;
ceres_trend.trend_lw_clr_err = trend_ceres_lw_clr_err;
ceres_trend.trend_sw         = trend_ceres_sw;
ceres_trend.trend_sw_err     = trend_ceres_sw_err;
ceres_trend.trend_sw_clr     = trend_ceres_sw_clr;
ceres_trend.trend_sw_clr_err = trend_ceres_sw_clr_err;

iSave = input('save ceres trends (-1/+1) : ');
if iSave > 0
  saver = ['save ceres_trends_' num2str(stopY-startY,'%02d') 'year.mat ceres_trend'];
  eval(saver);
end
