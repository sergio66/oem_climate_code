addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB

clear all; 

months = [1 4 7 10];
months = [1:12];
iCnt = 0;
for yy = 2002:2021
  mmS = 1; mmE = 12;
  if yy == 2002
    mmS = 09;
  elseif yy == 2021
    mmE = 08;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    iaFound(iCnt) = 0;
    fname = ['/asl/s1/sergio/AIRS_L3/EVAN_MANNING_MONTHLY_CCR/stratrad_airs.cc.v7.' num2str(yy) '-' num2str(months(mm),'%02d') '.nc'];
    if exist(fname)
      iaFound(iCnt) = 1;
      saveyear(iCnt)  = yy;
      savemonth(iCnt) = months(mm);
      saveday(iCnt)   = 15;

      fprintf(1,'reading %3i %s \n',iCnt,fname);
      a = read_netcdf_lls(fname);
    
      %                           %% nchan satzen lon  lat L/O/A  Asc/Desc
      %size(a.l1b_airs.rad)       %% 2378     10   18   18   3       2
      %size(a.l1b_airs.rad_nobs)  %% 2378     10   18   18   3       2
    
      rads = squeeze(a.l1b_airs.rad(:,:,:,:,3,2));
      nobs  = squeeze(a.l1b_airs.rad_nobs(:,:,:,:,3,2));
    
      rads = squeeze(nanmean(rads,[2 3])); %% [2 3] = satzen/lon
      nobs = squeeze(nansum(nobs,[2 3]));  %% [2 3] = satzen/lon
    
      saverads(iCnt,:,:) = rads;
      savenobs(iCnt,:,:) = nobs;
    end
  end
end

saverads = permute(saverads,[2 3 1]);
savenobs = permute(savenobs,[2 3 1]);

tCCR = rad2bt(a.l1b_airs.wnum,saverads);

tglobalCCR = squeeze(nanmean(tCCR,2));
iCntMax = iCnt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(';kzgj')
days = (1:iCntMax)/12 * 001; %% so this is years
days = (1:iCntMax)/12 * 365; %% so this is days
good = find(iaFound > 0);

addpath /home/sergio/MATLABCODE/FIND_TRENDS
for ii = 1 : 2378
  data = double(tglobalCCR(ii,good));
  boo = find(isfinite(data));
  if length(boo) > 20
    [c d] = polyfit(days(good(boo))/365,data(boo),1);
    trend(ii) = c(1);
    const(ii) = c(2);
  else
    trend(ii) = NaN;
    const(ii) = NaN;
  end
end

for ii = 1 : 2378
  data = double(tglobalCCR(ii,good));
  boo = find(isfinite(data));
  if length(boo) > 20
    [c d] = Math_tsfit_lin_robust(days(good(boo)),data(boo),4);
    bettertrend(ii) = c(2);
    betterconst(ii) = c(1);
  else
    bettertrend(ii) = NaN;
    betterconst(ii) = NaN;
 end
end

l1ctrends = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR_zonalavg/strow_rates_with_unc.mat');

figure(1); plot(a.l1b_airs.wnum,const,a.l1b_airs.wnum,betterconst); title('N month mean BTobs');
  plotaxis2; axis([640 1640 200 280])
figure(2); plot(a.l1b_airs.wnum,trend,a.l1b_airs.wnum,bettertrend,l1ctrends.f,l1ctrends.airsobs); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('CCR straight line fit','CCR Annual Cycle fit','LIC for STM','location','best','fontsize',10); ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')
