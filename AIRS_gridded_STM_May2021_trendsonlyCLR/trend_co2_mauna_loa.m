function [trend,unc,endyy] = trend_co2_mauna_loa()
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /home/sergio/MATLABCODE/TIME

moo = load('co2_mauna_loaa.txt');

iS = find(moo(:,1) == 2002 & moo(:,2) == 9);

for ii = 2003 : 2022
  iE = find(moo(:,1) == ii & moo(:,2) == 8);  
  boo = iS : iE;
  yyx = moo(boo,1);
  mmx = moo(boo,2);
  ddx = ones(size(yyx)) * 15;
  daysSince2002 = change2days(yyx,mmx,ddx,2002);
  data = moo(boo,4);
  [B, stats, err] = Math_tsfit_lin_robust(daysSince2002,data,4);

  trend(ii-2002) = B(2);
  unc(ii-2002)   = stats.se(2);
  endyy(ii-2002) = ii;
end

wah = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_Nyears_2021.mat');
%% mauno loa = 19.4721 N, 155.5922 W
boo   = find(wah.rlat0 > 19.47,1);
booM1 = boo-1;
woo   = wah.co2trend(boo,:);
wooM1 = wah.co2trend(booM1,:);
woo = (woo*0.1) + (wooM1*0.9);

errorbar(endyy,trend,unc); title('Mauna Loaa CO2 trend'); xlabel('Year ending'); ylabel('ppmv/yr'); xlim([min(endyy) max(endyy)])
errorbar(endyy-2002,trend,unc); title('Mauna Loaa CO2 trend'); xlabel('Num Years'); ylabel('ppmv/yr'); xlim([min(endyy-2002) max(endyy-2002)])

plot(1:19,woo,'x-',1:20,trend,'o-'); hl = legend('ESRL code','Mauna Loa','location','best'); title('CO2 trends'); xlabel('Number of years since 2002/09')
%plot(1:19,trend(1:19)./woo,'x-'); 

