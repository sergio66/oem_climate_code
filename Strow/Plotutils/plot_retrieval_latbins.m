ix = driver.iibin;
%---------------------------------------------------------------------------
% Plot Rates and Fit
figure(1)
g1 = driver.jacobian.chanset;
plot(f(g1),driver.rateset.rates(g1),'k-');
hold on;
plot(f(g1),driver.oem.fit(g1),'b-');
plot(f(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r-')
grid;
axis([500 3000 -0.15 +0.15]);
title('AIRS'); 
hl=legend('Obs','Fit','Residual');
set(hl,'fontsize',14)
%---------------------------------------------------------------------------
% Plot profiles, fit and ERA, MERRA, ERA-allsky
load('Data/allsky_rates.mat');
ptemp_rate_cloud = ptemp_rate(1:36,1:97);
ptemp_rate_cloud_err = ptemp_rate_err(1:36,1:97); 
gas1_rate_cloud = gas1_rate(1:36,1:97); 
gas1_rate_cloud_err = gas1_rate_err(1:36,1:97); 

merrarates = load('/asl/s1/rates/Clear/Merra_rates/overocean__lays_spanday01_profilerates_Oct5_2013_robust.mat');
load Data/airsL3_v6_rates_stats_March2014.mat

l3temprate = thestats.ptemprate; 
l3waterrate = thestats.waterrate; 
l3tempstd = thestats.ptempratestd; 
l3waterstd = thestats.waterratestd; 

load('/asl/s1/rates/Clear/Profile_rates/overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012.mat');

if length(driver.oem.finalrates) == 200
  water = driver.oem.finalrates(7:103);
  watersigs = driver.oem.finalsigs(7:103);
  temp = driver.oem.finalrates(104:200);
  tempsigs = driver.oem.finalsigs(104:200);
elseif length(driver.oem.finalrates) == 103 & length(driver.jacobian.Q1jacindex) < 1
  temp = driver.oem.finalrates(7:103);
  tempsigs = driver.oem.finalsigs(7:103);
  water = temp * 0;
  watersigs = temp*0;
elseif length(driver.oem.finalrates) > 103 & length(driver.jacobian.Q1jacindex) < 97
  bonk = length(driver.jacobian.Q1jacindex);
  temp = driver.oem.finalrates((7:103)+bonk);
  tempsigs = driver.oem.finalsigs((7:103)+bonk);
  water = temp * 0; 
  water(driver.jacobian.Q1jacindex) = driver.oem.finalrates(7:7+length(driver.jacobian.Q1jacindex)-1);
  watersigs = temp * 0;  
  watersigs(driver.jacobian.Q1jacindex) = driver.oem.finalsigs(7:7+length(driver.jacobian.Q1jacindex)-1);
else
  return
end

plevs = load('Data/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

% Normalize averaging kernel:
for i=1:97 
   sum_water=0.0; 
   sum_temp=0.0; 
   for j=1:97 
      sum_water=sum_water+driver.oem.ak_water(j,i); 
      sum_temp=sum_temp+driver.oem.ak_temp(j,i); 
   end 
   for j=1:97
      ak_water(j,i)=driver.oem.ak_water(j,i);%/sum_water; 
      ak_temp(j,i)=driver.oem.ak_temp(j,i);%/sum_temp; 
   end 
end 

%  Apply averaging kernel to waterrate:
for i=1:97
   sum_era=0;
   sum_cloud=0; 
   sum_merra=0; 
   sum_era_temp=0; 
   sum_cloud_temp=0; 
   sum_merra_temp=0; 
   for j=1:97
      sum_era=sum_era+waterrate(ix,j)*ak_water(i,j);
      sum_cloud=sum_cloud+gas1_rate_cloud(ix,j)*ak_water(i,j); 
      sum_merra=sum_merra+merrarates.water_allpars(ix,j,2)*ak_water(i,j); 
      sum_era_temp=sum_era_temp+ptemprate(ix,j)*ak_temp(i,j); 
      sum_cloud_temp=sum_cloud_temp+ptemp_rate_cloud(ix,j)*ak_temp(i,j); 
      sum_merra_temp=sum_merra_temp+merrarates.ptemp_allpars(ix,j,2)*ak_temp(i,j); 
   end
   waterrate_ak(ix,i)=sum_era;                        
   waterrate_merra_ak(ix,i)=sum_merra; 
   waterrate_cloud_ak(ix,i)=sum_cloud; 
   temprate_ak(ix,i)=sum_era_temp; 
   temprate_cloud_ak(ix,i)=sum_cloud_temp; 
   temprate_merra_ak(ix,i)=sum_merra_temp; 
end

figure(2)
  subplot(121)
  shadedErrorBarYLog10(water,plays,watersigs,'bo-');
  hold on
  shadedErrorBarYLog10(waterrate_ak(ix,:),plays,waterratestd(ix,:),'rx-');
%   shadedErrorBarYLog10(waterrate_merra_ak(ix,:),plays,squeeze(merrarates.water_allerrors(ix,:,2)),'gs-'); 
  shadedErrorBarYLog10(waterrate_cloud_ak(ix,:),plays,gas1_rate_cloud_err(ix,:),'g^-'); 
  % L3 Q 
  shadedErrorBarYLog10(l3waterrate(ix,:),Qlevs,0.001*l3waterstd(ix,:),'ko-'); 
  hold off; 
%  hl = title('AIRS(b) ERA(r) MERRA(g) ERA-CLD(blk) H2O fr/yr'); 
  hl = title('AIRS(b) ERA(r) AIRS:L3(blk) ERA:allsky(g) H2O fr/yr'); 
  set(hl,'fontsize',12); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]); 
  rms_wtr_strat=rms(water(1:49)'-waterrate(ix,1:49))./rms(waterrate(ix,1:49)); 
  rms_wtr_trop=rms(water(50:97)'-waterrate(ix,50:97))./rms(waterrate(ix,50:97));
  subplot(122)
  shadedErrorBarYLog10(temp,plays,tempsigs,'bo-');
  hold on
  shadedErrorBarYLog10(temprate_ak(ix,:),plays,ptempratestd(ix,:),'rx-'); 
%   shadedErrorBarYLog10(temprate_merra_ak(ix,:),plays,squeeze(merrarates.ptemp_allerrors(ix,:,2)),'gs-'); 
  shadedErrorBarYLog10(temprate_cloud_ak(ix,:),plays,ptemp_rate_cloud_err(ix,:),'g^-'); 
  % L3 T
  shadedErrorBarYLog10(l3temprate(ix,:),Tlevs,0.001*l3tempstd(ix,:),'ko-')
  hold off; 
%  hl = title('AIRS(b) ERA(r) MERRA(g) ERA-CLD(blk) T (K/yr)'); 
  hl = title('AIRS(b) ERA(r) AIRS:L3(blk) ERA:allsky(g) T (K/yr)'); 
  set(hl,'fontsize',12); 
  set(gca,'ydir','reverse'); grid; axis([-0.2 +0.15 1 3]);
  rms_tmp_strat=rms(temp(1:49)'-ptemprate(ix,1:49))./rms(ptemprate(ix,1:49)) ;
  rms_tmp_trop=rms(temp(50:97)'-ptemprate(ix,50:97))./rms(ptemprate(ix,50:97));
%---------------------------------------------------------------------------
