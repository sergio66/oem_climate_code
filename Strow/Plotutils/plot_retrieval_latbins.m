load Data/airs_f
figure(1);
  g1 = driver.jacobian.chanset_used;
  plot(f(g1),driver.rateset.rates(g1),f(g1),driver.rateset.rates(g1),'k',f(g1),driver.oem.fit(g1),'r.-'); grid; 
    axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','DEFINITELY OBS','fits'); set(hl,'fontsize',10)
figure(2);
  g1 = driver.jacobian.chanset_used;
  plot(f(g1),driver.rateset.rates(g1),f(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r'); grid; axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','fits'); set(hl,'fontsize',10)
% if length(driver.oem.finalrates) == 200
%   figure(3)
%   plot(driver.oem.finalrates(7:103),1:97,driver.oem.finalrates(104:200),1:97,'r')
%   set(gca,'ydir','reverse');
%   title('AIRS (b) : WV frac/yr (r) T K/yr'); grid
% end

% this is plotting
ecmcloud_file=['Data/allsky_rates.mat' 
];
load(ecmcloud_file);
ptemp_rate_cloud=ptemp_rate(1:36,1:97);
ptemp_rate_cloud_err=ptemp_rate_err(1:36,1:97); 
gas1_rate_cloud=gas1_rate(1:36,1:97); 
gas1_rate_cloud_err=gas1_rate_err(1:36,1:97); 

merrafile = '/asl/s1/rates/clear/Oct2013_MERRA/overocean__lays_spanday01_profilerates_Oct5_2013_robust.mat';
merrarates = load(merrafile);

load Data/airsL3_v6_rates_stats_March2014.mat
l3temprate=thestats.ptemprate; 
l3waterrate=thestats.waterrate; 
l3tempstd=thestats.ptempratestd; 
l3waterstd=thestats.waterratestd; 

ecmfile = '/asl/s1/rates/clear/Aug2013/';
ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012.mat'];
load(ecmfile); 


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
  water = temp * 0;      water(driver.jacobian.Q1jacindex) = driver.oem.finalrates(7:7+length(driver.jacobian.Q1jacindex)-1);
  watersigs = temp * 0;  watersigs(driver.jacobian.Q1jacindex) = driver.oem.finalsigs(7:7+length(driver.jacobian.Q1jacindex)-1);
else
  return
end

plevs = load('Data/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

%% shadedErrorBarYLog(temp,plays,tempsigs*1e4,'bo-',1);

%figure(2); clf
%  subplot(121)
%  shadedErrorBarYLog(water,plays,watersigs,'bo-',1);
%  hold on
%  shadedErrorBarYLog(waterrate(ix,:),plays,waterratestd(ix,:),'rx-',1);
%  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
%  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

%  subplot(122)
%  shadedErrorBarYLog(temp,plays,tempsigs,'bo-',1);
%  hold on
%  shadedErrorBarYLog(ptemprate(ix,:),plays,ptempratestd(ix,:),'rx-',1);
%  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10); 
%  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);

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
figure(5); clf
  subplot(121)
  shadedErrorBarYLog10(water,plays,watersigs,'bo-');
  hold on
  shadedErrorBarYLog10(waterrate_ak(ix,:),plays,waterratestd(ix,:),'rx-');
%  shadedErrorBarYLog10(squeeze(merrarates.water_allpars(ix,:,2)),plays,squeeze(merrarates.water_allerrors(ix,:,2)),'gs-');
  shadedErrorBarYLog10(waterrate_merra_ak(ix,:),plays,squeeze(merrarates.water_allerrors(ix,:,2)),'gs-'); 
%  shadedErrorBarYLog10(gas1_rate_cloud(ix,:),plays,gas1_rate_cloud_err(ix,:),'k^-'); 
  shadedErrorBarYLog10(waterrate_cloud_ak(ix,:),plays,gas1_rate_cloud_err(ix,:),'k^-'); 
% L3 rates  shadedErrorBarYLog10(l3waterrate(ix,:),Qlevs,l3waterstd(ix,:),'mv-'); 
  hold off; hl = title('AIRS(b) ERA(r) MERRA(g) ERA-CLD(blk) H2O fr/yr'); set(hl,'fontsize',9); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]); 
  rms_wtr_strat=rms(water(1:49)'-waterrate(ix,1:49))./rms(waterrate(ix,1:49)); 
  rms_wtr_trop=rms(water(50:97)'-waterrate(ix,50:97))./rms(waterrate(ix,50:97));

  subplot(122)
  shadedErrorBarYLog10(temp,plays,tempsigs,'bo-');
  hold on
%  shadedErrorBarYLog10(ptemprate(ix,:),plays,ptempratestd(ix,:),'rx-');
  shadedErrorBarYLog10(temprate_ak(ix,:),plays,ptempratestd(ix,:),'rx-'); 
%  shadedErrorBarYLog10(squeeze(merrarates.ptemp_allpars(ix,:,2)),plays,squeeze(merrarates.ptemp_allerrors(ix,:,2)),'gs-');
  shadedErrorBarYLog10(temprate_merra_ak(ix,:),plays,squeeze(merrarates.ptemp_allerrors(ix,:,2)),'gs-'); 
%  shadedErrorBarYLog10(ptemp_rate_cloud(ix,:),plays,ptemp_rate_cloud_err(ix,:),'k^-');
  shadedErrorBarYLog10(temprate_cloud_ak(ix,:),plays,ptemp_rate_cloud_err(ix,:),'k^-'); 
%  shadedErrorBarYLog10(thestats.ptemprate(ix,:),Tlevs,ptempratestd(ix,1:24),'mv-'); 
%  shadedErrorBarYLog10(thestats.ptemprate(ix,:),Tlevs,thestats.ptempratestd(ix,:),'mv-');
% L3 temp profile rates  shadedErrorBarYLog10(l3temprate(ix,:),Tlevs,l3tempstd(ix,:),'mv-')
  hold off; hl = title('AIRS(b) ERA(r) MERRA(g) ERA-CLD(blk) T (K/yr)'); set(hl,'fontsize',9); 
  set(gca,'ydir','reverse'); grid; axis([-0.2 +0.15 1 3]);
  rms_tmp_strat=rms(temp(1:49)'-ptemprate(ix,1:49))./rms(ptemprate(ix,1:49)) ;
  rms_tmp_trop=rms(temp(50:97)'-ptemprate(ix,50:97))./rms(ptemprate(ix,50:97));
%'test 9' 
%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
end

if strcmp(user_name(1:end-1),'sergio')
  hz = p2h(plays)/1000; %% change to km
  firstderW = gradient(water,hz); secondderW = gradient(firstderW,hz);
  firstderT = gradient(temp,hz);  secondderT = gradient(firstderT,hz);
  p50  = find(plays >= 50);
  p100 = find(plays >= 100);

  W1x             = zerocross(firstderW);   T1x             = zerocross(firstderT); 
  [MaxW1x,MinW1x] = localmaxmin(firstderW); [MaxT1x,MinT1x] = localmaxmin(firstderT);

  second_deriv = sqrt([nanmean(secondderW(p50).*secondderW(p50)) nanmean(secondderT(p50).*secondderT(p50))]);
  fprintf(1,'  ---->>>> 50mb-1000mb average second deriv (squared) of W, T = %8.6f %8.6f \n',second_deriv);

  % figure(3); clf; plot(secondderW,1:97,secondderT,1:97,'r')
  % figure(3); clf; semilogy(secondderW,plays,secondderT,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  % figure(3); clf; semilogy(firstderW,plays,secondderW,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  % figure(3); clf; semilogy(firstderT,plays,secondderT,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  figure(6); clf; 
  subplot(121); semilogy(temp,plays)
    set(gca,'ydir','reverse'); grid; axis([-0.1 +0.1 20 1000]);
    hl = legend('T(z)','Location','northeast'); set(hl,'Fontsize',8);
  subplot(122); semilogy(firstderT,plays,...
                         secondderT,plays,'r',...
                         zeros(size(T1x)),plays(ceil(T1x)),'bo',...
                         firstderT(MaxT1x),plays(MaxT1x),'bd',...
                         firstderT(MinT1x),plays(MinT1x),'bs')

    set(gca,'ydir','reverse'); grid; axis([-0.1 +0.1 20 1000]);
    hl = legend('dT/dz','d2T/dz2','Location','northeast'); set(hl,'Fontsize',8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(4); clf; imagesc(log10(abs(log10(abs(driver.oem.inv_se))))); caxis([-1 2]); colorbar

