for ii = 1 : 4; figure(ii); clf; end

iLoad = +1;   %% older
iLoad = -1;   %% newer
if iLoad == 1
  load('/home/sergio/MATLABCODE/ROSES_2013_GRANT/NewClearPDFs/era_geo_rates_allsky_AIRIBRAD.mat');
else
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Mar29_2019_Clr/Desc/all_latbins_rates.mat'); 
end

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;

load sarta_chans_for_l1c.mat
f = f(ichan);

%---------------------------------------------------------------------------
% Plot Rates and Fit
figure(1); clf
g1 = driver.jacobian.chanset;
plot(f(g1),driver.rateset.rates(g1),'k-');
hold on;
plot(f(g1),driver.oem.fit(g1),'b-');
plot(f(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r-')
grid;
axis([min(f(g1)) max(f(g1)) -0.15 +0.15]);
title('AIRS'); 
hl=legend('Obs','Fit','Residual');
set(hl,'fontsize',14)
hold off

%---------------------------------------------------------------------------
water = driver.oem.finalrates(driver.jacobian.water_i);
watersigs = driver.oem.finalsigs(driver.jacobian.water_i);
temp = driver.oem.finalrates(driver.jacobian.temp_i);
tempsigs = driver.oem.finalsigs(driver.jacobian.temp_i);
ozone = driver.oem.finalrates(driver.jacobian.ozone_i);
ozonesigs = driver.oem.finalsigs(driver.jacobian.ozone_i);

plevs = load('Data/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

for ii = 1 : length(driver.jacobian.wvjaclays_used)
  junk = driver.jacobian.wvjaclays_used{ii}-5;
  playsRET(ii) = mean(plays(junk));
end

if iLoad == 1
  waterrate_ak = thestats.akwaterrate(:,1:97);
  waterratestd = thestats.waterratestd(:,1:97)/2;
elseif iLoad == -1
  waterrate_ak = thestats.waterrate(:,1:97);
  waterratestd = thestats.waterratestd(:,1:97)/2;
end
figure(2); clf
  shadedErrorBarYLog10(water,playsRET,watersigs,'bo-');
  hold on
  semilogy(water,log10(playsRET),'bo-');
  shadedErrorBarYLog10(waterrate_ak(ix,:),plays,waterratestd(ix,:),'rx-');
  hold off; 
  hl = title('AIRS(b) ERA(r) H2O fr/yr'); 
  set(hl,'fontsize',12); 
  set(gca,'ydir','reverse'); grid; axis([-0.01 +0.01 1 3]); 

if iLoad == 1
  ptemprate_ak = thestats.akptemprate(:,1:97);
  ptempratestd = thestats.ptempratestd(:,1:97)/2;
elseif iLoad == -1
  ptemprate_ak = thestats.ptemprate(:,1:97);
  ptempratestd = thestats.ptempratestd(:,1:97)/2;
end
figure(3); clf
  shadedErrorBarYLog10(temp,playsRET,tempsigs,'bo-');
  hold on
  semilogy(temp,log10(playsRET),'bo-');
  shadedErrorBarYLog10(ptemprate_ak(ix,:),plays,ptempratestd(ix,:),'rx-'); 
  hold off; 
  hl = title('AIRS(b) ERA(r) (K/yr)'); 
  set(hl,'fontsize',12); 
  set(gca,'ydir','reverse'); grid; axis([-0.2 +0.15 1 3]);

if iLoad == 1
  o3rate_ak = thestats.akozonerate(:,1:97);
  o3ratestd = thestats.ozoneratestd(:,1:97)/2;
elseif iLoad == -1
  o3rate_ak = thestats.ozonerate(:,1:97);
  o3ratestd = thestats.ozoneratestd(:,1:97)/2;
end
figure(4); clf
  shadedErrorBarYLog10(ozone,playsRET,ozonesigs,'bo-');
  hold on
  semilogy(ozone,log10(playsRET),'bo-');
  shadedErrorBarYLog10(o3rate_ak(ix,:),plays,o3ratestd(ix,:),'rx-'); 
  hold off; 
  title('AIRS (b) ERA (r) Ozone fr/yr')
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]);   
%---------------------------------------------------------------------------

for ii = 1 : 4
  figure(ii); plotaxis2;
end
