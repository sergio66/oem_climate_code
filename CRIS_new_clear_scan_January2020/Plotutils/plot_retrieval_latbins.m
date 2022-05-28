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
if length(driver.oem.finalrates) == 200
  water = driver.oem.finalrates(7:103);
  watersigs = driver.oem.finalsigs(7:103);
  temp = driver.oem.finalrates(104:200);
  tempsigs = driver.oem.finalsigs(104:200);
elseif length(driver.oem.finalrates) == 296
  water = driver.oem.finalrates(6:102);
  watersigs = driver.oem.finalsigs(6:102);
  temp = driver.oem.finalrates(103:199);
  tempsigs = driver.oem.finalsigs(103:199);
  ozone = driver.oem.finalrates(200:296);
  ozonesigs = driver.oem.finalsigs(200:296);
elseif length(driver.oem.finalrates) == 300
  water = driver.oem.finalrates(10:106);
  watersigs = driver.oem.finalsigs(10:106);
  temp = driver.oem.finalrates(107:203);
  tempsigs = driver.oem.finalsigs(107:203);
  ozone = driver.oem.finalrates(204:300);
  ozonesigs = driver.oem.finalsigs(204:300);
end

plevs = load('Data/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

if iLoad == 1
  waterrate_ak = thestats.akwaterrate(:,1:97);
  waterratestd = thestats.waterratestd(:,1:97)/2;
elseif iLoad == -1
  waterrate_ak = thestats.waterrate(:,1:97);
  waterratestd = thestats.waterratestd(:,1:97)/2;
end
figure(2); clf
  shadedErrorBarYLog10(water,plays,watersigs,'bo-');
  hold on
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
  shadedErrorBarYLog10(temp,plays,tempsigs,'bo-');
  hold on
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
  shadedErrorBarYLog10(ozone,plays,ozonesigs,'bo-');
  hold on
  shadedErrorBarYLog10(o3rate_ak(ix,:),plays,o3ratestd(ix,:),'rx-'); 
  hold off; 
  title('AIRS (b) ERA (r) Ozone fr/yr')
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]);   
%---------------------------------------------------------------------------

for ii = 1 : 4
  figure(ii); plotaxis2;
end
