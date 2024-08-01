for ii = 1 : 4; figure(ii); clf; end

iLoad = +1;   %% older
iLoad = -1;   %% newer
iLoad = 0;    %% neweest

if iLoad == 1
  era = load('/home/sergio/MATLABCODE/ROSES_2013_GRANT/NewClearPDFs/era_geo_rates_allsky_AIRIBRAD.mat');
elseif iLoad == -1
  era = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Mar29_2019_Clr/Desc/all_latbins_rates.mat'); 
elseif iLoad == 0
  era5 = load('../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_desc_64latbins.mat'); ahwoo = driver.iLat;   %% and use era5.trend64
  era5 = load('../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat'); ahwoo = driver.iibin;            %% and use era5.trend
end

if topts.dataset < 30
  hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
  vchan2834 = hdfread(hdffile,'freq');
  f = vchan2834;

  load sarta_chans_for_l1c.mat
  f = f(ichan);
else
  f = aux.f;
end

%---------------------------------------------------------------------------
% Plot Rates and Fit
figure(1); clf

g1 = driver.jacobian.chanset;
plot(f(g1),driver.rateset.rates(g1),'k-','linewidth',2);
hold on;
plot(f(g1),driver.oem.fit(g1),'b-','linewidth',2);
plot(f(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r-','linewidth',2);
plot(f(g1),+driver.rateset.unc_rates(g1),'color',[1 1 1]*0.75);
plot(f(g1),-driver.rateset.unc_rates(g1),'color',[1 1 1]*0.75);
grid;
axis([min(f(g1)) max(f(g1)) -0.1 +0.1]);
plotaxis2;
hl=legend('Obs','Fit','Residual','location','best');
set(hl,'fontsize',10)
hold off
if topts.dataset < 30
  title('AIRS'); xlabel('Wavenumber cm-1')
else
  title('AMSU'); xlabel('Grequency GHz')
end

if isfield(driver.oem,'spectral_deltan00')

  plot(f(g1),driver.rateset.rates(g1),'k-','linewidth',2);
  hold on;
  plot(f(g1),driver.oem.fit(g1),'b-','linewidth',2);
  plot(f(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r-','linewidth',2);
  plot(f(g1),driver.oem.spectral_deltan00(g1),'g.-');

  plot(f(g1),+driver.rateset.unc_rates(g1),'color',[1 1 1]*0.75);
  plot(f(g1),-driver.rateset.unc_rates(g1),'color',[1 1 1]*0.75);
  grid;
  axis([min(f(g1)) max(f(g1)) -0.1 +0.1]);
  title('AIRS'); 
  hold off

  plotaxis2;
  hl=legend('Obs','Fit','Residual','No TraceGas Obs','location','best');
  set(hl,'fontsize',10)
end

if topts.dataset < 30
  title('AIRS'); xlabel('Wavenumber cm-1')
else
  title('AMSU'); xlabel('Grequency GHz')
end

%---------------------------------------------------------------------------
water = driver.oem.finalrates(driver.jacobian.water_i);
watersigs = driver.oem.finalsigs(driver.jacobian.water_i);
temp = driver.oem.finalrates(driver.jacobian.temp_i);
tempsigs = driver.oem.finalsigs(driver.jacobian.temp_i);
ozone = driver.oem.finalrates(driver.jacobian.ozone_i);
ozonesigs = driver.oem.finalsigs(driver.jacobian.ozone_i);

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plays = plevsA./plevsB;
%plays = plevs(4:100); 
plays = flipud(plays);

clear *_ak* *std*

for ii = 1 : length(driver.jacobian.wvjaclays_used)
  junk = driver.jacobian.wvjaclays_used{ii}-6;
  playsRET(ii) = mean(plays(junk));
end

clear waterrate_ak0 waterratestd0 ptemprate_ak0 ptempratestd0 o3rate_ak0 o3rate_std0

if abs(iLoad) <= 1
  if iLoad == 1
    waterrate_ak0 = era.thestats.akwaterrate(:,1:97);
    waterratestd0 = era.thestats.waterratestd(:,1:97)/2;
  elseif iLoad == -1
    waterrate_ak0 = era.thestats.waterrate(:,1:97);
    waterratestd0 = era.thestats.waterratestd(:,1:97)/2;
  elseif iLoad == 0
    waterrate_ak0 = ones(ix,1)*era5.trend_gas_1(1:100,ahwoo)';
    waterratestd0 = ones(ix,1)*era5.trend_gas_1_err(1:100,ahwoo)';
  end
  for ii = 1 : length(driver.jacobian.wvjaclays_used)
    junk = driver.jacobian.wvjaclays_used{ii}-6;
    waterrate_ak1(:,ii) = mean(waterrate_ak0(:,junk)');
    waterratestd1(:,ii) = mean(waterratestd0(:,junk)');
  end
  ak = driver.oem.ak_water;
  waterrate_akF(ix,:) = (ak * waterrate_ak1(ix,:)')';
  waterratestdF(ix,:) = (ak * waterratestd1(ix,:)')';
  playsRETW = playsRET(1:length(water));
  figure(2); clf
    shadedErrorBarYLog10(water,playsRETW,watersigs,'bo-');
    hold on
    semilogy(water,log10(playsRETW),'bo-');
    shadedErrorBarYLog10(waterrate_ak1(ix,:),playsRETW,waterratestd1(ix,:),'rd-');
    shadedErrorBarYLog10(waterrate_akF(ix,:),playsRETW,waterratestdF(ix,:),'gx-');
    hold off; 
    ax = axis; line([ax(1) ax(2)],log10([aux.spres  aux.spres]),'color','g','linewidth',2)
    ax = axis; line([ax(1) ax(2)],log10([aux.trop_P aux.trop_P]),'color','c','linewidth',2)
    title('(b)UMBC (r)ERA (g)AK*ERA WVfrac(z) 1/yr','fontsize',12);
    set(gca,'ydir','reverse'); grid; axis([-0.01 +0.01 2 3]); 
  
  if iLoad == 1
    ptemprate_ak0 = era.thestats.akptemprate(:,1:97);
    ptempratestd0 = era.thestats.ptempratestd(:,1:97)/2;
  elseif iLoad == -1
    ptemprate_ak0 = era.thestats.ptemprate(:,1:97);
    ptempratestd0 = era.thestats.ptempratestd(:,1:97)/2;
  elseif iLoad == 0
    ptemprate_ak0 = ones(ix,1)*era5.trend_ptemp(1:100,ahwoo)';
    ptempratestd0 = ones(ix,1)*era5.trend_ptemp_err(1:100,ahwoo)';
  end
  for ii = 1 : length(driver.jacobian.wvjaclays_used)
    junk = driver.jacobian.wvjaclays_used{ii}-6;
    ptemprate_ak1(:,ii) = mean(ptemprate_ak0(:,junk)');
    ptempratestd1(:,ii) = mean(ptempratestd0(:,junk)');
  end
  ak = driver.oem.ak_temp;
  ptemprate_akF(ix,:) = (ak * ptemprate_ak1(ix,:)')';
  ptempratestdF(ix,:) = (ak * ptempratestd1(ix,:)')';
  figure(3); clf
    shadedErrorBarYLog10(temp,playsRETW,tempsigs,'bo-');
    hold on
    semilogy(temp,log10(playsRETW),'bo-');
    shadedErrorBarYLog10(ptemprate_ak1(ix,:),playsRETW,ptempratestd1(ix,:),'rd-'); 
    shadedErrorBarYLog10(ptemprate_akF(ix,:),playsRETW,ptempratestdF(ix,:),'gx-'); 
    hold off; 
    ax = axis; line([ax(1) ax(2)],log10([aux.spres  aux.spres]),'color','g','linewidth',2)
    ax = axis; line([ax(1) ax(2)],log10([aux.trop_P aux.trop_P]),'color','c','linewidth',2)
    hl = title('(b)UMBC (r)ERA (g)AK*ERA T(z) (K/yr)','fontsize',12); 
    set(gca,'ydir','reverse'); grid; axis([-0.2 +0.15 1 3]);
  
  if iLoad == 1
    o3rate_ak0 = era.thestats.akozonerate(:,1:97);
    o3ratestd0 = era.thestats.ozoneratestd(:,1:97)/2;
  elseif iLoad == -1
    o3rate_ak0 = era.thestats.ozonerate(:,1:97);
    o3ratestd0 = era.thestats.ozoneratestd(:,1:97)/2;
  elseif iLoad == 0
    o3rate_ak0 = ones(ix,1)*era5.trend_gas_3(1:100,ahwoo)';
    o3ratestd0 = ones(ix,1)*era5.trend_gas_3_err(1:100,ahwoo)';
  end
  for ii = 1 : length(driver.jacobian.wvjaclays_used)
    junk = driver.jacobian.wvjaclays_used{ii}-6;
    o3rate_ak1(:,ii) = mean(o3rate_ak0(:,junk)');
    o3ratestd1(:,ii) = mean(o3ratestd0(:,junk)');
  end
  ak = driver.oem.ak_ozone;
  o3rate_akF(ix,:) = (ak * o3rate_ak1(ix,:)')';
  o3ratestdF(ix,:) = (ak * o3ratestd1(ix,:)')';
  figure(4); clf
    shadedErrorBarYLog10(ozone,playsRETW,ozonesigs,'bo-');
    hold on
    semilogy(ozone,log10(playsRETW),'bo-');
    shadedErrorBarYLog10(o3rate_ak1(ix,:),playsRETW,o3ratestd1(ix,:),'rd-'); 
    shadedErrorBarYLog10(o3rate_akF(ix,:),playsRETW,o3ratestdF(ix,:),'gx-'); 
    hold off; 
    ax = axis; line([ax(1) ax(2)],log10([aux.spres  aux.spres]),'color','g','linewidth',2)
    ax = axis; line([ax(1) ax(2)],log10([aux.trop_P aux.trop_P]),'color','c','linewidth',2)
    title('(b)UMBC (r)ERA (g)AK*ERA O3(z) fr/yr','fontsize',12)
    set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]);   
  %---------------------------------------------------------------------------

else
  figure(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 2 : 4
  figure(ii); plotaxis2;
end
