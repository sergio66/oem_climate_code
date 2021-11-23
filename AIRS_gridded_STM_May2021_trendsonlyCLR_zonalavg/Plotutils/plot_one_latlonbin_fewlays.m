addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

iDir = 0;
if iDir == 0
  if driver.NorD == +1
    outputdir = ['Output/Quantile' num2str(driver.iQuantile,'%02d') '/'];
  elseif driver.NorD == -1
    outputdir = ['Output_Day/Quantile' num2str(driver.iQuantile,'%02d') '/'];
  end
end

for ii = 1 : 6
  fig_screen(ii);
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = pN./pD;
plays = flipud(plays(1:100));

ii = 1;
  filex = [outputdir '/test' num2str(driver.iLat) '.mat'];
  driver = load(filex);

for ii = 1 : length(driver.jacobian.wvjaclays_used)
  junk = driver.jacobian.wvjaclays_used{ii}-5;
  playsRET(ii) = mean(plays(junk));
end

if driver.NorD == +1
  era5 = load('../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_trends_desc.mat');
elseif driver.NorD == -1
  era5 = load('../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat');
end

for ii = 1
  filex = [outputdir '/test' num2str(driver.iLat) '.mat'];
  a = load(filex);

  chanset(:,ii) = zeros(1,2645);
  g = a.jacobian.chanset;
  chanset(g,ii) = 1;
 
  traceNstemp(:,ii) = a.oem.finalrates(1:6);
  wv_ret(:,ii)      = a.oem.finalrates(driver.jacobian.water_i);
  temp_ret(:,ii)    = a.oem.finalrates(driver.jacobian.temp_i);
  ozone_ret(:,ii)   = a.oem.finalrates(driver.jacobian.ozone_i);

  traceNstemp_sigs(:,ii) = a.oem.finalsigs(1:6);
  wv_ret_sigs(:,ii)      = a.oem.finalsigs(driver.jacobian.water_i);
  temp_ret_sigs(:,ii)    = a.oem.finalsigs(driver.jacobian.temp_i);
  ozone_ret_sigs(:,ii)   = a.oem.finalsigs(driver.jacobian.ozone_i);

  tz_ak40(ii,:,:) = a.oem.ak_temp;
  wv_ak40(ii,:,:) = a.oem.ak_water;
  o3_ak40(ii,:,:) = a.oem.ak_ozone;

  input_rates(:,ii) = a.rateset.rates;
  fit_to_rates(:,ii) = a.oem.fit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; semilogy(wv_ret,playsRET,era5.trend_gas_1(1:97,driver.iibin),plays(1:97),'linewidth',2); 
  set(gca,'ydir','reverse'); ylim([10 1000]); title('WV frac rate'); hl=legend('UMBC','ERA5','location','best');
figure(2); clf; semilogy(temp_ret,playsRET,era5.trend_ptemp(1:97,driver.iibin),plays(1:97),'linewidth',2); 
  set(gca,'ydir','reverse'); ylim([10 1000]); title('T(z) rate'); hl=legend('UMBC','ERA5','location','best');
fprintf(1,'iibin = %4i UMBC CO2rate/stemprate  = %8.6f ppm/yr %8.6f K/yr   ERA stemprate = %8.6f K/yr \n',driver.iibin,traceNstemp([1 6],1),era5.trend_stemp(driver.iibin))
