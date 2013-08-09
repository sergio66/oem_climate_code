%%   MakeERA_ratespectra  : just uses linear trends to add jacobians together
%%   MakeERA_ratespectra2 : puts in "time series" of gas_1,gas_3, ptemp and stemp, based on linear trends
%%   MakeERA_ratespectra3 : puts in "time series" of gas_1,gas_3, ptemp and stemp, based on linear trends AND ADDS IN variability

addpath ../Latbins_rates/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ecmrates = input('enter (1) 2002-2011 or (2) 2007-2011  or (3) 2007-2012  ERA rates : ');
ecmfile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/';
if ecmrates == 1
  co2str = '_07_2002_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era2_lays_profilerates_May23_2011_robust.mat'];
elseif ecmrates == 2
  %ecmfile = [ecmfile ...
  %    'overocean_gsx_1dayV1_ecmwf2_lays_spanday16_profilerates_Aug18_2011_robust.mat'];
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era_lays_spanday01_profilerates_Aug26_2011_robust.mat'];
elseif ecmrates == 3
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Aug28_2012_robust_span_07_2007_07_2012.mat'];
end
load(ecmfile);

disp('tracegas times : (1) 8yr 07/2002-07/2010 (2) 8 yr 09/2003-09/2011 (3) 9 yr 09/2002-09/2011');
disp('tracegas times : (1) 8yr 09/2002-08/2010 (2) 7 yr 09/2003-08/2010')
tracegasrates = input('Enter tracegas rates : ');
if tracegasrates == 1
  co2str = '_07_2002_07_2010';
  co2str = '_09_2002_08_2010';
elseif tracegasrates == 2
  co2str = '_09_2003_09_2011';
  co2str = '_09_2003_08_2010';
elseif tracegasrates == 3
  co2str = '_09_2003_09_2011';
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsA = plevs(1:100)-plevs(2:101);
playsB = log(plevs(1:100)./plevs(2:101));
plays = playsA./playsB;
plays = plays(4:100);
plays = flipud(plays);

load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era3_lays_spanday01_save_lat_Sep7_2011.mat
tropics = find(abs(save_lat) <= 30);

figure(1); clf
  subplot(121)
  shadedErrorBarY(nanmean(waterrate(tropics,:)),plays,nanstd(waterrate(tropics,:)),'rx-',1);
  hold off; title('ERA(r) Water frac/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]);

  subplot(122)
  shadedErrorBarY(nanmean(ptemprate(tropics,:)),plays,nanstd(ptemprate(tropics,:)),'rx-',1);
  hold off; title('ERA(r) Temp K/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);


figure(2); clf
  pcolor(save_lat,plays,double(waterrate')); shading interp; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA WV frac/yr ')
figure(3); clf
  pcolor(save_lat,plays,double(ptemprate')); shading interp; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA T K/yr ')

figure(2); caxis([-0.02 +0.02]); colorbar
figure(3); caxis([-0.1 +0.1]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../../oem_pkg/Test/M_TS_jac_all.mat

%% see str1 and str2
CO2 = 2.2;  rate(1) = CO2;
O3  = 0.01; rate(2) = O3;
N2O = 1.0;  rate(3) = N2O;
CH4 = 5.0;  rate(4) = CH4;
CFC = 02;   rate(5) = CFC;

amp = 0.1;
for ii = 1 : 36
  m_ts_jac = squeeze(M_TS_jac_all(ii,:,:));
  thefit = zeros(1,2378);
  coeffs = zeros(1,200);
  coeffs(1:5) = rate;
  coeffs(006)     = stemprate(ii) + amp*(rand(1,1)-0.5)*stempratestd(ii);
  coeffs(007:103) = waterrate(ii,:) + amp*(rand(1,97)-0.5).*waterratestd(ii,:);
  coeffs(104:200) = ptemprate(ii,:) + amp*(rand(1,97)-0.5).*ptempratestd(ii,:);
  for ix = 1 : length(coeffs)
    thefit = thefit + coeffs(ix)/qrenorm(ix)*m_ts_jac(:,ix)';
  end
  simrate(ii,:) = thefit;
  b_obs(ii,:,2)  = thefit;
  b_obs(ii,:,10) = zeros(size(thefit));
end;

addpath /strowdata1/shared/sergio/MATLABCODE/oem_pkg_run/Latbins_rates/
g = dogoodchan;
airs = load('../../oem_pkg/Test/fitout_8year_v4_robust.mat');
figure(5); plot(f(g),squeeze(airs.b_obs(:,g,2))); title('ACTUAL OBS')

figure(6); plot(f(g),squeeze(airs.b_obs(:,g,2)) - simrate(:,g)); title('OBS - CAL') 

b_err_obs = airs.b_err_obs;
figure(4); plot(f,simrate); title('SIMS');
strx = 'see /home/sergio/MATLABCODE/oem_pkg_run/MakeERA_ratespectra/generate_era_rate_spectra.m';
save syntheticERArates.mat b_err_obs b_obs rate stemprate waterrate ptemprate stempratestd waterratestd ptempratestd amp strx ecmrates
