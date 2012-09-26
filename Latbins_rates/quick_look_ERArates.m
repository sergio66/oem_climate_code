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

figure(3); clf
  plot(nanmean(waterrate(tropics,:)),plays,'b--','linewidth',2); hold on
  plot(nanmean(ptemprate(tropics,:)),plays,'r--','linewidth',2); 
plot(nanmean(waterrate(tropics,:))+nanstd(waterrate(tropics,:)),plays,'c--','linewidth',2);
plot(nanmean(waterrate(tropics,:))-nanstd(waterrate(tropics,:)),plays,'c--','linewidth',2);
plot(nanmean(ptemprate(tropics,:))+nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
plot(nanmean(ptemprate(tropics,:))-nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
hold off
  set(gca,'ydir','reverse'); grid; title('ERA (b) WV frac/yr (r) Temp K/yr')
  %hl = legend('ERA WV','location','east'); set(hl,'fontsize',10)
  %end

figure(3); clf
  subplot(121)
  plot(nanmean(waterrate(tropics,:)),plays,'b--','linewidth',2); hold on
  plot(nanmean(waterrate(tropics,:))+nanstd(waterrate(tropics,:)),plays,'c--','linewidth',2);
  plot(nanmean(waterrate(tropics,:))-nanstd(waterrate(tropics,:)),plays,'c--','linewidth',2);
  hold off; title('ERA Water frac/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]);

  subplot(122)
  plot(nanmean(ptemprate(tropics,:)),plays,'r--','linewidth',2); hold on
  plot(nanmean(ptemprate(tropics,:))+nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
  plot(nanmean(ptemprate(tropics,:))-nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
  hold off; title('ERA Temp frac/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);

  %set(gca,'ydir','reverse'); grid; title('ERA (b) WV frac/yr (r) Temp K/yr')
  %hl = legend('ERA WV','location','east'); 
  %end



%hold on
%hold off

%figure(4)
%hold off
%  set(gca,'ydir','reverse'); grid; title('T K/yr')
%  hl = legend('OEM','ERA','location','east'); set(hl,'fontsize',10)

%{

chanset = jacobian.chanset;
%g = dogoodchan;
figure(5);
plot(f,input_rates(18,:),'b',f,fitted_rates(18,:),'r',...
      f,input_rates(18,:)-fitted_rates(18,:),'k',...
      f(chanset),input_rates(18,chanset)-fitted_rates(18,chanset),'ko')
  hl = legend('input','fit','bias=input-fit','location','north'); set(hl,'fontsize',10)
  axis([640 2780 -0.10 +0.10]); grid

wvrates.oem = wv(18,:);
wvrates.oem_d = dwv(18,:);
wvrates.era = nanmean(waterrate(tropics,:));
wvrates.plays = plays;

trates.oem = t(18,:);
trates.oem_d = dt(18,:);
trates.era = nanmean(ptemprate(tropics,:));
trates.plays = plays;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'rateset = %s \n',rateset.datafile);
fprintf(1,'jacset  = %s \n',jacobian.filename);

%}
