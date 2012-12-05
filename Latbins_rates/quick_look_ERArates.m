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

noaa_str = '2002_09_2012_08.mat';
nTLS = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/NOAA_STAR/TLS_2002_09_2012_08.mat');
nTMT = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/NOAA_STAR/TMT_2002_09_2012_08.mat');
nTUT = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/NOAA_STAR/TUT_2002_09_2012_08.mat');

% these are 30 yr trends
% rssTLT = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU_MONTHLY_GRIDS/zonal_trends_tlt_v03_3.txt');
% rssTMT = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU_MONTHLY_GRIDS/zonal_trends_tmt_v03_3.txt');
% rssTTT = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU_MONTHLY_GRIDS/zonal_trends_ttt_v03_3.txt');
% rssTLS = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU_MONTHLY_GRIDS/zonal_trends_tls_v03_3.txt');
% rssTTS = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU_MONTHLY_GRIDS/zonal_trends_tts_v03_3.txt');

rssTRP = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/rss_amsu_rateNEW_trp_2002_9_2012_8.mat');
rssNML = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/rss_amsu_rateNEW_nml_2002_9_2012_8.mat');
rssSML = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/rss_amsu_rateNEW_sml_2002_9_2012_8.mat');
rssNPL = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/rss_amsu_rateNEW_npl_2002_9_2012_8.mat');
rssSPL = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/rss_amsu_rateNEW_spl_2002_9_2012_8.mat');

rssPRESS_tls = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU/std_atmosphere_wt_function_chan_tls.txt');
rssPRESS_tts = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU/std_atmosphere_wt_function_chan_tts.txt');
rssPRESS_tlt = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU/std_atmosphere_wt_function_chan_tlt_ocean.txt');
rssPRESS_tmt = load('/home/sergio/MATLABCODE/RATES_GLOBALVIEW/RSS_AMSU/std_atmosphere_wt_function_chan_tmt_ocean.txt');

xPRESS_tls = find(rssPRESS_tls(:,6) == max(rssPRESS_tls(:,6))); xPRESS_tls = rssPRESS_tls(xPRESS_tls(1),4)/100;
xPRESS_tls = sum(rssPRESS_tls(:,6) .* rssPRESS_tls(:,4))/sum(rssPRESS_tls(:,6))/100;
  plot(rssPRESS_tls(:,6),rssPRESS_tls(:,4)/100);
%ret

xPRESS_tts = find(rssPRESS_tts(:,6) == max(rssPRESS_tts(:,6))); xPRESS_tts = rssPRESS_tts(xPRESS_tts(1),4)/100;
xPRESS_tts = sum(rssPRESS_tts(:,6) .* rssPRESS_tts(:,4))/sum(rssPRESS_tts(:,6))/100;
  plot(rssPRESS_tts(:,6),rssPRESS_tts(:,4)/100);
%ret

xPRESS_tmt = find(rssPRESS_tmt(:,6) == max(rssPRESS_tmt(:,6))); xPRESS_tmt = rssPRESS_tmt(xPRESS_tmt(1),4)/100;
xPRESS_tmt = sum(rssPRESS_tmt(:,6) .* rssPRESS_tmt(:,4))/sum(rssPRESS_tmt(:,6))/100;
  plot(rssPRESS_tmt(:,6),rssPRESS_tmt(:,4)/100);
%ret

xPRESS_tlt = find(rssPRESS_tlt(:,6) == max(rssPRESS_tlt(:,6))); xPRESS_tlt = rssPRESS_tlt(xPRESS_tlt(1),4)/100;
xPRESS_tlt = sum(rssPRESS_tlt(:,6) .* rssPRESS_tlt(:,4))/sum(rssPRESS_tlt(:,6))/100;
  plot(rssPRESS_tlt(:,6),rssPRESS_tlt(:,4)/100);
  xPRESS_tlt = 760;
%ret

xPRESS_ttt = 0.5*(xPRESS_tmt + xPRESS_tls);

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsA = plevs(1:100)-plevs(2:101);
playsB = log(plevs(1:100)./plevs(2:101));
plays = playsA./playsB;
plays = plays(4:100);
plays = flipud(plays);

load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era3_lays_spanday01_save_lat_Sep7_2011.mat
tropics = find(abs(save_lat) <= 30);

figure(6); clf
  subplot(121)
  shadedErrorBarY(nanmean(water(tropics,:)),plays,nanstd(water(tropics,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(waterrate(tropics,:)),plays,nanstd(waterrate(tropics,:)),'rx-',1);
  hold off; title('AIRS(b) ERA(r) Water frac/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

  subplot(122)
  shadedErrorBarY(nanmean(temp(tropics,:)),plays,nanstd(temp(tropics,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(ptemprate(tropics,:)),plays,nanstd(ptemprate(tropics,:)),'rx-',1);

  oink = find(abs(nTLS.lats) < 30);
  
  hl=errorbar_x(nanmean(nTMT.trend(oink)),xPRESS_tmt,nanstd(nTMT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTUT.trend(oink)),xPRESS_ttt,nanstd(nTUT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTLS.trend(oink)),xPRESS_tls,nanstd(nTLS.trend(oink)),'ko'); set(hl,'linewidth',2)

  %these were 30 yr trends
  %oink = find(abs(rssTLT(:,1)) < 30);   errorbar_x(nanmean(rssTLT(oink,2))/10,xPRESS_tlt,nanstd(rssTLT(oink,2)),'k*')
  %oink = find(abs(rssTMT(:,1)) < 30);   errorbar_x(nanmean(rssTMT(oink,2))/10,xPRESS_tmt,nanstd(rssTMT(oink,2)),'k*')
  %oink = find(abs(rssTTT(:,1)) < 30);   errorbar_x(nanmean(rssTTT(oink,2))/10,xPRESS_ttt,nanstd(rssTTT(oink,2)),'k*')
  %oink = find(abs(rssTLS(:,1)) < 30);   errorbar_x(nanmean(rssTLS(oink,2))/10,xPRESS_tls,nanstd(rssTLS(oink,2)),'k*')
  %oink = find(abs(rssTTS(:,1)) < 30);   errorbar_x(nanmean(rssTTS(oink,2))/10,xPRESS_tts,nanstd(rssTTS(oink,2)),'k*')

  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(1),xPRESS_tls,rssTRP.amsu_rss_error_robust(1),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(2),xPRESS_tts,rssTRP.amsu_rss_error_robust(2),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(3),xPRESS_tmt,rssTRP.amsu_rss_error_robust(3),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(4),xPRESS_tlt,rssTRP.amsu_rss_error_robust(4),'g*'); set(hl,'linewidth',2)

  hold off; title('AIRS(b) ERA(r) Temp K/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);

figure(4); clf
  pcolor(save_lat,plays,double(waterrate')); shading interp; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA WV frac/yr ')
figure(5); clf
  pcolor(save_lat,plays,double(ptemprate')); shading interp; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA T K/yr ')

figure(1); caxis([-0.02 +0.02]); colorbar
figure(4); caxis([-0.02 +0.02]); colorbar

figure(2); caxis([-0.1 +0.1]); colorbar  
figure(5); caxis([-0.1 +0.1]); colorbar

figure(8);
  shadedErrorBar(save_lat,params(:,6),params_sigs(:,6),'bo-',1);
  hold on
  shadedErrorBar(save_lat,stemprate,stempratestd,'ro-',1);
  hold off
grid
title('AIRS(b) ERA(r) STemp K/yr'); %set(hl,'fontsize',10); grid

%figure(4)
%hold off
%  set(gca,'ydir','reverse'); grid; title('T K/yr')
%  hl = legend('OEM','ERA','location','east'); set(hl,'fontsize',10)

figure(9); clf
  shadedErrorBarY(nanmean(temp(tropics,:)),log10(plays),nanstd(temp(tropics,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(ptemprate(tropics,:)),log10(plays),nanstd(ptemprate(tropics,:)),'rx-',1);

  oink = find(abs(nTLS.lats) < 30);
  
  lxPRESS_tmt = log10(xPRESS_tmt);
  lxPRESS_ttt = log10(xPRESS_ttt);
  lxPRESS_tlt = log10(xPRESS_tlt);
  lxPRESS_tts = log10(xPRESS_tts);
  lxPRESS_tls = log10(xPRESS_tls);

  hl=errorbar_x(nanmean(nTMT.trend(oink)),lxPRESS_tmt,nanstd(nTMT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTUT.trend(oink)),lxPRESS_ttt,nanstd(nTUT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTLS.trend(oink)),lxPRESS_tls,nanstd(nTLS.trend(oink)),'ko'); set(hl,'linewidth',2)

  %these were 30 yr trends
  %oink = find(abs(rssTLT(:,1)) < 30);   errorbar_x(nanmean(rssTLT(oink,2))/10,lxPRESS_tlt,nanstd(rssTLT(oink,2)),'k*')
  %oink = find(abs(rssTMT(:,1)) < 30);   errorbar_x(nanmean(rssTMT(oink,2))/10,lxPRESS_tmt,nanstd(rssTMT(oink,2)),'k*')
  %oink = find(abs(rssTTT(:,1)) < 30);   errorbar_x(nanmean(rssTTT(oink,2))/10,lxPRESS_ttt,nanstd(rssTTT(oink,2)),'k*')
  %oink = find(abs(rssTLS(:,1)) < 30);   errorbar_x(nanmean(rssTLS(oink,2))/10,lxPRESS_tls,nanstd(rssTLS(oink,2)),'k*')
  %oink = find(abs(rssTTS(:,1)) < 30);   errorbar_x(nanmean(rssTTS(oink,2))/10,lxPRESS_tts,nanstd(rssTTS(oink,2)),'k*')

  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(1),lxPRESS_tls,rssTRP.amsu_rss_error_robust(1),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(2),lxPRESS_tts,rssTRP.amsu_rss_error_robust(2),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(3),lxPRESS_tmt,rssTRP.amsu_rss_error_robust(3),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(4),lxPRESS_tlt,rssTRP.amsu_rss_error_robust(4),'g*'); set(hl,'linewidth',2)

  hold off; title('AIRS(b) ERA(r) Temp K/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 1 3]);


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
