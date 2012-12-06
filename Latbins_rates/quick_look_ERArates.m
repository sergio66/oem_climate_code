figure(6); 

ecmrates = input('enter (1) 2002 - 2012 or (2) 2002-2011 or (3) 2007-2011  or (4) 2007-2012  ERA rates : ');
ecmfile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/';
if ecmrates == 1
  co2str = '_07_2002_07_2010';
  ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012.mat'];
elseif ecmrates == 2
  co2str = '_07_2002_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era2_lays_profilerates_May23_2011_robust.mat'];
elseif ecmrates == 3
  %ecmfile = [ecmfile ...
  %    'overocean_gsx_1dayV1_ecmwf2_lays_spanday16_profilerates_Aug18_2011_robust.mat'];
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era_lays_spanday01_profilerates_Aug26_2011_robust.mat'];
elseif ecmrates == 4
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Aug28_2012_robust_span_07_2007_07_2012.mat'];
end
load(ecmfile);

%disp('tracegas times : (1) 8yr 07/2002-07/2010 (2) 8 yr 09/2003-09/2011 (3) 9 yr 09/2002-09/2011');
%disp('tracegas times : (1) 8yr 09/2002-08/2010 (2) 7 yr 09/2003-08/2010')
%tracegasrates = input('Enter tracegas rates : ');
tracegasrates = 1;
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
iWhich = input('Enter (-4)SH (-3)SP (-2)SML (-1)ST (0)Tropics (+1)NT (+2)NML (3)NP (4)NH (5)All (10)individual bin (20)my subset : ');
if iWhich == 0
  lat_index = find(abs(save_lat) <= 30);
elseif iWhich == -4
  lat_index = find(save_lat <= 0);
elseif iWhich == +4
  lat_index = find(save_lat > 0);
elseif iWhich == -3
  lat_index = find(save_lat <= -60);
elseif iWhich == +3
  lat_index = find(save_lat >= 60);
elseif iWhich == -2
  lat_index = find(save_lat > -60 & save_lat < -30);
elseif iWhich == +2
  lat_index = find(save_lat > +30 & save_lat < +60);
elseif iWhich == -1
  lat_index = find(save_lat > -30 & save_lat < 0);
elseif iWhich == +1
  lat_index = find(save_lat > 00 & save_lat < 30);
elseif iWhich == +5
  lat_index = 1:36;
elseif iWhich == 10
  figure(9); clf; plot(save_lat,'o-'); grid
  lat_index = input('Enter individual bin : ');
elseif iWhich == 20
  figure(9); clf; plot(save_lat,'o-'); grid
  ijunk = input('Enter [latbinstart latbinstop] : ');
  lat_index = find(save_lat >= min(ijunk) & save_lat <= max(ijunk))
end

oink = find(abs(nTLS.lats) < 30);
oink = lat_index;
  
if length(lat_index) > 1
  nanmean_era_plots
else
  individual_era_plots
end
