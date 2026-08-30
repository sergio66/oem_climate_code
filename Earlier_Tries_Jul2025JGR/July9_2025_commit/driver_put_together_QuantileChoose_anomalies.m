addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

%% see also /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
%% https://www.arm.gov/publications/proceedings/conf05/extended_abs/mlawer_ej.pdf
RRTM_bands0 = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
wx = 10:1:3000; rx = ttorad(wx,300); 
fluxx = trapz(wx,rx)/1000;
for ii = 2 : length(RRTM_bands0)
  junk = find(wx > RRTM_bands0(ii-1) & wx <= RRTM_bands0(ii));
  fluxx(ii) = trapz(wx(junk),rx(junk))/1000;
end

RRTM_bands = RRTM_bands0(4:end);

dirCode = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/';
dirData = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/driver_bt1231_timeseries_YY0MM0DD0_YYEMMEDDE_0001_4608.m')
disp('see Mitchell_Goldberg-Dissertation.pdf for list of chans : wget https://aosc.umd.edu/sites/default/files/dissertations-theses/Mitchell%20Goldberg-Dissertation.pdf')

  disp('Fig 5.8        667.766 cm-1       1 mb')
  disp('Fig 5.9        667.715 cm-1       2 mb')
  disp('Fig 5.9        667.270 cm-1      15 mb')
  disp('Fig 5.8        667.018 cm-1      25 mb')
  disp('Fig 4.4        681.457 cm-1      90 mb')
  disp('Fig 4.5        704.436 cm-1     350 mb')
  disp('Fig 4.6        723.029 cm-1     700 mb')
  disp('Fig 4.7        801.099 cm-1     850 mb')
  disp(' ')
  disp('Fig 4.8        1519.07 cm-1     315 mb')
  disp('Fig 4.9        1598.49 cm-1     490 mb  MidTropWV')
  disp('Fig 5.10       1519.07 cm-1     315 mb  UTWV')
  disp('Fig 5.1        1520.87 cm-1     315 mb???')
  disp(' ')
  disp('Fig 5.4        1040.03 cm-1     80 mb???')
  disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_anomalies_quantiles.m')
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_anomalies_quantiles.m')
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_anomalies_quantiles.m')

load h2645structure.mat

iNumYears = 20;    iNumAnomTimeSteps = 454; % 365/16 = 22.8125; floor(iNumYears*22.8125 - 2) = 454!!!! (about 23 timesteps/year, 2 are missing because of missing AIRS data)
iNumYears = 21.75; iNumAnomTimeSteps = 496; % 365/16 = 22.8125; floor(iNumYears*22.8125 - 2) = 498!!!! (about 23 timesteps/year, 2 are missing because of missing AIRS data)
iNumYears = 22;    iNumAnomTimeSteps = 500; % 365/16 = 22.8125; floor(iNumYears*22.8125 - 2) = 500!!!! (about 23 timesteps/year, 2 are missing because of missing AIRS data)

iQAX = input('Enyter iQAX [1,3,4] 3=default : ');
if length(iQAX) == 0
  iQAX = 3;
end

iCnt = 0;
for jj = 1 : 64
  fprintf(1,'latbin %2i of 64 \n',jj);
  for ii = 1 : 72
    iCnt = iCnt + 1;
    if iNumYears == 20
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'];
    elseif iNumYears == 21.75
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_Anomaly_TimeSteps498.mat'];
    elseif iNumYears == 22
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_Anomaly_TimeSteps502.mat'];
    else
      error('unknown iNumYears')
    end
    iaFound(iCnt) = 0;   
    if exist(fname)
      iaFound(iCnt) = +1;
    end
  end
end
fprintf(1,'found %4i of 4608 anomaly files that look like %s \n',sum(iaFound),fname);
disp(' ')

%% Readme_Anomaly_TWP_lat35_lon66
iTWP = (35-1)*72 + 66;

iHuge = -1;
iHuge = +1;
iHuge = input('Enter (-1/default) save one channel, all tiles, all timesteps or (+1) save all channels, all tiles, all timesteps (huge memory) or (+2) save f < 1650 cm-1, all timesteps (+3) 450 good channels : ');
if length(iHuge) == 0
  iHuge = -1;
end

do_XX_YY_from_X_Y
coslat = cos(reshape(YY,4608,1)*pi/180)*ones(1,iNumAnomTimeSteps);

if iHuge > 0
  iDorA = input('+1 Desc/default or -1 Asc : ');
  if length(iDorA) == 0
    iDorA = +1;
  end
end 
iaFound = zeros(size(iaFound));

if iHuge == +1
  disp('wow : lotsa lotsa memory! making the btanom and fluxanom ....')
  if iDorA > 0
    btanomD = nan(4608,2645,iNumAnomTimeSteps); %% too much memory!
    fluxanomD = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66D = nan(2645,iNumAnomTimeSteps);
  else
    btanomA = nan(4608,2645,iNumAnomTimeSteps); %% too much memory!
    fluxanomA = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66A = nan(2645,iNumAnomTimeSteps);
  end
  monitor_memory_whos
  ind = 1 : 2645;
elseif iHuge == +2
  disp('wow : lotsa lotsa memory! making the btanom and fluxanom ....')
  f2645 = instr_chans2645;
  f2645 = h.vchan;
  ind = find(f2645 < 1650); [mean(diff(ind)) std(diff(ind))];
  if iDorA > 0
    btanomD = nan(4608,length(ind),iNumAnomTimeSteps);  %% lot of memory!
    fluxanomD = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66D = nan(2645,iNumAnomTimeSteps);
  else
    btanomA = nan(4608,length(ind),iNumAnomTimeSteps);  %% lot of memory!
    fluxanomA = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66A = nan(2645,iNumAnomTimeSteps);
  end
  monitor_memory_whos
elseif iHuge == +3
  disp('wow : lotsa memory! making the btanomD .... ')
  f2645 = instr_chans2645;
  f2645 = h.vchan;
  ind = get_climateQA_LW_noCFC(h.vchan);
  ichan = h.ichan(ind);

  [mean(diff(ind)) std(diff(ind))];
  if iDorA > 0
    btanomD = nan(4608,length(ind),iNumAnomTimeSteps);  %% lot of memory!
    fluxanomD = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66D = nan(2645,iNumAnomTimeSteps);
  else
    btanomA = nan(4608,length(ind),iNumAnomTimeSteps);  %% lot of memory!
    fluxanomA = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
    btanomTWP_35_66A = nan(2645,iNumAnomTimeSteps);
  end
  monitor_memory_whos
else
  btanomD = nan(4608,iNumAnomTimeSteps);
  btanomA = nan(4608,iNumAnomTimeSteps);
  fluxanomD = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
  fluxanomA = nan(4608,iNumAnomTimeSteps,length(RRTM_bands));
  btanomTWP_35_66D = nan(2645,iNumAnomTimeSteps);
  btanomTWP_35_66A = nan(2645,iNumAnomTimeSteps);
end

bt1231_D = nan(1,4608);
bt1231_A = nan(1,4608);
btChID_D = nan(1,4608);
btChID_A = nan(1,4608);

disp('reading in anomalies')

iCnt = 0;
iQuant = 3;
iQuant = 5;
iQuant = input('Enter quantile : 1 .. 5 (03 is deafult) : ');
if length(iQuant) == 0
  iQuant = 3;
end

iChIX  = 0213;  %% 0704 cm-1     find(h.vchan >= 0704,1)
iChIX  = 1145;  %% 1040 cm-1     find(h.vchan >= 1040,1)
iChIX  = 2025;  %% 1519 cm-1     find(h.vchan >= 1519,1)
iChIX  = 1520;  %% 1231 cm-1     find(h.vchan >= 1231,1)
iChIX  = 1511;  %% 1226.82 cm-1  find(h.vchan >= 1226.65,1)
iChIX  = 1861;  %% 1419 cm-1     find(h.vchan >= 1419,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
[h,ha,p,pa] = rtpread('/asl/rtp/airs/airs_l1c_v674/clear/2023/ecmwf_airicrad_day022_clear.rtp');

% ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/tile_fits_quantiles.m
% for qi = qi1 : qi2
%    % r1231 = squeeze(d.rad_quantile_desc(:,1520,qi));
%    % r1228 = squeeze(d.rad_quantile_desc(:,1513,qi));
%    r1231 = squeeze(d.rad_quantile_desc(k_desc,1520,qi));
%    r1228 = squeeze(d.rad_quantile_desc(k_desc,1513,qi));
%    bt1231 = rad2bt(fairs(1520),r1231);
%    bt1228 = rad2bt(fairs(1513),r1228);
% end

boo = [1511 1520];

addpath /home/sergio/MATLABCODE
load /home/sergio/KCARTA/L2SComparisons/l2s_kc122_H16_605_2830.mat
[fc,qc] = convolve_airs(double(w),double(d),double(h.ichan));
semilogy(w,d(:,1),fc,qc(:,1),fc(boo),qc(boo,1),'ro'); xlim([1220 1240])

plot(p.rlon,p.rlat,'.')
lala = find(p.rlat >= 0,1);

robs = nanmean(p.robs1,2);
robs = p.robs1(:,lala);

plot(h.vchan,rad2bt(h.vchan,robs),h.vchan(boo),rad2bt(h.vchan(boo),robs(boo)),'ro'); xlim([1220 1235])

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iHuge < 0
  junk = -1;
  while junk < 0
    junk = input('Enter wavenumber to look for : ');
    iChIX = find(h.vchan >= junk,1);
    fprintf(1,'iChIX = %4i corresponds to SARTA ID %4i %6.2f cm-1 \n',iChIX,h.ichan(iChIX),h.vchan(iChIX))
    junk = input('Proceed (-1/+1[default]) : ');
    if length(junk) == 0
      junk = +1;
    end
    if junk < 0
      return
    end
  end
end

iChID = iChIX;
iChIA = iChIX;

%% fname = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin64/LonBin72/iQAX_3_fits_LonBin72_LatBin64_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'
%% Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_anomalies_quantiles.m
%%   which loads in   fn_summary = sprintf('../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin%1$02d/LonBin%2$02d/iQAX_3_summarystats_LatBin%1$02d_LonBin%2$02d_timesetps_001_%3$03d_V1.mat',lati,loni,i16daysSteps);

fprintf(1,'iChIX = 1520 and %04i  --> %8.3f    %8.3f cm-1 \n',iChIX,f2645(1520),f2645(iChIX))

%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
i1231 = find(ind == 1520);
iXCh  = find(ind == iChIX);
f_ind = f2645(ind);
fprintf(1,'f2645(1520)  --> h.vchan(i1231)     : %8.6f %8.6f \n',f2645(1520),h.vchan(ind(i1231)))
fprintf(1,'f2645(iChIX) --> h.vchan(ind(iXCh)) : %8.6f %8.6f \n',f2645(iChIX),h.vchan(ind(iXCh)))
%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
%% this is ireelevant if using ALL channels of f < 620 cm-1 ... but is needed in case using only 480 great channels;
 
for jj = 1 : 64
  fprintf(1,' \n latbin %2i of 64  will loop over 72 lonbins .... \n',jj);
  figure(1); clf; pcolor(reshape(bt1231_D,72,64)'); title('BT 1231'); shading interp; colorbar; xlabel('LonBin 01-72'); ylabel('LatBin 01-64'); colormap jet; pause(0.1);
  figure(2); clf; pcolor(reshape(btChID_D,72,64)'); title('BT ChID'); shading interp; colorbar; xlabel('LonBin 01-72'); ylabel('LatBin 01-64'); colormap jet; pause(0.1);
  datime = 2002.75 + (1:iNumAnomTimeSteps)*16/365;

  if iDorA == +1
    figure(3); clf; junkjunk = squeeze(btanomD(:,i1231,:));  junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); 
               junkjunk = squeeze(nanmean(junkjunk,1)); pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom D');
    figure(4); clf; junkjunk = squeeze(btanomD(:,iXCh,:)); junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); 
               junkjunk = squeeze(nanmean(junkjunk,1)); pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom D');
  else
    figure(5); clf; junkjunk = squeeze(btanomA(:,i1231,:));  junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); 
               junkjunk = squeeze(nanmean(junkjunk,1)); pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom A');
    figure(6); clf; junkjunk = squeeze(btanomA(:,iXCh,:)); junkjunk = reshape(junkjunk,72,64,iNumAnomTimeSteps); 
               junkjunk = squeeze(nanmean(junkjunk,1)); pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom A');
  end
  pause(0.1) 
  for ii = 1 : 72
    if mod(ii,10) == 0 
      fprintf(1,'+')
    else
      fprintf(1,'.')
    end
    iCnt = iCnt + 1;
    if iNumYears == 20
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'];
    elseif iNumYears == 21.75
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_Anomaly_TimeSteps498.mat'];
    elseif iNumYears == 22
      fname = ['LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/iQAX_' num2str(iQAX) '_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d')];
      fname = [dirData fname '_V1_Anomaly_TimeSteps502.mat'];
    else
      error('unknown iNumYears')
    end
    iaFound(iCnt) = 0;   
    if exist(fname)
      iaFound(iCnt) = +1;
      a = load(fname,'bt_anom_desc','bt_desc','bt_anom_asc','bt_asc','rad_anom_asc','rad_anom_desc');
      junkD = a.bt_anom_desc;
      junkA = a.bt_anom_asc;
      junkradA = squeeze(a.rad_anom_asc(iQuant,:,:));   %BUG THIS IS ZERO
      junkradD = squeeze(a.rad_anom_desc(iQuant,:,:));
      if iHuge == +1 | iHuge == +2 | iHuge == +3
        if iDorA > 0
          btanomD(iCnt,:,:) = squeeze(junkD(iQuant,ind,:));
          fluxanomD(iCnt,:,1) = trapz(h.vchan,junkradD)/1000;
          if iCnt == iTWP
            btanomTWP_35_66D = squeeze(junkD(iQuant,:,:));
          end
          for flfl = 1 : length(RRTM_bands)-1
            junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
            fluxanomD(iCnt,:,flfl+1) = trapz(h.vchan(junk),junkradD(junk,:))/1000;
          end
        else
          btanomA(iCnt,:,:) = squeeze(junkA(iQuant,ind,:));
          fluxanomA(iCnt,:,1) = trapz(h.vchan,junkradA)/1000;
          if iCnt == iTWP
            btanomTWP_35_66A = squeeze(junkA(iQuant,:,:));
          end
          for flfl = 1 : length(RRTM_bands)-1
            junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
            fluxanomA(iCnt,:,flfl+1) = trapz(h.vchan(junk),junkradA(junk,:))/1000;
          end
        end
      else
        btanomD(iCnt,:) = squeeze(junkD(iQuant,iXCh,:));
        btanomA(iCnt,:) = squeeze(junkA(iQuant,iXCh,:));
        if iCnt == iTWP
          btanomTWP_35_66A = squeeze(junkA(iQuant,:,:));
          btanomTWP_35_66D = squeeze(junkD(iQuant,:,:));
        end
        fluxanomD(iCnt,:,1) = trapz(h.vchan,junkradD)/1000;
        fluxanomA(iCnt,:,1) = trapz(h.vchan,junkradA)/1000;
        for flfl = 1 : length(RRTM_bands)-1
          junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
          fluxanomD(iCnt,:,flfl+1) = trapz(h.vchan(junk),junkradD(junk,:))/1000;
          fluxanomA(iCnt,:,flfl+1) = trapz(h.vchan(junk),junkradA(junk,:))/1000;
        end
      end
      bt1231_D(iCnt) = a.bt_desc(i1231,iQuant);
      btChID_D(iCnt) = a.bt_desc(iXCh,iQuant);
      bt1231_A(iCnt) = a.bt_asc(i1231,iQuant);
      btChID_A(iCnt) = a.bt_asc(iXCh,iQuant);
    end  %% iaFound
  end    %% loop ii
end      %% loop jj

%{
a = load(fname,'yy_desc','mm_desc','dd_desc','hh_desc','rtime_desc');
yyD = a.yy_desc;
mmD = a.mm_desc;
ddD = a.dd_desc;
hhD = a.hh_desc;
rtimeD = a.rtime_desc;

a = load(fname,'yy_asc','mm_asc','dd_asc','hh_asc','rtime_asc');
yyA = a.yy_asc;
mmA = a.mm_asc;
ddA = a.dd_asc;
hhA = a.hh_asc;
rtimeA = a.rtime_asc;

comment = 'see driver_put_together_QuantileChoose_anomalies.m';
if iDorA > 0
  saver = ['save btanomTWP_35_66_Q' num2str(iQuant,'%02d') '_numyears_' num2str(ceil(iNumYears)) '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat btanomTWP_35_66D yyD mmD ddD hhD rtimeD comment'];
elseif iDorA < 0
  saver = ['save btanomTWP_35_66_Q' num2str(iQuant,'%02d') '_numyears_' num2str(ceil(iNumYears)) '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat btanomTWP_35_66A yyA mmA ddA hhA rtimeA comment'];
else
  saver = ['save btanomTWP_35_66_Q' num2str(iQuant,'%02d') '_numyears_' num2str(ceil(iNumYears)) '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat btanomTWP_35_66A btanomTWP_35_66D yyD mmD ddD hhD rtimeD comment'];
end
ix1519 = find(f2645 >= 1519,1)
ix1558 = find(f2645 >= 1558,1)
ix1419 = find(f2645 >= 1419,1)
ix1231 = find(f2645 >= 1231,1)
yymmjunk = yyD + (mmD-1)/12 + (ddD-1)/12/30;
plot(yymmjunk,btanomTWP_35_66D([1520 1861],:)); legend('1231','1419');
plot(yymmjunk,smooth(btanomTWP_35_66D(1520,:),23),yymmjunk,smooth(btanomTWP_35_66D(1861,:),23),yymmjunk,smooth(btanomTWP_35_66D(2025,:),23)); plotaxis2; legend('1231','1419','1519');
plot(yymmjunk,smooth(btanomTWP_35_66D(1520,:),5),yymmjunk,smooth(btanomTWP_35_66D(1861,:),5),yymmjunk,smooth(btanomTWP_35_66D(2084,:),5),'linewidth',2); plotaxis2; legend('1231','1419','1558');

mamoo = load('anomaly_tile_2515_timeseries_Q03.mat')
mamoo.yymm = mamoo.yy + (mamoo.mm-1)/12 + (mamoo.dd-1)/12/30;
plot(mamoo.yymm,smooth(mamoo.btavgAnomFinal(1520,:),5),mamoo.yymm,smooth(mamoo.btavgAnomFinal(1861,:),5),mamoo.yymm,smooth(mamoo.btavgAnomFinal(2084,:),5),'linewidth',2); plotaxis2; legend('1231','1419','1558');

%% this is HUGE
%  mazoo = load('/asl/s1/sergio/JUNK/anomaly_ALL_Q03.mat');
% sum(sum(squeeze(mazoo.btanom(iTWP-1,:,:))-mamoo.btavgAnomFinal))   %% iTWP = 2515, play with 2514 or 2516
% sum(sum(squeeze(mazoo.btanom(iTWP,:,:))-mamoo.btavgAnomFinal))
% sum(sum(squeeze(mazoo.btanom(iTWP-1,:,:))-mamoo.btavgAnomFinal))
% sum(sum(squeeze(mazoo.btanom(iTWP+1,:,:))-mamoo.btavgAnomFinal)) YAY

wazoo = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/wgtfcn_jac.mat');
pcolor(wazoo.fout,1:97,wazoo.jout'); shading interp; xlim([640 1640])
for pp = 1 : 2378
  moo = find(
eval(saver)
%}
error('save btanomTWP_35_66A')

fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDorA > 0
  put_together_QuantileChoose_anomalies_plots_AND_save_D
else
  put_together_QuantileChoose_anomalies_plots_AND_save_A
end

disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  
disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  
disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iHuge > 0
  %iX = input('Do you wish to do global and (+1/default) N zonal averages or (-1) one tile : ');
  iX = 1;
  if length(iX) == 0
    iX = 1;
  end
  if iX > 0
    stand_alone_make_globalavg_and_N_average_anomalies
  else
    make_globalavg_and_onetile_anomaly
    %make_onetile_anomaly
  end
end

