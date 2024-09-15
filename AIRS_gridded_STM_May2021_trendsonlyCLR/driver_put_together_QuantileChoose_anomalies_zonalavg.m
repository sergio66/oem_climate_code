addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

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

%iHuge = -1;
%iHuge = +1;
%iHuge = input('Enter (-1/default) save one channel, all tiles, all timesteps or (+1) save all channels, all tiles, all timesteps (huge memory) or (+2) save f < 1650 cm-1, all timesteps (+3) 450 good channels : ');
%if length(iHuge) == 0
%  iHuge = -1;
%end
iHuge = +1; %% all channels

do_XX_YY_from_X_Y
coslat = cos(reshape(YY,4608,1)*pi/180)*ones(1,iNumAnomTimeSteps);

%if iHuge > 0
%  iDorA = input('+1 Desc/default or -1 Asc : ');
%  if length(iDorA) == 0
%    iDorA = +1;
%  end
%end 
iDorA = 0; %% both day and night

iaFound = zeros(size(iaFound));

if iHuge == +1
  disp('wow : lotsa lotsa memory! making the btanom and fluxanom ....')
  if iDorA > 0
    btanomD   = zeros(2645,iNumAnomTimeSteps,64); %% too much memory!
    fluxanomD = zeros(64,iNumAnomTimeSteps,length(RRTM_bands));
  elseif iDorA < 0
    btanomA   = zeros(2645,iNumAnomTimeSteps,64); %% too much memory!
    fluxanomA = zeros(64,iNumAnomTimeSteps,length(RRTM_bands));
  else
    btanomD   = zeros(2645,iNumAnomTimeSteps,64); %% too much memory!
    fluxanomD = zeros(64,iNumAnomTimeSteps,length(RRTM_bands));
    btanomA   = zeros(2645,iNumAnomTimeSteps,64); %% too much memory!
    fluxanomA = zeros(64,iNumAnomTimeSteps,length(RRTM_bands));
  end
  monitor_memory_whos
  ind = 1 : 2645;
end

bt1231_D = zeros(1,64);
bt1231_A = zeros(1,64);
btChID_D = zeros(1,64);
btChID_A = zeros(1,64);

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

iChID = iChIX;
iChIA = iChIX;

%% fname = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin64/LonBin72/iQAX_3_fits_LonBin72_LatBin64_V1_200200090001_202200080031_Anomaly_TimeStepsX457.mat'
%% Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_anomalies_quantiles.m
%%   which loads in   fn_summary = sprintf('../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin%1$02d/LonBin%2$02d/iQAX_3_summarystats_LatBin%1$02d_LonBin%2$02d_timesetps_001_%3$03d_V1.mat',lati,loni,i16daysSteps);

  f2645 = instr_chans2645;
  f2645 = h.vchan;

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
  figure(1); clf; plot(bt1231_D'); title('BT 1231'); 
  figure(2); clf; plot(btChID_D'); title('BT ChID'); 
  datime = 2002.75 + (1:iNumAnomTimeSteps)*16/365;

  if iDorA == +1
    figure(3); clf; junkjunk = squeeze(btanomD(i1231,:,:))';  pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom D');
    figure(4); clf; junkjunk = squeeze(btanomD(iXCh,:,:))';   pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom D');
  elseif iDorA == -1 
    figure(3); clf; junkjunk = squeeze(btanomA(i1231,:,:))';  pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom A');
    figure(4); clf; junkjunk = squeeze(btanomA(iXCh,:,:))';   pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom A');
  else
    figure(3); clf; junkjunk = squeeze(btanomD(i1231,:,:))';  pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom D');
    figure(4); clf; junkjunk = squeeze(btanomD(iXCh,:,:))';   pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom D');
    figure(5); clf; junkjunk = squeeze(btanomA(i1231,:,:))';  pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT 1231 anom A');
    figure(6); clf; junkjunk = squeeze(btanomA(iXCh,:,:))';   pcolor(datime,1:64,junkjunk); shading interp; colormap(usa2); colorbar; caxis([-1 +1]*5); title('BT ChID anom A');
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
      junkradA = squeeze(a.rad_anom_asc(iQuant,:,:));   
      junkradD = squeeze(a.rad_anom_desc(iQuant,:,:));
      if iHuge == +1 | iHuge == +2 | iHuge == +3
        if iDorA > 0
          btanomD(:,:,jj) = squeeze(junkD(iQuant,ind,:))/72 + btanomD(:,:,jj);
          fluxanomD(jj,:,1) = trapz(h.vchan,junkradD)/1000/72 + fluxanomD(jj,:,1);  %% total
          for flfl = 1 : length(RRTM_bands)-1
            junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
            fluxanomD(jj,:,flfl+1) = trapz(h.vchan(junk),junkradD(junk,:))/1000/72 + fluxanomD(jj,:,flfl+1);  %% rrtm components            
          end
        elseif iDorA < 0
          btanomA(:,:,jj) = squeeze(junkA(iQuant,ind,:))/72 + btanomA(:,:,jj);
          fluxanomA(jj,:,1) = trapz(h.vchan,junkradA)/1000/72 + fluxanomA(jj,:,1);  %% total
          for flfl = 1 : length(RRTM_bands)-1
            junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
            fluxanomA(jj,:,flfl+1) = trapz(h.vchan(junk),junkradA(junk,:))/1000/72 + fluxanomA(jj,:,flfl+1);  %% rrtm components
          end
        elseif iDorA == 0
          btanomA(:,:,jj) = squeeze(junkA(iQuant,ind,:))/72 + btanomA(:,:,jj);
          fluxanomA(jj,:,1) = trapz(h.vchan,junkradA)/1000/72 + fluxanomA(jj,:,1);  %% total
          btanomD(:,:,jj) = squeeze(junkD(iQuant,ind,:))/72 + btanomD(:,:,jj);
          fluxanomD(jj,:,1) = trapz(h.vchan,junkradD)/1000/72 + fluxanomD(jj,:,1);  %% total
          for flfl = 1 : length(RRTM_bands)-1
            junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
            fluxanomD(jj,:,flfl+1) = trapz(h.vchan(junk),junkradD(junk,:))/1000/72 + fluxanomD(jj,:,flfl+1);  %% rrtm components
            fluxanomA(jj,:,flfl+1) = trapz(h.vchan(junk),junkradA(junk,:))/1000/72 + fluxanomA(jj,:,flfl+1);  %% rrtm components
          end
          %% plot(1:500,squeeze(sum(fluxanomA(1,:,2:14),3)),'b.-',1:500,squeeze(fluxanomA(1,:,1)),'r',1:500,squeeze(fluxanomA(1,:,2:14)),'k')
          %% plot(RRTM_bands,squeeze(fluxanomD(1,500,1:14)),'ro-',RRTM_bands,squeeze(fluxanomA(1,500,1:14)),'bx-'); plotaxis2;  %% remember, the first one is the SUM!!!!!
          %% semilogy(RRTM_bands,abs(squeeze(fluxanomD(1,500,1:14))),'ro-',RRTM_bands,abs(squeeze(fluxanomA(1,500,1:14))),'bx-'); plotaxis2;
        end
      end
      bt1231_D(jj) = a.bt_desc(i1231,iQuant)/72 + bt1231_D(jj);
      btChID_D(jj) = a.bt_desc(iXCh,iQuant)/72  + btChID_D(jj);
      bt1231_A(jj) = a.bt_asc(i1231,iQuant)/72 + bt1231_A(jj);
      btChID_A(jj) = a.bt_asc(iXCh,iQuant)/72  + btChID_A(jj);
    end  %% iaFound
  end    %% loop ii
end      %% loop jj
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDorA > 0
  btanomD = permute(btanomD,[3 1 2]);
  put_together_QuantileChoose_anomalies_plots_AND_save_D_zonalavg
elseif iDorA > 0
  btanomA = permute(btanomA,[3 1 2]);
  put_together_QuantileChoose_anomalies_plots_AND_save_A_zonalavg
else
  btanomD = permute(btanomD,[3 1 2]);
  btanomA = permute(btanomA,[3 1 2]);
  put_together_QuantileChoose_anomalies_plots_AND_save_D_zonalavg
  put_together_QuantileChoose_anomalies_plots_AND_save_A_zonalavg
end

disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  
disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  
disp('just load in the huge mat file with anomalies, and rerun make_globalavg_and_N_average_anomalies.m');  

yymm = yyA + (mmA-1)/12 + (ddA-1)/30/12;
figure(1); plot(yymm,squeeze(nanmean(fluxanomD,1)))
figure(1); plot(yymm,squeeze(nanmean(fluxanomD(:,:,1),1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_fluxband_anomaly_timeseries

error('gslkjglkgjlksj')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iHuge > 0
  iX = 1;
  if iX > 0
    %% iQuant = 3;
    stand_alone_make_globalavg_and_N_average_anomalies_zonalavg
  end
end

