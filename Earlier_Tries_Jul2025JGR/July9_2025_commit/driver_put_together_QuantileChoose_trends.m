%thedir0 = '/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv2/LatBin' num2str(JOB,'%02d') '/';
%iLon = 1;
%thefilein  = [thedir0 '/' LonBin' num2str(iLon,'%02d') '/cfbins_LatBin' num2str(JOB,'%02d') '_LonBin' num2str(iLon,'%02d') '_V2.mat'];

addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

% if ~exist('b_obs')
%   %load convert_strowrates2oemrates_allskygrid_obsonly.mat
%   load convert_sergio_clearskygrid_obsonly_Q16.mat
%   b_obs = b_desc;
% end

disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_fits_quantiles.m');
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_fits_quantiles.m');
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_fits_quantiles.m');

disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/ReadmeQuick to figure out how to make trends from YYYY1/MM1/DD1 to YYYY2/MM2/DD2')
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/ReadmeQuick to figure out how to make trends from YYYY1/MM1/DD1 to YYYY2/MM2/DD2')
disp('Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/ReadmeQuick to figure out how to make trends from YYYY1/MM1/DD1 to YYYY2/MM2/DD2')
disp(' ')

if ~exist('h')
  load h2645structure.mat
end

airs_noise = instr_chans2645('airs',2);

load latB64.mat
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
figure(1); scatter_coast(X(:),Y(:),50,landfrac); colorbar; title('landfrac');  caxis([0 1])

b_asc = nan(72,64,2645);
b_desc = nan(72,64,2645);

iType = +01;   %% strows Q16 18 year trends         Strow ran for me in Feb 2021, 2002/09 to 2020/08
iType = -01;   %% sergio Q16 18 year TEST trends    I     ran for me in Aug 2021, 2002/09 to 2020/08, to check I get same results as Strow +01
iType = +02;   %% sergio Q16 19 year trends         I     ran for me in Aug 2021, 2002/09 to 2021/07
iType = +03;   %% sergio Extreme 19 year trends     I     ran for me in Aug 2021, 2002/09 to 2020/08
iType = -03;   %% sergio Mean 19 year trends        I     ran for me in Aug 2021, 2002/09 to 2020/08
iType = +04;   %% sergio Q16 19 year trends         I     ran for me in Aug 2021, 2002/09 to 2021/08 FULL
iType = +05;   %% sergio Q16 12 year trends         I     ran for me in Aug 2022, 2002/09 to 2014/08 FULL
iType = +06;   %% sergio Q16 07 year trends         I     ran for me in Aug 2022, 2012/05 to 2019/04 FULL, same as Suomi NPP NSR
iType = +07;   %% sergio Q16 20 year trends         I     ran for me in Sep 2022, 2002/09 to 2022/08 FULL, 20 years
iType = +08;   %% sergio Q16 06 year trends         I     ran for me in Sep 2022, 2015/01 to 2021/12 FULL, 06 years, OCO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iType = +09;   %% sergio iQAX_3 Q05 20 year trends  I     ran for me in Oct 2022, 2002/09 to 2022/08 FULL, 20 years, with newer defn of quantiles
iType = +10;   %% sergio iQAX_3 Q05 05 year trends  I     ran for me in June 2023, 2002/09 to 2007/08 FULL, 05 years, with newer defn of quantiles
iType = +11;   %% sergio iQAX_3 Q05 10 year trends  I     ran for me in June 2023, 2002/09 to 2012/08 FULL, 10 years, with newer defn of quantiles
iType = +12;   %% sergio iQAX_3 Q05 15 year trends  I     ran for me in June 2023, 2002/09 to 2015/08 FULL, 15 years, with newer defn of quantiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iType = +13;   %% sergio iQAX_4 Q03 08 year trends  I     ran for me in Oct 2023, 2002/09 to 2010/08 08 years, with alternate defn of quantiles, beginning 8 years of mission, SW drifting
iType = +14;   %% sergio iQAX_3 Q03 04 years trends for AIRS RTA report 2018/09 - 2022/08
iType = +15;   %% sergio iQAX_3 Q03 14 years trends for Sarah/Cathy     2008/01 - 2022/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iType = +16;   %% sergio iQAX_3 Q03 last 4 years trends                 2020/07 - 2024/06
iType = +17;   %% sergio iQAX_3 Q03 21.0   years trends                 2003/07 - 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iType = +30;   %% sergio iQAX_1 Q01 20 years trends ASMU full average   2002/09 - 2022/08

disp('Choices DataSet to use ')
disp(' Look at  ../../oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_fits_quantiles.m, and set_iQAX');
disp(' datasets 1-8 are mostly interquantile')
disp('     iQAX = +1; %% quantile   quants = [0 0.01 0.02 0.03 0.04 0.05 0.10 0.25 0.50 0.75 0.9 0.95 0.96 0.97 0.98 0.99 1.00];')
disp(' datasets 9 is Q(x-->1.0)')
disp('     iQAX = +3; %% quantile   quants = [0 0.50 0.9 0.95 0.97 1.00];');

disp(' <---------------------------------------------------------------------------------------> ')
disp('                       (+1) Strow  Quantile Mar 2021 2002/09 to 2020/08 Full 18 years');
disp('                       (-1) Sergio Quantile Aug 2021 2002/09 to 2020/08 Full 18 years');
disp('                        (2) Sergio Quantile Aug 2021 2002/09 to 2021/07 ');
disp('                        (4) Sergio Quantile Aug 2021 2002/09 to 2021/08 Full 19 years **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                        (3) Sergio Extreme Aug 2021 2002/09 to 2021/07 ');
disp('                       (-3) Sergio Mean Aug 2021 2002/09 to 2021/07 ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                        (5) Sergio Quantile Aug 2022 2002/09 to 2014/08 Full 12+ years **** ');
disp('                        (6) Sergio Quantile Sep 2022 2012/05 to 2019/04 Suomi NPP years **** ');
disp('                        (7) Sergio Quantile Sep 2022 2002/09 to 2022/08 Full 20+ years **** ');
disp('                        (8) Sergio Quantile Sep 2022 2015/01 to 2021/12 Full 06+ OCO2 years **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                        (9) Sergio Quantile Sep 2022 2002/09 to 2022/08 Full 20+ years, new quantile defn **** ');
disp('                       (10) Sergio Quantile Jun 2023 2002/09 to 2007/08 Full 05+ years, new quantile defn **** ');
disp('                       (11) Sergio Quantile Jun 2023 2002/09 to 2012/08 Full 10+ years, new quantile defn **** ');
disp('                       (12) Sergio Quantile Jun 2023 2002/09 to 2017/08 Full 15+ years, new quantile defn **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                       (13) Sergio Quantile Sep 2023 2002/09 to 2010/08 First 08+ years, three quantiles defn **** ');
disp('                       (14) Sergio Quantile Jan 2024 2018/09 to 2022/08 last  04- years, new quantiles defn **** ');
disp('                       (15) Sergio Quantile Jan 2024 2008/01 to 2022/12 Mid   14 years, new quantiles defn **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                       (16) Sergio Quantile Aug 2024 2020/07 to 2024/06 last  04- years, new quantiles defn **** ');
disp('                       (17) Sergio Quantile Aug 2024 2003/07 to 2024/06 21 years,        new quantiles defn **** ');
disp(' <---------------------------------------------------------------------------------------> ')
disp('                       (30) Sergio Quantile Aug 2024 2002/09 to 2022/08 20 years AMSU **** ');
iType = input('Enter DataSet to use (+1,-1,+2,+4,+5,+6,+7,+8  or  +9,+10,+11,+12   or +3,-3   or +13,+14,+15,+16      or +30 AMSU) : ');

if iType <= 8
  quants =  [0.0100 0.0200 0.0300 0.0400 0.0500 0.1000 0.2500 0.5000 0.7500 0.9000 0.9500 0.9600 0.9700 0.9800 0.9900 1];
elseif iType >=9 & iType <= 12 | iType >= 14
  quants =  [0.5000 0.7000 0.9000 0.9500 0.9700 1];
elseif iType == 13
  quants =  [0.0000 0.0300 0.5000 0.9700 1];
elseif iType == 30
  quants =  [0.0000 1];
end

if iType == 30
  b_asc  = nan(72,64,14);
  b_desc = nan(72,64,14);
end

if iType ~= 3 & iType < 9
  iQuantile = 16;  %% hottest, used for AIRS STM May 21
  iQuantile = 08;
  iQuantile = 04;
  iQuantile = input('Enter iQuantile to make (1-16, 0 = avg, 50 = hottest 5) : ');
elseif iType >= 9 & iType <= 12 | iType >= 14 & iType <= 29
  iQuantile = 03;  %% Q0.95, used for AIRS STM May 22
  iQuantile = input('Enter iQuantile to make (1-5, 1 = mean (Q0.50) , 5 = clearest (Q0.99) : ');
elseif iType == 13
  iQuantile = 03;  %% Q0.97, used for AIRS STM Oct 2023
  iQuantile = input('Enter iQuantile to make (1-3, 1 = average, 2 = Q0.03, 3 = Q0.97 : ');
elseif iType == 30
  %iQuantile = input('Enter iQuantile to make (1 only');
  iQuantile = 01;  %% Q0.97, used for AIRS STM Oct 2023
end

iAllorSeasonal = 1;
if iQuantile == 09
  iAllorSeasonal = input('Enter ALL (+1/default) or (-1) DJF (-2) MAM (-3) JJA (-4) SON : ');
  if length(iAllorSeasonal) == 0
    iAllorSeasonal = 1;
  end
end

iKeepPlotting = -1;

if iType == 13
  iQAX = 3;
elseif iType >= 9 & iType <= 15
  iQAX = 3;
else
  iQAX = [];
end

if iType == 1
  fnamePROCESS = ['convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == -1
  fnamePROCESS = ['iType_-1_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 2
  fnamePROCESS = ['iType_2_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 4
  fnamePROCESS = ['iType_4_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 5
  fnamePROCESS = ['iType_5_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 6
  fnamePROCESS = ['iType_6_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 7
  fnamePROCESS = ['iType_7_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 8
  fnamePROCESS = ['iType_8_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
%%%%
elseif iType == 9
  fnamePROCESS = ['iType_10_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 10
  fnamePROCESS = ['iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 11
  fnamePROCESS = ['iType_11_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 12
  fnamePROCESS = ['iType_12_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
%%%%
elseif iType == 13
  fnamePROCESS = ['iType_13_iQAX_4_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 14
  fnamePROCESS = ['iType_14_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 15
  fnamePROCESS = ['iType_15_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 16
  fnamePROCESS = ['iType_16_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
elseif iType == 17
  fnamePROCESS = ['iType_17_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d')];
%%%%
elseif iType == 3
  fnamePROCESS = ['iType_3_extreme_convert_sergio_clearskygrid_obsonly.mat'];
elseif iType == -3
  fnamePROCESS = ['iType_-3_mean_convert_sergio_clearskygrid_obsonly.mat'];
%%%%%
elseif iType == 30
  fnamePROCESS = ['iType_30_AMSU_iQAX_' num2str(iQuantile,'%02d')];
end

if iAllorSeasonal == -1
  fnamePROCESS = [fnamePROCESS '_DJF'];
elseif iAllorSeasonal == -2
  fnamePROCESS = [fnamePROCESS '_MAM'];
elseif iAllorSeasonal == -3
  fnamePROCESS = [fnamePROCESS '_JJA'];
elseif iAllorSeasonal == -4
  fnamePROCESS = [fnamePROCESS '_SON'];
end

fnamePROCESS = [fnamePROCESS '.mat'];

if exist(fnamePROCESS)
 fnamePROCESS
 error('output file exists!')
end

saver = ['save ' fnamePROCESS ' b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];
saver = ['save ' fnamePROCESS ' b_* X Y landfrac salti h lagcor* mean_BT airs_noiseTtrue'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iType == 30
  disp('easy to do AMSU!!!')

    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0    = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/'];
    thedirERA5 = ['/NOTDONE//home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/'];

  thefilein = [thedir0 '/trends_anomalies_AMSU_20year_72x64.mat'];

  x = load(thefilein);
  for iLon = 1 : 72
    for iLat = 1 : 64
      %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%%  
      b_asc(iLon,iLat,:) = x.trend_T(:,iLon,iLat);
      b_desc(iLon,iLat,:) = x.trend_T(:,iLon,iLat);
      b_err_asc(iLon,iLat,:) = x.trend_T_err(:,iLon,iLat);
      b_err_desc(iLon,iLat,:) = x.trend_T_err(:,iLon,iLat);
      lagcor_obs_anom_asc(iLon,iLat,:)  = NaN;
      lagcor_obs_anom_desc(iLon,iLat,:) = NaN;
      mean_BT(iLon,iLat,:)  = NaN;
    end
  end

  figure(1); scatter_coast(X',Y',50,squeeze(b_desc(:,:,5))'); colorbar; colormap(usa2); caxis([-1 +1]*0.15); shading interp; title('d Ch5/dt');

  if ~exist(fnamePROCESS)
    airs_noiseTtrue = [];
    eval(saver);
    saver
  else 
    fprintf(1,'%s already exists \n',fnamePROCESS);
  end

  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iLat = 1 : 64
  if iType == +1
    thedir0 = ['/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv1_fits/LatBin' num2str(iLat,'%02d') '/'];   %% STROW STUFF March 2021
  elseif iType == -1
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 2
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 3
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == -3
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps429.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon_v3/Mean/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 4
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
  elseif iType == 5
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0 = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iType == 6 | iType == 7 | iType == 8 | iType == 9
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0    = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
    thedirERA5 = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/SimulateTimeSeries/ERA5_ConstTracegas/'];
  elseif iType >= 10 & iType <= 15
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0    = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
    thedirERA5 = ['/NOTDONE/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/SimulateTimeSeries/ERA5_ConstTracegas/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iType >= 16 & iType <= 29
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    thedir0    = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon/LatBin' num2str(iLat,'%02d') '/'];
    thedirERA5 = ['/NOTDONE/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/SimulateTimeSeries/ERA5_ConstTracegas/'];
  elseif iType == 30
    error('already done')
    % Look at /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/compare_bt1231trends_Q16_vs_extreme.m
    %fname = ['../DATAObsStats_StartSept2002_CORRECT_LatLon_v3/Extreme/LatBin' num2str(jj,'%02d') '/LonBin' num2str(ii,'%02d') '/extreme_fits_LonBin' num2str(ii,'%02d') '_LatBin' num2str(jj,'%02d') '_V1_TimeSteps433.mat'];
    thedir0    = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/'];
    thedirERA5 = ['/NOTDONE//home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/'];
  end

  for iLon = 1 : 72
    fprintf(1,'+')
    if iType == +1
      %% full 18 year by Strow
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1.mat'];
    elseif iType == -1
      %% full 18 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps412.mat'];
    elseif iType == 2
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == 3
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 'LonBin' num2str(iLon,'%02d') '/extreme_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == -3
      %% paritally incomplete 19 year by Sergio
      thefilein = [thedir0 'LonBin' num2str(iLon,'%02d') '/mean_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps429.mat'];
    elseif iType == 4
      %% full 19 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps433.mat'];
      thefileERA5 = [thedirERA5 '/reconstruct_era5_const_tracegas_spectra_geo_rlat' num2str(iLat,'%02d') '.mat'];
    elseif iType == 5
      %% full 19 year by Sergio
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200500010001_201400120031_TimeStepsX228.mat'];
    elseif iType == 6
      %% full Suomi CrIS NSR 7 year by Sergio
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin01/LonBin37/fits_LonBin37_LatBin01_V1_201200050001_201900040030_TimeStepsX159
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_201200050001_201900040030_TimeStepsX159.mat'];
    elseif iType == 7
      %% full AIRS 20 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457.mat'];
    elseif iType == 8
      %% full OCO2 06 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/fits_LonBin*_LatBin*_V1_201500010001_202100120031_TimeStepsX160.mat
      thefilein   = [thedir0 '/LonBin' num2str(iLon,'%02d') '/fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_201500010001_202100120031_TimeStepsX160.mat'];
      thefileERA5 = [thedirERA5 '/reconstruct_era5_const_tracegas_spectra_geo_rlat' num2str(iLat,'%02d') '_2014_09_2021_08.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif iType == 9
      %% full AIRS 20 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457.mat'];
      %% or if seasonal ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin64/LonBin72/iQAX_3_fits_LonBin72_LatBin64_V1_TimeSteps457_MAM.mat
      if iAllorSeasonal == -1
        thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457_DJF.mat'];
      elseif iAllorSeasonal == -2
        thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457_MAM.mat'];
      elseif iAllorSeasonal == -3
        thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457_JJA.mat'];
      elseif iAllorSeasonal == -4
        thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps457_SON.mat'];
      end
    elseif iType == 10
      %% full AIRS 05 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200200090001_200700080031_TimeStepsX114.mat'];
    elseif iType == 11
      %% full AIRS 10 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200200090001_201200080031_TimeStepsX228.mat'];
    elseif iType == 12
      %% full AIRS 15 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200200090001_201700080031_TimeStepsX342.mat'];
    elseif iType == 13
      %% AIRS 08 year, beginning of mission when SW drifting
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_4_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200200090001_201000080031_TimeStepsX183.mat'];
    elseif iType == 14
      %% AIRS 04 year, end of mission when orbit is beginning to drift
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_201800090001_202200080031_TimeSteps_366_457_X091.mat'];
    elseif iType == 15
      %% AIRS 14 year, overlap with IASI
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps457.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_200800010001_202200120031_TimeSteps_122_464_X342.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif iType == 16
      %% AIRS 4 year
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps498.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_202000070001_202400060030_TimeSteps_408_499_X091.mat'];
    elseif iType == 17
      %% AIRS 21 year starting 2003/07
      %% ../DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin*/LonBin*/iQAX_3_fits_LonBin*_LatBin*_V1_TimeSteps498.mat
      thefilein = [thedir0 '/LonBin' num2str(iLon,'%02d') '/iQAX_3_fits_LonBin' num2str(iLon,'%02d') '_LatBin' num2str(iLat,'%02d') '_V1_TimeSteps498.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif iType == 30
      thefilein = [thedir0 '/trends_anomalies_AMSU_20year_72x64.mat'];
    end

    iBoo = (iLat-1)*72 + iLon;
    if ~exist(thefilein)
      fprintf(1,'oops iLat/iLon = %2i %2i tilenumber %4i file %s DNE \n',iLat,iLon,iBoo,thefilein);
    end

    x = load(thefilein);
    if iType == 8
      xERA5 = load(thefileERA5);
    end
    
    if iType == 30  
      error('already done')
      %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%% AMSU %%%  
      b_asc(iLon,iLat,:) = x.trend_T(:,iLon,iLat);
      b_desc(iLon,iLat,:) = x.trend_T(:,iLon,iLat);
      b_err_asc(iLon,iLat,:) = x.trend_T_err(:,iLon,iLat);
      b_err_desc(iLon,iLat,:) = x.trend_T_err(:,iLon,iLat);
      lagcor_obs_anom_asc(iLon,iLat,:)  = NaN;
      lagcor_obs_anom_desc(iLon,iLat,:) = NaN;
      mean_BT(iLon,iLat,:)  = NaN;

    elseif abs(iType) ~= 3 & iType < 29 
      %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%%
      %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%%
      %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%% AIRS %%%

      %% QUANTILES QUANTILES QUANTILES
      %% these are rads
      %b_asc(iLon,iLat,:) = x.b_asc(:,iQuantile,2);
      %b_desc(iLon,iLat,:) = x.b_desc(:,iQuantile,2);
      %b_err_asc(iLon,iLat,:) = x.berr_asc(:,iQuantile,2);
      %b_err_desc(iLon,iLat,:) = x.berr_desc(:,iQuantile,2);
  
      %% these are bt
      if iQuantile <= 16
        b_asc(iLon,iLat,:) = x.dbt_asc(:,iQuantile);
        b_desc(iLon,iLat,:) = x.dbt_desc(:,iQuantile);
        b_err_asc(iLon,iLat,:) = x.dbt_err_asc(:,iQuantile);
        b_err_desc(iLon,iLat,:) = x.dbt_err_desc(:,iQuantile);
        lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc(:,iQuantile);
        lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc(:,iQuantile);
    
        if exist('xERA5')
          b_cal_desc(iLon,iLat,:) = xERA5.thesave.xtrendSpectral(:,iLon);
          b_cal_err_desc(iLon,iLat,:) = xERA5.thesave.xtrendSpectral_unc(:,iLon);
        end

        mean_rad(iLon,iLat,:) = squeeze(x.b_desc(:,iQuantile,1)); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
        mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,squeeze(x.b_desc(:,iQuantile,1)));
      elseif iQuantile == 50
        iQuantileX = 12:16;
        b_asc(iLon,iLat,:) = nanmean(x.dbt_asc(:,iQuantileX),2);
        b_desc(iLon,iLat,:) = nanmean(x.dbt_desc(:,iQuantileX),2);
        b_err_asc(iLon,iLat,:) = nanmean(x.dbt_err_asc(:,iQuantileX),2);
        b_err_desc(iLon,iLat,:) = nanmean(x.dbt_err_desc(:,iQuantileX),2);
        lagcor_obs_anom_asc(iLon,iLat,:)  = nanmean(x.lag_asc(:,iQuantileX),2);
        lagcor_obs_anom_desc(iLon,iLat,:) = nanmean(x.lag_desc(:,iQuantileX),2);
    
        mean_rad(iLon,iLat,:) = nanmean(squeeze(x.b_desc(:,iQuantileX,1)),2); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
        mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,nanmean(squeeze(x.b_desc(:,iQuantileX,1)),2));
      end    
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    elseif iType == +3
      %% these are extremes      
      b_asc(iLon,iLat,:) = x.dbt_asc;
      b_desc(iLon,iLat,:) = x.dbt_desc;
      b_err_asc(iLon,iLat,:) = x.dbt_err_asc;
      b_err_desc(iLon,iLat,:) = x.dbt_err_desc;
      lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc';
      lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc';
  
      mean_rad(iLon,iLat,:) = x.b_desc(:,1); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
      mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,mean_rad(iLon,iLat,:));
  
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    elseif iType == -3
      %% these are means      
      b_asc(iLon,iLat,:) = x.dbt_asc;
      b_desc(iLon,iLat,:) = x.dbt_desc;
      b_err_asc(iLon,iLat,:) = x.dbt_err_asc;
      b_err_desc(iLon,iLat,:) = x.dbt_err_desc;
      lagcor_obs_anom_asc(iLon,iLat,:)  = x.lag_asc';
      lagcor_obs_anom_desc(iLon,iLat,:) = x.lag_desc';
  
      mean_rad(iLon,iLat,:) = x.b_desc(:,1); %% 10 params = <mean> + <linear trend> t + sum(i=1,4) Ci sin(wi t) + Di cos (wi t)
      mean_BT(iLon,iLat,:)  = rad2bt(h.vchan,mean_rad(iLon,iLat,:));
  
      junk = squeeze(mean_BT(iLon,iLat,:));
      junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
      airs_noiseTtrue(iLon,iLat,:) = junk;

    end

  end
  fprintf(1,'%3i out of 64 \n',iLat);

  if iKeepPlotting > 0
    junk = squeeze(b_desc(:,:,1520));
    scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fnamePROCESS)
  eval(saver);
  saver
else 
  fprintf(1,'%s already exists \n',fnamePROCESS);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps/

figs_put_together_QuantileChoose_trends

figure(1); clf; scatter_coast(X',Y',50,squeeze(b_desc(:,:,1520))'); colorbar; colormap(usa2); caxis([-1 +1]); shading interp; title('dBT1231/dt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /home/sergio/KCARTA/MATLAB
plot(h.vchan,reshape(b_err_desc,72*64,2645)')
plot(h.vchan,reshape(mean_BT,72*64,2645)')
airs_noise = instr_chans2645('airs',2);
for iLat = 1 : 64
  for iLon = 1 : 72
    junk = squeeze(mean_BT(iLon,iLat,:));
    junk = nedt_T0_T1(h.vchan,airs_noise,250*ones(h.nchan,1),real(junk));
    airs_noiseTtrue(iLon,iLat,:) = junk;
  end
end

addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
[g2378,g2645] = compare_goodchans_2378_2645();   %% shows g is index into f-ABCD and NOT chanID
junk = reshape(airs_noiseTtrue,72*64,2645)'/sqrt(120);   %% about 12000 obs/tile/16 days .. so 1% of this is 120 ... noise goes down by sqrt(N)
plot(h.vchan(g2645),junk(g2645,:))
ylim([0 0.3])

plot(h.vchan,nanmean(reshape(b_err_desc,72*64,2645)',2),h.vchan(g2645),nanmean(junk(g2645,:),2)); 
  ylim([0 0.3]); 
  hl = legend('from b_err','from 1/sqrt(N)','location','best','fontsize',10);
%}
