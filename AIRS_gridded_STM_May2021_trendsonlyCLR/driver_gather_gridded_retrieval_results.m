restoredefaultpath
addpath0

% addpath /home/sergio/MATLABCODE
% addpath /home/sergio/MATLABCODE/TROPOPAUSE
% addpath /home/sergio/MATLABCODE/COLORMAP
% addpath /home/sergio/MATLABCODE/COLORMAP/LLS
% addpath /home/sergio/MATLABCODE/PLOTTER
% addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
% addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
% addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
% addpath /home/sergio/MATLABCODE/NANROUTINES/
% addpath /home/sergio/MATLABCODE/SHOWSTATS/
% addpath ../FIND_NWP_MODEL_TRENDS
% addpath /asl/matlib/science
% addpath /asl/matlib/aslutil
% addpath /asl/matlib/h4tools
% addpath /asl/matlib/maps

addpath ../FIND_NWP_MODEL_TRENDS                          %% to get eg get_ERA5_trends_thermodynamic_and_spectral.m
addpath /home/sergio/git/matlabcode/PLOTTER/TILEDPLOTS/   %% to get eg aslmap_2x2tiledlayout.m
addpath /home/sergio/git/matlabcode/SHOWSTATS/            %% to get myhist2d.m
addpath /home/sergio/git/matlabcode/NANROUTINES/          %% to get nanpolyfit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')

if exist('llsmap5.mat')
  load llsmap5  
  %if length(llsmap5) == 64
  %  llsmap5 = llsmap5(2:end,:);
  %end
else
  llsmap5 = usa2;
end

llsmap5NAN = [[0.0 0.0 0.0]; llsmap5];
llsmap5_0 = llsmap5;
llsmap5   = llsmap5NAN;

disp(' ')
disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNorD = input('night (+1) or day (-1)  [+1 = night = default] ');
if length(iNorD) == 0
  iNorD = +1;
end
iDorA = iNorD;

iOCBset = input('obs cal or bias (0,+1,-1) [default 0] : ');
if length(iOCBset) == 0
  iOCBset = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% ABC this may help in do_compute_save_ERA5_feedbacks.m
%% ABC this may help in do_compute_save_ERA5_feedbacks.m
%% ABC this may help in do_compute_save_ERA5_feedbacks.m

%save -v7.3 nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat nwp_spectral_trends_cmip6_era5_airsL3_umbc
save('nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat','-struct','nwp_spectral_trends_cmip6_era5_airsL3_umbc','-v7.3');
vars_cmip6_era5_airsL3_umbc = whos('-file','nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');
era5rates = load('nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat','era5_100_layertrends');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iAK = input('do AvgKernels ???? (-1/+1 = yes = default) ');
if length(iAK) == 0
  iAK = +1;
end
if iAK > 0
  %% see Plotutils/plot_retrieval_latbins_fewlays
  %% SEARCH for  << CAN CUT AND PASTE THIS TO RERUN >>
  
  iNumYears = input('  Enter number of years from 2002-X (so we can load in appropriate ERA5 trend file !!!! [default = 23, can also give =-N for data ending in 2022 (end of mission)] ');  
  if length(iNumYears) == 0
    iNumYears = 20;
    iNumYears = 23;    
  end
  %junk = ['../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
  %junk = ['../FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
  %era5 = load(junk);
  
  if iNumYears > 23
    disp('ooer have not done ERA5 trends for > 23 years, just set to 23')
  end
  era5 = get_ERA5_trends_thermodynamic_and_spectral(min(iNumYears,23),iNorD);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iJunk = input('clear 50 figs???? (-1 = no = default) ');
if length(iJunk) == 0
  iJunk = -1;
end
if iJunk > 0
  for ii = 1 : 50
    figure(ii); colormap jet
  end
else
  figure(6);  clf; %% this puts d/dt ST plot with correct title
  figure(28); clf; %% this puts d/dt ST plot with correct title
  figure(29); clf; %% this puts d/dt ST plot with correct title
  figure(30); clf; %% this puts d/dt ST plot with correct title
end

iNumLay = 20;
iNumLay = input('Enter number of expected fat layers (20 [5 thick AIRS layers], having good luck with (default 49 [2 thick AIRS layers]) : ');
if length(iNumLay) == 0
  iNumLay = 20;
  iNumLay = 49;
end
iNavgLay = floor(100/iNumLay);
iWVind = 07:26; iWVind = (1:iNumLay)+6+0*iNumLay;
iTind  = 27:46; iTind  = (1:iNumLay)+6+1*iNumLay;
iO3ind = 47:66; iO3ind = (1:iNumLay)+6+2*iNumLay;

disp('quants for dataset 1-8            = [0 0.01 0.02 0.03 0.04 0.05 0.10 0.25 0.50 0.75 0.9 0.95 0.96 0.97 0.98 0.99 1.00]');
disp('quants for dataset 9,10,11,12     = [0.50 0.80 0.90 0.95 0.97 1.00]');
disp('quants for dataset 13             = [0.03 0.50 0.97]');
disp('quants for dataset 14,15,16,17    = [0.50 0.80 0.90 0.95 0.97 1.00]');
disp('quants for dataset 30 [AMSU only] = [1]');
dataset = input('Enter \n (+1) Strow 2002/09-2020/08 Q1-16 \n (-1) Sergio 2002/09-2020/08 Q1-16 \n (2) Sergio 2002/09-2021/07 OLD  Q1-16 \n (3) Sergio 2002/09-2021/08 Extreme \n (-3) Sergio 2002/09-2021/08 Mean \n (4) Sergio 2002/09-2021/08 FULL  Q1-16 \n (5) Sergio 2002/09-2014/08 CMIP6  Q1-16 \n (6) Sergio 2012/05-2019/04 CrIS NSR overlap \n (7) Sergio 2002/09-2022/08 20 YEARS \n (8) Sergio 2015/01-2021/12 OCO2 overlap \n (9,10,11,12) Sergio 2002/09-2022/2007/2012/2017/08 20 YEARS new quants iQAX=3 \n (18,19) Sergio 2002/09-2025/08 23 YEARS (20) for Ryan 2003/01-2022/12  \n  new quants iQAX=3  19,20 uses new chip disks \n  :::  [19 = Default] : ');
if length(dataset) == 0
  dataset = 19;
end

if dataset == 1 | dataset == -1
  iNumYears = 18;
elseif dataset == 2
  iNumYears = 19;
elseif abs(dataset) == 3
  iNumYears = 19;
elseif dataset == 4
  iNumYears = 19;
elseif dataset == 5
  iNumYears = 12;
elseif dataset == 6
  iNumYears = 07;
elseif dataset == 7
  iNumYears = 20;
elseif dataset == 8
  iNumYears = 07;
elseif dataset == 9
  iNumYears = 20;
elseif dataset == 10
  iNumYears = 05;
elseif dataset == 11
  iNumYears = 10;
elseif dataset == 12
  iNumYears = 15;
elseif dataset == 16
  iNumYears = 04;
elseif dataset == 17
  iNumYears = 22;
elseif dataset == 18
  iNumYears = 23;
elseif dataset == 19
  iNumYears = 23;
elseif dataset == 30  %% AMSU
  iNumYears = 20;
end

if dataset == 30
  iQuantile = 01;
elseif dataset == -3
  iQuantile = 00;
elseif dataset >= 9 & dataset < 30
  iQuantile = 05; 
  if iOCBset == 0
    iQuantile = input('Dataset = 9,10,11,12,18,19,20  iOCBset = 0 (obs)  ==> Which quantile 1..5   [3 = Default] : ');
    if length(iQuantile) == 0
      iQuantile = 3;
    end
  elseif iOCBset == 1
    iQuantile = input('Dataset = 9,10,11,12,18,19,20   iOCBset = 1 (cal)  ==> Which quantile 1..16   [3 = Default] : ');
    if length(iQuantile) == 0
      iQuantile = 3;
    end
  else
    iOCBset
    error('huh??? iOCBset = 0,1 only!!!')
  end
elseif dataset ~= 3 & dataset < 9
  iQuantile = 16;  %% AIRS STM 2021, hottest
  iQuantile = 08;  %% 
  iQuantile = input('Which quantile -1 for extremes [(1--16) (99 for orig Q16, done for AIRS STM)]   [16 = Default] : ');
  if length(iQuantile) == 0
    iQuantile = 16;
  end
else
  iQuantile = [];
end

if iOCBset == 0
  if dataset == 30
    %% AMSU
    data_trends = load(['iType_30_AMSU_iQAX_01.mat']);
  elseif dataset == 3
    data_trends = load(['iType_3_extreme_convert_sergio_clearskygrid_obsonly.mat']);
  elseif dataset == -3
    data_trends = load(['iType_-3_mean_convert_sergio_clearskygrid_obsonly.mat']);
  elseif dataset == 1
    if iQuantile >= 1 & iQuantile <= 16
      data_trends = load(['convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    elseif iQuantile == 99
      data_trends = load(['convert_sergio_clearskygrid_obsonly_Q16.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == -1
    if iQuantile >= 1 & iQuantile <= 16
      data_trends = load(['iType_-1_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    elseif iQuantile == 99
      data_trends = load(['XYZconvert_sergio_clearskygrid_obsonly_Q16.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 2
    if iQuantile >= 1 & iQuantile <= 16
      data_trends = load(['iType_2_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    elseif iQuantile == 99
      data_trends = load(['XYZconvert_sergio_clearskygrid_obsonly_Q16.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 4
    if iQuantile >= 1 & (iQuantile <= 16 | iQuantile == 50)
      data_trends = load(['iType_4_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 5
    if iQuantile >= 1 & (iQuantile <= 16 | iQuantile == 50)
      data_trends = load(['iType_5_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 6
    if iQuantile >= 1 & (iQuantile <= 16 | iQuantile == 50)
      data_trends = load(['iType_6_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 7
    if iQuantile >= 1 & (iQuantile <= 16 | iQuantile == 50)
      data_trends = load(['iType_7_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 8
    if iQuantile >= 1 & (iQuantile <= 16 | iQuantile == 50)
      data_trends = load(['iType_8_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 9
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 10
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_10_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 11
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_11_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 12
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_12_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 16
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_16_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 17
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_17_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  elseif dataset == 18 | dataset == 19 | dataset == 20
    if iQuantile >= 1 & (iQuantile <= 5 | iQuantile == 50)
      data_trends = load(['iType_' num2str(dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
    else
      fprintf(1,'oops error trying to read in data trends for iOCBset = %2i dataset = %2i iQuantile = %2i \n',iOCBset,dataset,iQuantile);  error('yuk yuk')
    end
  else
    dataset
    error('huh unknown dataset')
  end
elseif iOCBset == 1
  %% OLD
  % for iibin = 1 : 64
  %   strlatbin = num2str(iibin,'%02d');
  %   datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];
  %   datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/'  strY '/reconstruct_era5_spectra_geo_rlat' strlatbin '_2002_09_2022_08.mat'];
  %   junk = load(datafile);
  %   junkind = (1:72) + (iibin-1)*72;
  %   %data_trends.b_desc(1:72,iibin,:)     = real(junk.thesave.xtrend)';
  %   %data_trends.b_err_desc(1:72,iibin,:) = real(junk.thesave.xtrendErr)';
  %   data_trends.b_desc(1:72,iibin,:)     = real(junk.thesave.xtrendSpectral)';
  %   data_trends.b_err_desc(1:72,iibin,:) = real(junk.thesave.xtrendSpectral_unc)';
  % end

  %% SEARCH for  << CAN CUT AND PASTE THIS TO RERUN >>
  if dataset == 9 | dataset == 20
    %% set 9  is 2002/09 to 2022/08
    %% set 20 is 2003/01 to 2022/12
    strY = '20yrs';
  elseif dataset == 19
    strY = '23yrs';
  end
  if dataset == 9 | dataset == 20
    datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/'  strY '/ERA5_spectraltrends_2002_09_2022_08.mat'];
  elseif dataset == 19
    datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/'  strY '/ERA5_spectraltrends_2002_09_2025_08.mat'];
  end
  junk = load(datafile);
  data_trends.b_desc = real(junk.era5_rates_desc);
  data_trends.b_err_desc = real(junk.era5_rates_desc_unc);
  junk = load('h2645structure.mat');
  data_trends.h.ichan = junk.h.ichan;
  data_trends.h.vchan = junk.h.vchan;
  do_XX_YY_from_X_Y
  figure(6); clf; scatter_coast(XX,YY,50,data_trends.b_desc(1520,:)); colormap(usa2); caxis([-0.15 +0.15]); title('dBT1231/dt'); pause(0.1)
end

iNumChan = 2645;
if dataset == 30
  iNumChan = 13;
end

fnamelastloaded = 'none';
if ~exist('iaFound')
  %% SEARCH for  << CAN CUT AND PASTE THIS TO RERUN >>
  %% along with eg
  %{
  %% get ERA5 thermodynamic rates  
  iNumYears = 20;
  era5 = get_ERA5_trends_thermodynamic_and_spectral(min(iNumYears,23),iNorD);
  
  and if you want, load in spectral raates
  strY = '20yrs'; dataset = 9;     datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/'  strY '/ERA5_spectraltrends_2002_09_2022_08.mat'];
  junk = load(datafile);
  data_trends.b_desc = real(junk.era5_rates_desc);
  data_trends.b_err_desc = real(junk.era5_rates_desc_unc);
  junk = load('h2645structure.mat');
  data_trends.h.ichan = junk.h.ichan;
  data_trends.h.vchan = junk.h.vchan;
  do_XX_YY_from_X_Y
  %}
  
  clear results*
  iaFound = zeros(1,4608);
  existfname = zeros(1,4608);

  save_cov_set.cov_set = nan(13,4608);
  save_cov_set.fmat = nan(6,4608);

  cdofs = nan(4608,66);
  lencdofs = nan(1,4608);
  thedofs = nan(1,4608);

  results = nan(4608,6);
  resultsWV = nan(4608,iNumLay);
  resultsT  = nan(4608,iNumLay);
  resultsO3 = nan(4608,iNumLay);

  resultsunc = nan(4608,6);
  resultsWVunc = nan(4608,iNumLay);
  resultsTunc  = nan(4608,iNumLay);
  resultsO3unc = nan(4608,iNumLay);

  spectral_deltan00 = nan(iNumChan,4608);
  rates             = nan(iNumChan,4608);
  fits              = nan(iNumChan,4608);
  componentfits     = nan(5,iNumChan,4608);
  nedt              = nan(iNumChan,4608);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dataset < 30
  %moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','b_err_desc');
  moonoise = load('iType_19_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat','b_err_desc');  
elseif dataset == 30
  moonoise = load('iType_30_AMSU_iQAX_01.mat','b_err_desc');
  moonoise.b_err_desc = moonoise.b_err_desc(:,:,1:13);
end
b_err_desc = moonoise.b_err_desc; clear moonoise;
b_err_desc = permute(b_err_desc,[3 1 2]);
b_err_desc = reshape(b_err_desc,iNumChan,72*64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iWarning = 0;
clear fname

iCnt = 0;
for ii = 1 : 72 : 72*64
  if iOCBset == 0
    if dataset ~= 3
      if iNorD > 0
        fname = ['/asl/s1/sergio/Tiles4608/Output_WORKS_May18_2021_Great_AIRS_STM/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here after July 2021
        fname = ['Output/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
      elseif iNorD < 0
        fname = ['Output_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021
      end
    elseif dataset == 3
      if iNorD > 0
        fname = ['Output/Extreme/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
      elseif iNorD < 0
        fname = ['Output_Day/Extreme/test' num2str(ii) '.mat']; %% stored here before before July 2021
      end
    elseif dataset == -3
      if iNorD > 0
        fname = ['Output/Quantile00/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
      elseif iNorD < 0
        fname = ['Output_Day/Quantile00/test' num2str(ii) '.mat']; %% stored here before before July 2021
      end
    end
  elseif iOCBset == 1
    if iNorD > 0
      fname = ['Output_CAL/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['Output_Day_CAL/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021
    end
  end

  if exist(fname) > 0
    iCnt = iCnt + 1;
    junkdir = dir(fname);
    fprintf(1,'%s %s \n',[junkdir.folder     fname],junkdir.date);
  else
    fprintf(1,'%s DNE \n',fname);  
  end
end
fprintf(1,'found %4i of 64 (subset) files \n',iCnt);
if iCnt == 0
  fprintf(1,'with iOCBset = %2i dataset = %2i iNorD = %2i seem to have found nothing nada zilch in %s \n',iOCBset,dataset,iNorD,fname )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iJunk = input('correct dates/names etc etc???  Proceed or quit (+1 default/-1) : ');
if length(iJunk) == 0
  iJunk = +1;
end
if iJunk < 0
  return
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);

clear pavg

iLoopXY = +1; %% loop over lonbins, plot zonal avg once after all 64 lats,72 lons tried to be read
iLoopXY = -1; %% loop over latbins, plot zonal avg each time after all 64 lats tried to be read, so 72 plots DEFAULT

iDoAgain = +1;
if iLoopXY > 0
  read_loop_latbins
else
  read_loop_lonbins
end

%%%%%%%%%%
%% quick plots
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
aslmap(6,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt so far');     caxis([-1 +1]*0.15); colormap(llsmap5)
%figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
%figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,1:iNumLay,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.015)
figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dT/dt');      colormap(llsmap5); caxis([-1 +1]*0.15)
figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(rlat,pavg,waha);  shading interp; colorbar('horizontal'); set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.01)
  figure(29); set(gca,'yscale','log'); ylim([10 1000]);       figure(30); set(gca,'yscale','linear'); ylim([100 1000]);
figure(31); clf; waha = reshape(iaFound,72,64);                                                   pcolor(rlon,rlat,waha'); shading flat;   colorbar; set(gca,'ydir','normal');  title('read in so far');  xlabel('Longitude'); ylabel('Latitude'); colormap(jet);

%% show ERA5
figure(1); figure(2); figure(4); 
disp('all 4608 read in .... RET to  continue'); pause

%%%%%%%%%%

disp('WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster')
disp('WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster')
disp('WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster')

a = load(fnamelastloaded);
fprintf(1,'last loaded file %s has xb(1:6) = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n',fnamelastloaded,a.oem.xb(1:6))

%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 12;  figure(ii); clf; end;

wah = resultsTunc'; wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2));  
  figure(1); clf; pcolor(rlat,pavg,wah); title('\sigma T'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp; 
  set(gca,'yscale','log'); ylim([1 1000])
wah = resultsWVunc'; wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  figure(2); clf; pcolor(rlat,pavg,wah); title('\sigma WV'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp
wah = resultsO3unc'; wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  figure(3); clf; pcolor(rlat,pavg,wah); title('\sigma O3'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp

% load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
% load latB64.mat
% rlat65 = latB2; rlon73 = -180 : 5 : +180;
% rlon = -180 : 5 : +180;  rlat = latB2; 
% rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
% rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
% [Y,X] = meshgrid(rlat,rlon);
% X = X; Y = Y;
do_XX_YY_from_X_Y

show_S_N
disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
% [salti, landfrac] = usgs_deg10_dem(Y(:),X(:));

% see more /home/sergio/git/rtpmake/CLUST_RTPMAKE/COMMON_SETTINGS/set_landfrac_using_L1B_L1C_or_usgs.m
addpath /home/sergio/git/matlabcode/DEM_DigitalELeveationModel
[salti,landfrac,gebco] = gdemm_dem_and_imerg_lf(Y(:),X(:));

Ylat = Y(:);
Xlon = X(:);
% save landfrac_mask4608.mat landfrac Ylat Xlon rlat65 rlon73 rlon rlat

if dataset < 30
  junk = load('h2645structure.mat');
  f    = junk.h.vchan;
  
  i1419 = find(f >= 1419,1);
  i1231 = find(f >= 1231,1);
  i0900 = find(f >= 0900,1);
  
  clf;; scatter_coast(Xlon,Ylat,50,results(:,1)); title('d/dt CO2');  caxis([1.5 2.5]); caxis([2.0 2.5])
  aslmap(4,rlat65,rlon73,smoothn((reshape(results(:,1),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('d/dt CO2');  caxis([1.5 2.5])
  
  figure(29); waha = squeeze(nanmean(reshape(resultsT,72,64,iNumLay),1)); waha = waha';        pcolor(waha); shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dT/dt'); colormap(llsmap5); caxis([-1 +1]*0.15)
  figure(30); waha = squeeze(nanmean(reshape(resultsWV,72,64,iNumLay),1)); waha = waha';       pcolor(waha); shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.015)
  
  if dataset == 8
    aslmap(6,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
    eracal = load('saved_co2_era_cal_retrieval_2002_2021.mat');
    eracal = load('saved_co2_era_cal_retrieval_2015_2021.mat');
    aslmap(2,rlat65,rlon73,smoothn((reshape(eracal.results(:,1),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('ERA d/dt CO2');  caxis([-1 +1]*2); 
    aslmap(3,rlat65,rlon73,smoothn((reshape(results(:,1),72,64)'),1),[-90 +90],[-180 +180]);        colormap(jet);      title('UMBC d/dt CO2');  caxis([0 4]);
    aslmap(4,rlat65,rlon73,smoothn((reshape(results(:,1)-eracal.results(:,1),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('UMBC-ERA d/dt CO2');  caxis([0 4]);
    oco2 = load('oco2_timeseries.mat');
    aslmap(35,rlat65,rlon73,smoothn(oco2.co2_trend',1),[-90 +90],[-180 +180]); colormap(jet);  title('d/dt CO2 from OCO2 2015-2021');  caxis([2.0 3.0])  
      caxis([2.35 2.55])
    disp('showed the CO2 trends ... ret to continue'); pause
  end
  
  aslmap(5,rlat65,rlon73,smoothn((reshape(rates(i1231,:),72,64)'),1), [-90 +90],[-180 +180]); title('dBT1231/dt'); caxis([-1 +1]*0.15); colormap(llsmap5)
  aslmap(6,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
  
  aslmap(7,rlat65,rlon73,smoothn((reshape(rates(i1419,:),72,64)'),1), [-90 +90],[-180 +180]); title('dBT1419/dt'); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(8); junk = smoothn((reshape(rates(i1419,:),72,64)')); junk = reshape(rates(i1419,:),72,64)'; plot(rlat,nanmean(junk,2)); title('Zonal dBT1419/dt'); grid; xlim([-90 +90])
end
  
figure(9); scatter_coast(Xlon,Ylat,50,thedofs); jett = jet(64); jett(1,:) = 1; colormap(jett); title('ALL DOFS'); caxis([0 max(thedofs)]); caxis([0 30])
figure(10); boo = find(lencdofs == 66); sumc = mean(cdofs(boo,:),1); 
  plot(cumsum(sumc)); fl = ceil(sum(sumc)); line([6 6],[0 fl],'color','k'); line([26 26],[0 fl],'color','k'); line([46 46],[0 fl],'color','k');
  text(2,10,'TG','fontsize',10);   text(16,10,'WV','fontsize',10);   text(36,10,'T','fontsize',10);   text(56,10,'O3','fontsize',10); xlim([0 66])

for ii = 1 : iNumLay
  ind = (1:iNavgLay) + (ii-1)*iNavgLay;
  pjunk20(ii) = mean(plays(ind));
end
cssumc = cumsum(sumc);

figure(11); semilogy(cssumc(iWVind)-cssumc(7),pjunk20,cssumc(iTind)-cssumc(27),pjunk20,cssumc(iO3ind)-cssumc(47),pjunk20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([0.1 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')
figure(11); plot(cssumc(iWVind)-cssumc(7),pjunk20,cssumc(iTind)-cssumc(27),pjunk20,cssumc(iO3ind)-cssumc(47),pjunk20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([50 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')
figure(11); plot(cssumc(iWVind)-cssumc(7),fliplr(pjunk20),cssumc(iTind)-cssumc(27),fliplr(pjunk20),cssumc(iO3ind)-cssumc(47),fliplr(pjunk20),'linewidth',2)
  set(gca,'ydir','reverse'); ylim([50 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')

pflip20 = fliplr(pjunk20);
wvsumc     = sumc(iWVind);         semilogy(wvsumc,pjunk20); set(gca,'ydir','reverse'); ylim([0.01 1000])
wvsumcflip = fliplr(sumc(iWVind)); semilogy(wvsumc,pjunk20,'o-',wvsumcflip,pflip20,cumsum(wvsumcflip),pflip20); set(gca,'ydir','reverse'); ylim([0.01 1000])

wvsumcflip = cumsum(fliplr(sumc(iWVind)));
tsumcflip  = cumsum(fliplr(sumc(iTind)));
o3sumcflip = cumsum(fliplr(sumc(iO3ind)));
figure(11); plot(wvsumcflip,pflip20,tsumcflip,pflip20,o3sumcflip,pflip20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([50 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')
figure(11); semilogy(wvsumcflip,pflip20,tsumcflip,pflip20,o3sumcflip,pflip20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([0.050 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')

if dataset == 30
  load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/hAMSU.mat');
  f = hAMSU.vchan;
  f = f(1:13);
end
figure(12); clf; plot(f,nanmean(rates,2),'b',f,nanstd(rates,[],2),'c--',f,nanmean(rates-fits,2),'r',f,nanstd(rates-fits,[],2),'m--');
  plotaxis2; hl = legend('mean obs','std obs','mean(obs-fits)','std(obs-fits)','fontsize',10,'location','best');

disp('ret to continue to spectral chisqr'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'checking noise : nansum(nansum(nedt(jacobian.chanset,:)-b_err_desc(jacobian.chanset,:))) = %8.6f \n',nansum(nansum(nedt(jacobian.chanset,:)-b_err_desc(jacobian.chanset,:))))

% iIgnoreChans_N2O = -1; %% retrieve using N2O chans
% iIgnoreChans_N2O = +1; %% ignore using N2O chans
settings.iIgnoreChans_CH4 = -1;
settings.iIgnoreChans_N2O = -1;
settings.iIgnoreChans_SO2 = -1;
chanset = jacobian.chanset;

read_fileMean17years

plotopt.iUpperWavenumLimit = 1620;
plotopt.rlon = pMean17years.rlon;
plotopt.rlat = pMean17years.rlat;
if ~isfield(settings,'iInstr')
  settings.iInstr = 1;
end
if dataset ~= 30
  [raaBadFoav,indBadFov,chisqrX,chisqrR] = plot_spectral_region_chisqr(rates(chanset,:),0*rates(chanset,:),0*rates(chanset,:),fits(chanset,:),f(chanset,:),nedt(chanset,:),-1,settings,plotopt);
  figure(11); ylim([-1 +1]*0.1/2)
  figure(12); ylim([-1 +1]*5)
  for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*10); end
end
if dataset == 30
  junk = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/AMSU_12channels_20years/SARTA_CALCS/amsu_sarta_trends.mat');
  era5.era5_spectral_rates2645 = era5.era5_spectral_rates;
  era5.era5_spectral_rates     = junk.trend_sarta;
  chisqrX = zeros(1,4608);
  chisqrR.iAll    = chisqrX;
  chisqrR.i15um   = chisqrX;
  chisqrR.iWindow = chisqrX;
  chisqrR.iWV     = chisqrX;
end

figure(21); clf
if iNumYears >= 0
  iNumYearsX = roundN(iNumYears,5);  %% only have spectra for 5,10,15,20 ie 2022-2007,2012,2017,2022
  if iNumYearsX == 0
    iNumYearsX = 5;
  end
else
  iNumYearsX = -4; %% only have spectra for -4   ie 2018-2022
end
if iNumYears ~= iNumYearsX
  fprintf(1,' WARNING : get_ERA5_trends_thermodynamic_and_spectral.m using ERA5 trends from %2i years instead of %2i years \n',iNumYearsX,iNumYears);
end
if isfield(oem,'spectral_deltan00')
  plot(f,nanmean(rates'),'b',f,nanmean(spectral_deltan00'),'g',f,nanmean(era5.era5_spectral_rates'),'r',...
       f,nanstd(rates'),'b--',f,nanstd(spectral_deltan00'),'g--',f,nanstd(era5.era5_spectral_rates'),'r--'); plotaxis2;
  axis([640 1640 -0.4 +0.4])
  if iNumYears >= 0
    hl = legend('actual data','what was fitted',['ERA5 ' num2str(iNumYearsX) 'yrs, vary tracegas'],'location','best');
  elseif iNumYears < 0
    hl = legend('actual data','what was fitted',['ERA5 ' num2str(iNumYearsX) 'yrs, const tracegas'],'location','best');
  end
  title('After doing data-sum(jac(i)*startxb(i)')
  ylabel('solid : mean; dashed : std [K/yr]');
  xlabel('Wavenumber cm-1');
end

%[pMean17years.salti,pMean17years.landfrac] = usgs_deg10_dem(pMean17years.rlat,pMean17years.rlon);
[pMean17years.salti,pMean17years.landfrac] = gdemm_dem_and_imerg_lf(pMean17years.rlat,pMean17years.rlon);
figure(06); clf; aslmap(06,rlat65,rlon73,smoothn((reshape(results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt UMBC');     caxis([-1 +1]*0.15); colormap(llsmap5)
figure(50); clf; aslmap(50,rlat65,rlon73,smoothn((reshape(era5.trend_stemp',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt ERA5'); caxis([-1 +1]*0.15); colormap(llsmap5)
junk = find(pMean17years.landfrac == 0); fprintf(1,'SKT rates ocean : ERA5 = %8.4f    UMBC = %8.4f K/yr \n',nanmean(era5.trend_stemp(junk)),nanmean(results(junk,6)))
junk = find(pMean17years.landfrac == 1); fprintf(1,'SKT rates land  : ERA5 = %8.4f    UMBC = %8.4f K/yr \n',nanmean(era5.trend_stemp(junk)),nanmean(results(junk,6)))
junk = find(pMean17years.landfrac >= 0); fprintf(1,'SKT rates all   : ERA5 = %8.4f    UMBC = %8.4f K/yr \n',nanmean(era5.trend_stemp(junk)),nanmean(results(junk,6)))
figure(51); clf; plot(rlat,nanmean(reshape(results(:,6),72,64),1),'b',rlat,nanmean(reshape(era5.trend_stemp(junk),72,64),1),'r','linewidth',2); plotaxis2; xlabel('latitude'); ylabel('dSKT/dt K/yr'); hl = legend('umbc','era5','location','best');
  disp('ret to continue'); pause

%{
rlat4608 = p.rlat; rlon4608 = p.rlon; landfrac4608 = p.landfrac;
era5spectral = era5.era5_spectral_rates;
save strow_2028_2022_JPLreport.mat f rates spectral_deltan00 era5spectral iNumYearsX
clear era5spectral *4608
%}

maskLF = ones(size(chisqrX));
figure(52); clf; 
aslmap(52,rlat65,rlon73,smoothn((reshape(maskLF.*chisqrX,72,64)') ,1), [-90 +90],[-180 +180]); title('Chisqr');     caxis([0 +1]*0.2); colormap(jet)
clear plotoptions;
plotoptions.cx = [0 +1]*0.15; plotoptions.maintitle = '\chi^2'; plotoptions.cmap = jet(64);
plotoptions.str11 = 'ALL';   
plotoptions.str12 = 'T(z)';   
plotoptions.str21 = 'Window';   
plotoptions.str22 = '5*WV(z)';   
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;

z11 = maskLF.*chisqrR.iAll;
z12 = maskLF.*chisqrR.i15um;
z21 = maskLF.*chisqrR.iWindow;
z22 = maskLF.*chisqrR.iWV;

iFig = 52; figure(iFig); clf; 
aslmap_2x2tiledlayout(z11,z12,z21,z22*5,iFig,plotoptions);
disp('Fig 52 : chisqr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('ret to continue to gridded results'); pause

%% look at current AIRS UMBC trend retrievals
show_unc

%%% <<<< THE BIG MAJOR ONE >>>>
%%% <<<< THE BIG MAJOR ONE >>>>
%%% <<<< THE BIG MAJOR ONE >>>>
  plot_driver_gather_gridded_retrieval_results   
%%% <<<< THE BIG MAJOR ONE >>>>
%%% <<<< THE BIG MAJOR ONE >>>>
%%% <<<< THE BIG MAJOR ONE >>>>

%% give a chance for a quick save, can save this sit back .. then call eg 
%%   driver_compute_feedbacks_from_smallFATfile.m
savesmallFATfile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% then look at model data or just do OLR trends and feedbacks
disp('if you really only want feedbacks then all you have to do is run "simple_get_the_model_trends_do_feedbacks" which gets the model trends, then runs do_feedbacks');
iX = input('do all the complicated stuff (+1) or just the simple stuff/feedbacks (default, -1) : ');
if length(iX) == 0
  iX = -1;
end

if iX == +1
  look_at_other_model_data
elseif iX == -1
  simple_get_the_model_trends_do_feedbacks
end
disp('if you save a big fat file, you can then look at SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/wrapper_driver_gather_ALL_rates_AIRSL3_NWP_XMIP.m to make the global plots');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSave = input('Enter to save  (-1, default) small file (0) nothing (+1) big fat file : ');
if length(iSave) == 0
  iSave = -1;
end
if iSave == +1
  savebigFATfile
elseif iSave == -1
  savesmallFATfile
end
