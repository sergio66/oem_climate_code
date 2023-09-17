addpath /asl/matlib/maps/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

pavg = load('/home/sergio/MATLABCODE/airslevels.dat');
pavg = flipud(plevs2plays(pavg));

load llsmap5

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
% YY = Y(:);
% XX = X(:);
XX = X'; XX = XX(:); XX = XX';
YY = Y'; YY = YY(:); YY = YY';

iCosWgt = +1;
if iCosWgt > 0
  % YY = Y(:)'; 
  YY = Y'; YY = YY(:); YY = YY';
  YY = cos(YY*pi/180);
else
  YY = ones(1,4608);
end

N_01 = load('ERA5_atm_N_cld_data_2002_09_to_2003_08_trends_desc.mat');
N_02 = load('ERA5_atm_N_cld_data_2002_09_to_2004_08_trends_desc.mat');
N_03 = load('ERA5_atm_N_cld_data_2002_09_to_2005_08_trends_desc.mat');
N_04 = load('ERA5_atm_N_cld_data_2002_09_to_2006_08_trends_desc.mat');
N_05 = load('ERA5_atm_N_cld_data_2002_09_to_2007_08_trends_desc.mat');
N_06 = load('ERA5_atm_N_cld_data_2002_09_to_2008_08_trends_desc.mat');
N_07 = load('ERA5_atm_N_cld_data_2002_09_to_2009_08_trends_desc.mat');
N_08 = load('ERA5_atm_N_cld_data_2002_09_to_2010_08_trends_desc.mat');
N_09 = load('ERA5_atm_N_cld_data_2002_09_to_2011_08_trends_desc.mat');
N_10 = load('ERA5_atm_N_cld_data_2002_09_to_2012_08_trends_desc.mat');
N_11 = load('ERA5_atm_N_cld_data_2002_09_to_2013_08_trends_desc.mat');
N_12 = load('ERA5_atm_N_cld_data_2002_09_to_2014_08_trends_desc.mat');
N_13 = load('ERA5_atm_N_cld_data_2002_09_to_2015_08_trends_desc.mat');
N_14 = load('ERA5_atm_N_cld_data_2002_09_to_2016_08_trends_desc.mat');
N_15 = load('ERA5_atm_N_cld_data_2002_09_to_2017_08_trends_desc.mat');
N_16 = load('ERA5_atm_N_cld_data_2002_09_to_2018_08_trends_desc.mat');
N_17 = load('ERA5_atm_N_cld_data_2002_09_to_2019_08_trends_desc.mat');
N_18 = load('ERA5_atm_N_cld_data_2002_09_to_2020_08_trends_desc.mat');
N_19 = load('ERA5_atm_N_cld_data_2002_09_to_2021_08_trends_desc.mat');
N_20 = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat');

for ii = 1 : 20
  str = ['NX = N_' num2str(ii,'%02i') ';'];
  eval(str)

  clear junk
  junk.pavg = pavg;
     junk.stemprate = NX.trend_stemp .* maskLF;
     junk.ptemprate = NX.trend_ptemp .* maskLF100;
     junk.waterrate = NX.trend_gas_1 .* maskLF100;
     junk.RHrate    = NX.trend_RH .* maskLF100;
     junkamp = atmospheric_amplification(junk);
  mean_stemprate(ii)      = nansum(NX.trend_stemp .* YY)/sum(YY);
  mean_err_stemprate(ii)  = nanstd(NX.trend_stemp .* YY);  
  mean_err_stemprate2(ii) = nansum(NX.trend_stemp_err .* YY)/sum(YY);
 
  aslmap(1,rlat65,rlon73,smoothn((reshape(NX.trend_stemp,72,64))',1),[-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*1)
  title(['avg over N = ' num2str(ii) ' years'])

  maskLF = ones(1,4608);
  maskLF100 = ones(100,1) * ones(1,4608);

  era5ampT(ii,:) = junkamp.T_amp';
  era5ampRH(ii,:) = junkamp.RH_amp';
  era5ampWV(ii,:) = junkamp.WVfrac_percent_amp';

  pause(1)
end

iNum = 20;

figure(1); clf; errorbar(1:iNum,mean_stemprate,mean_err_stemprate, 'bo-','linewidth',2); plotaxis2; ylim([-1 +1])
figure(1); clf; errorbar(1:iNum,mean_stemprate,mean_err_stemprate2,'bo-','linewidth',2); plotaxis2; ylim([-1 +1]/5)
figure(2); semilogy(era5ampT,pavg(1:97),'linewidth',2); set(gca,'ydir','reverse'); plotaxis2; hl = legend(num2str((1:20)')); title('T amplification');
figure(2); h = plot(era5ampT,pavg(1:97),'linewidth',2); set(gca,'ydir','reverse');  plotaxis2; hl = legend(num2str((1:20)'),'fontsize',8);    title('T amplification');  set(h, {'color'}, num2cell(jet(20), 2));
figure(3); h = plot(era5ampWV,pavg(1:97),'linewidth',2); set(gca,'ydir','reverse'); plotaxis2; hl = legend(num2str((1:20)'),'fontsize',8);    title('WV amplification'); set(h, {'color'}, num2cell(jet(20), 2));
figure(4); h = plot(era5ampRH,pavg(1:97),'linewidth',2); set(gca,'ydir','reverse'); plotaxis2; hl = legend(num2str((1:20)'),'fontsize',8);    title('RH amplification'); set(h, {'color'}, num2cell(jet(20), 2));

i900 = find(pavg >= 900,1);
i800 = find(pavg >= 800,1);
i600 = find(pavg >= 600,1);
i300 = find(pavg >= 300,1);

iS = 3;
iS = 10;
iS = 5;
figure(5); plot(1:iNum,era5ampT(:,i900), 1:iNum,era5ampT(:,i800), 1:iNum,era5ampT(:,i600), 1:iNum,era5ampT(:,i300),'linewidth',2); xlabel('Years'); ylabel('T amplification');   hl = legend('900 mb','800 mb','600 mb','300 mb','fontsize',8,'location','best');
figure(6); plot(1:iNum,era5ampRH(:,i900),1:iNum,era5ampRH(:,i800),1:iNum,era5ampRH(:,i600),1:iNum,era5ampRH(:,i300),'linewidth',2); xlabel('Years'); ylabel('RH amplification'); hl = legend('900 mb','800 mb','600 mb','300 mb','fontsize',8,'location','best');;
figure(7); plot(1:iNum,era5ampWV(:,i900),1:iNum,era5ampWV(:,i800),1:iNum,era5ampWV(:,i600),1:iNum,era5ampWV(:,i300),'linewidth',2); xlabel('Years'); ylabel('WV amplification'); hl = legend('900 mb','800 mb','600 mb','300 mb','fontsize',8,'location','best');;
figure(5); plot(1:iNum,smooth(era5ampT(:,i800),iS), 1:iNum,smooth(era5ampT(:,i600),iS), 1:iNum,smooth(era5ampT(:,i300),iS),'linewidth',2); xlabel('Years'); ylabel('T amplification');   hl = legend('800 mb','600 mb','300 mb','fontsize',8,'location','best'); plotaxis2;
figure(6); plot(1:iNum,smooth(era5ampRH(:,i800),iS),1:iNum,smooth(era5ampRH(:,i600),iS),1:iNum,smooth(era5ampRH(:,i300),iS),'linewidth',2); xlabel('Years'); ylabel('RH amplification'); hl = legend('800 mb','600 mb','300 mb','fontsize',8,'location','best'); plotaxis2;
figure(7); plot(1:iNum,smooth(era5ampWV(:,i800),iS),1:iNum,smooth(era5ampWV(:,i600),iS),1:iNum,smooth(era5ampWV(:,i300),iS),'linewidth',2); xlabel('Years'); ylabel('WV amplification'); hl = legend('800 mb','600 mb','300 mb','fontsize',8,'location','best'); plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');
iAorOorL = input('Enter region : ');
if length(iAorOorL) == 0 
  iAorOorL = 0;
end

%% this is bascially copied from AIRS_gridded_STM_May2021_trendsonlyCLR/plot_driver_gather_gridded_retrieval_results. but with TWO modifications
clear maskLF
maskLF = zeros(1,4608);
maskLF = nan(1,4608);    %% MODIFICATION 1
if size(landfrac) ~= size(XX)
  XX = XX';
  YY = YY';
end
if iAorOorL == -7
  maskLF(abs(YY) > 60) = 1;
elseif iAorOorL == -6
  maskLF(abs(YY) > 30 & abs(YY) <= 60) = 1;
elseif iAorOorL == -5
  maskLF(abs(YY) < 30) = 1;
elseif iAorOorL == 0
  maskLF = ones(1,4608);
elseif iAorOorL == -1
  maskLF(landfrac == 1) = 1;
elseif iAorOorL == +1
  maskLF(landfrac == 0) = 1;
elseif iAorOorL == -2
  maskLF(landfrac == 1 & abs(YY) <= 30) = 1;
elseif iAorOorL == +2
  maskLF(landfrac == 0 & abs(YY) <= 30) = 1;
elseif iAorOorL == -3
  maskLF(landfrac == 1 & abs(YY) > 30 & abs(YY) <= 60) = 1;
elseif iAorOorL == +3
  maskLF(landfrac == 0 & abs(YY) > 30 & abs(YY) <= 60) = 1;
elseif iAorOorL == -4
  maskLF(landfrac == 1 & abs(YY) > 60) = 1;
elseif iAorOorL == +4
  maskLF(landfrac == 0 & abs(YY) > 60) = 1;
end
maskLFmatr = reshape(maskLF,72,64)';
mask = find(maskLF == 1);
maskLF100 = ones(100,1)*maskLF;  %% MODIFICATION 2

clear junk
  junk = era5;
  junk.pavg = pavg;
     junk.stemprate = junk.stemprate .* maskLF;
     junk.ptemprate = junk.ptemprate .* maskLF100;
     junk.waterrate = junk.waterrate .* maskLF100;
     junk.RHrate = junk.RHrate .* maskLF100;
     era5amp = atmospheric_amplification(junk);
  junkERA5 = junk;

