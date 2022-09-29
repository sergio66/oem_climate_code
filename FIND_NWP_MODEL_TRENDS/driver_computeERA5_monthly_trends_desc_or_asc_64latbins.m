%% see driver_computeERA5_monthly_trends_desc_or_asc.m
%% monthly, 18 years x 12 months/year = 216
%% monthly, 19 years x 12 months/year = 228

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/IDL_WV_ROUTINES/atmos_phys/MATLAB/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

load('llsmap5.mat');

if ~exist('iDorA')
  iDorA = -1; %% asc
  iDorA = +1; %% desc
end

clear iaFound

timespan = 18;  %% AIRS
timespan = 19;  %% AIRS
timespan = 20;  %% AIRS
timespan = 12;  %% AMPI6
timespan = 07;  %% CrIS NSR
timespan = 06;  %% OCO2
fprintf(1,'timespan = %2i years \n',timespan)

if timespan == 06
  %% OCO2
  savestr_version = 'Jan2015_Dec2021_OCO2_06yr';
  StartY = 2015; StartYM = 1;    %% start 01/2015
  StopY  = 2021; StopYM  = 12;   %% stop  12/2021  
  savestr_version = 'Sep2014_Aug2021_OCO2_06yr';
  StartY = 2014; StartYM = 09;   %% start 01/2015
  StopY  = 2021; StopYM  = 08;   %% stop  12/2021  
elseif timespan == 07
  %% CrIS NSR
  savestr_version = 'May2012_Apr2019_07yr';
  StartY = 2012; StartYM = 5;   %% start 05/2012
  StopY  = 2019; StopYM  = 4;   %% stop  04/2019  
elseif timespan == 12
  %% AMIP6
  savestr_version = 'Sept2002_Aug2014_12yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2014; StopYM  = 8;   %% stop  08/2014
%%
elseif timespan == 16
  %% HUH
  savestr_version = 'Sept2017_15yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2018; StopYM  = 8;   %% stop  08/2017  
elseif timespan == 18
  %% AIRS
  savestr_version = 'Sept2002_Aug2020_18yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2020; StopYM  = 8;   %% stop  08/2020
elseif timespan == 19
  %% AIRS
  savestr_version = 'Sept2002_Jul2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 7;   %% stop  08/2021  
  savestr_version = 'Sept2002_Aug2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 8;   %% stop  08/2021  
else
  error('huh check timespan')
end

iNumYears = timespan;
iaMax = iNumYears*12;

comment = 'see computeERA5_trends.m';
comment = 'see driver_computeERA5_monthly_trends_desc_or_asc.m';
if iNumYears == 18 
  if iDorA > 0
    load ERA5_atm_data_2002_09_to_2020_08_desc.mat
  else
    load ERA5_atm_data_2002_09_to_2020_08_asc.mat
  end
elseif iNumYears == 19 | iNumYears == 12 | iNumYears == 07 | iNumYears == 06
  if iDorA > 0
    load ERA5_atm_data_2002_09_to_2021_08_desc.mat
  else
    load ERA5_atm_data_2002_09_to_2021_08_asc.mat
  end
else
  error('can only handle 06,07,12,18,19 years right now')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = flipud(pN./pD);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);
iaIndexUse = zeros(size(dayOFtime));

dayOFtimeS = change2days(StartY,StartYM,01,2002);
if StopYM == 2
  dayOFtimeE = change2days(StopY,StopYM,  28,2002);
elseif length(intersect(StopYM,[1 3 5 78 10 12])) == 1
  dayOFtimeE = change2days(StopY,StopYM,  31,2002);
else
  dayOFtimeE = change2days(StopY,StopYM,  30,2002);
end
usetime = find(dayOFtime >= dayOFtimeS & dayOFtime < dayOFtimeE);
iaIndexUse(usetime) = 1;

iAorOorL = -1; %% land
iAorOorL = +1; %% ocean 
iAorOorL = 0;  %% all 
clear trend*

computeERA5_surface_trends_64latbins
computeERA5_atmos_trends_64latbins

if iNumYears == 06
  if iDorA > 0 & iAorOorL == 0
    save ERA5_atm_data_2014_09_to_2021_08_trends_desc_64latbins.mat comment trend*
  else
    error(';ksjg')
  end
elseif iNumYears == 07
  if iDorA > 0 & iAorOorL == 0
    save ERA5_atm_data_2012_05_to_2019_04_trends_desc_64latbins.mat comment trend*
  elseif iDorA < 0 & iAorOorL == 0
    save ERA5_atm_data_2012_05_to_2019_04_trends_asc_64latbins.mat comment trend*
  elseif iDorA > 0 & iAorOorL == 1
    save ERA5_atm_data_2012_05_to_2019_04_trends_desc_64latbins_ocean.mat comment trend*
  elseif iDorA < 0 & iAorOorL == 1
    save ERA5_atm_data_2012_05_to_2019_04_trends_asc_64latbins_ocean.mat comment trend*
  elseif iDorA > 0 & iAorOorL == -1
    save ERA5_atm_data_2012_05_to_2019_04_trends_desc_64latbins_land.mat comment trend*
  elseif iDorA < 0 & iAorOorL == -1
    save ERA5_atm_data_2012_05_to_2019_04_trends_asc_64latbins_land.mat comment trend*
  end
elseif iNumYears == 12
  if iDorA > 0
    save ERA5_atm_data_2002_09_to_2014_08_trends_desc_64latbins.mat comment trend*
  else
    save ERA5_atm_data_2002_09_to_2014_08_trends_asc_64latbins.mat comment trend*
  end
elseif iNumYears == 18
  if iDorA > 0
    save ERA5_atm_data_2002_09_to_2020_08_trends_desc_64latbins.mat comment trend*
  else
    save ERA5_atm_data_2002_09_to_2020_08_trends_asc_64latbins.mat comment trend*
  end
elseif iNumYears == 19
  if iDorA > 0
    save ERA5_atm_data_2002_09_to_2021_08_trends_desc_64latbins.mat comment trend*
  else
    save ERA5_atm_data_2002_09_to_2021_08_trends_asc_64latbins.mat comment trend*
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ../AIRS_gridded_STM_May2021_trendsonlyCLR/statistical_significance_ludescherPNAS_lag1_functions.m
% also see strow_override_defaults_latbins_AIRS_fewlays.m
% if driver.i16daytimestep < 0
%   nc_cor = nc_rates(driver);
%   driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;    %% THIS IS AN ARRAY
%   %driver.rateset.unc_rates = driver.rateset.unc_rates *sqrt(1/8); %% this accounts for counting ....
% end
%
% where from nc_rates.m
%    nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
%    nc = sqrt(nc);
%    nc_cor = real(nc');


plays100 = load('/home/sergio/MATLABCODE/airslevels.dat');
plays100 = plevs2plays(plays100); plays100 = flipud(plays100(1:100));

figure(1); pcolor(rlat,plays100,trend64_ptemp); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dT(z)/dt')
figure(2); pcolor(rlat,plays100,trend64_ptemp_err); colormap(llsmap5); caxis([0 +1]*0.05); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dT(z)/dt unc'); colormap(jet)
figure(3); pcolor(rlat,plays100,trend64_ptemp_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dT(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_ptemp_lag)./(1-trend64_ptemp_lag)));
figure(4); pcolor(rlat,plays100,trend64_ptemp_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*0.50); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dT(z)/dt unc.*corr'); colormap(jet)

figure(1); pcolor(rlat,plays100,trend64_gas_1); colormap(llsmap5); caxis([-1 +1]*0.01); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dWVfrac(z)/dt')
figure(2); pcolor(rlat,plays100,trend64_gas_1_err); colormap(llsmap5); caxis([0 +1]*0.001); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dWVfrac(z)/dt unc'); colormap(jet)
figure(3); pcolor(rlat,plays100,trend64_gas_1_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dWVfrac(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_gas_1_lag)./(1-trend64_gas_1_lag)));
figure(4); pcolor(rlat,plays100,trend64_gas_1_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*0.010); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dWVfrac(z)/dt unc.*corr'); colormap(jet)

figure(1); pcolor(rlat,plays100,trend64_RH); colormap(llsmap5); caxis([-1 +1]*0.5); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dRH(z)/dt')
figure(2); pcolor(rlat,plays100,trend64_RH_err); colormap(llsmap5); caxis([0 +1]*0.1); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dRH(z)/dt unc'); colormap(jet)
figure(3); pcolor(rlat,plays100,trend64_RH_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dRH(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_RH_lag)./(1-trend64_RH_lag)));
figure(4); pcolor(rlat,plays100,trend64_RH_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*1); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([1 1000]); title('dRH(z)/dt unc.*corr'); colormap(jet)

figure(1); pcolor(rlat,plays100,trend64_ptemp); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('dT(z)/dt')
figure(2); pcolor(rlat,plays100,trend64_gas_1); colormap(llsmap5); caxis([-1 +1]*0.01); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('dWVfrac(z)/dt')
figure(3); pcolor(rlat,plays100,trend64_RH); colormap(llsmap5); caxis([-1 +1]*0.5); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('dRH(z)/dt')

