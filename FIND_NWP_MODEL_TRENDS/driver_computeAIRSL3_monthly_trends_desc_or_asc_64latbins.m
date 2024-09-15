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

error('need to translate AIRSL3 to ERA5 fiels!! ugh!!!')
error('need to translate AIRSL3 to ERA5 fiels!! ugh!!!')
error('need to translate AIRSL3 to ERA5 fiels!! ugh!!!')

load('llsmap5.mat');

if ~exist('iDorA')
  iDorA = -1; %% asc
  iDorA = +1; %% desc
end

clear iaFound

timespan = 18;  %% AIRS
timespan = 19;  %% AIRS
timespan = 12;  %% AMPI6
timespan = 07;  %% CrIS NSR
timespan = 06;  %% OCO2
timespan = 20;  %% AIRS
timespan = 22;  %% AIRS
fprintf(1,'timespan = %2i years \n',timespan)

if timespan == 22
  %% AIRS
  savestr_version = 'Sept2002_Aug2024_22yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2024; StopYM  = 8;   %% stop  08/2024
else
  error('huh check timespan')
end

iNumYears = timespan;
iaMax = iNumYears*12;

comment = 'see driver_computeAIRSL3_monthly_trends_desc_or_asc_64latbins.m';
if iNumYears == 20 | iNumYears == 22
  if iDorA > 0
    load /asl/s1/sergio/AIRS_L3/airs_L3v7_extra_Sept2002_Aug2024_22yr.mat
  else
    load /asl/s1/sergio/AIRS_L3/airs_L3v7_extra_Sept2002_Aug2024_22yr.mat
  end
else
  error('can only handle 20 or 22 years right now')
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
if ~exist('dd')
 dd = 15 * ones(size(yy));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
all.yy = yy; 
all.mm = mm; 
all.dd = dd; 
error('need to translate AIRSL3 to ERA5 fiels!! ugh!!!')

%%%%%%%%%%%%%%%%%%%%%%%%%

dayOFtime = change2days(yy,mm,dd,2002);
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

computeERA5_surface_anoms_64latbins
computeERA5_atmos_anoms_64latbins

yymm = 2002 + dayOFtime/365.25;
figure(1); clf
plot(yymm,smooth(nanmean(anom64_olr_clr,1),5),'b',yymm,smooth(nanmean(anom64_olr,1),5),'k',yymm,smooth(nanmean(anom64_stemp,1),5),'r','linewidth',2); 
  plotaxis2; hl = legend('OLR CLR','OLR','SKT'); xlabel('Time');

figure(2); clf
plot(rlat,smooth(nanmean(anom64_olr_clr,2),5),'b',rlat,smooth(nanmean(anom64_olr,2),5),'k',rlat,smooth(nanmean(anom64_stemp,2),5),'r','linewidth',2); 
  plotaxis2; hl = legend('OLR CLR','OLR','SKT'); xlabel('Latitude');

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = flipud(meanvaluebin(plevs));

figure(3); clf
pcolor(yymm,1:100,squeeze(nanmean(anom64_ptemp,2))); shading interp; colorbar; caxis([-1 +1]*2); set(gca,'ydir','reverse'); title('mean T anomaly'); xlabel('Time'); ylabel('Layer')
pcolor(yymm,plays,squeeze(nanmean(anom64_ptemp,2))); shading interp; colorbar; caxis([-1 +1]*2); set(gca,'ydir','reverse'); title('mean T anomaly'); xlabel('Time'); ylabel('P [mb]')
  set(gca,'yscale','log'); ylim([10 1000])

figure(4); clf
pcolor(yymm,1:100,squeeze(nanmean(anom64_gas_1,2))); shading interp; colorbar; caxis([-1 +1]*2/10); set(gca,'ydir','reverse'); title('mean fracWV anomaly'); xlabel('Time'); ylabel('Layer')
pcolor(yymm,plays,squeeze(nanmean(anom64_gas_1,2))); shading interp; colorbar; caxis([-1 +1]*2/10); set(gca,'ydir','reverse'); title('mean fracWV anomaly'); xlabel('Time'); ylabel('P [mb]')
  set(gca,'yscale','log'); ylim([10 1000])

figure(5); clf
pcolor(yymm,1:100,squeeze(nanmean(anom64_RH,2))); shading interp; colorbar; caxis([-1 +1]*5); set(gca,'ydir','reverse'); title('mean RH anomaly'); xlabel('Time'); ylabel('Layer')
pcolor(yymm,plays,squeeze(nanmean(anom64_RH,2))); shading interp; colorbar; caxis([-1 +1]*5); set(gca,'ydir','reverse'); title('mean RH anomaly'); xlabel('Time'); ylabel('P [mb]')
  set(gca,'yscale','log'); ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the results

if iNumYears == 06
  if iDorA > 0 & iAorOorL == 0
    save AIRSL3_atm_data_2014_09_to_2021_08_trends_desc_64latbins.mat comment trend* anom*
  else
    error(';ksjg')
  end
elseif iNumYears == 07
  if iDorA > 0 & iAorOorL == 0
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_desc_64latbins.mat comment trend* anom*
  elseif iDorA < 0 & iAorOorL == 0
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_asc_64latbins.mat comment trend* anom*
  elseif iDorA > 0 & iAorOorL == 1
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_desc_64latbins_ocean.mat comment trend* anom*
  elseif iDorA < 0 & iAorOorL == 1
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_asc_64latbins_ocean.mat comment trend* anom*
  elseif iDorA > 0 & iAorOorL == -1
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_desc_64latbins_land.mat comment trend* anom*
  elseif iDorA < 0 & iAorOorL == -1
    save AIRSL3_atm_data_2012_05_to_2019_04_trends_asc_64latbins_land.mat comment trend* anom*
  end
end
return

if iNumYears == 12
  if iDorA > 0
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2014_08_trends_desc_64latbins.mat'];
  else
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2014_08_trends_asc_64latbins.mat'];
  end
elseif iNumYears == 18
  if iDorA > 0
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2020_08_trends_desc_64latbins.mat'];
  else
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2020_08_trends_asc_64latbins.mat'];
  end
elseif iNumYears == 19
  if iDorA > 0
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2021_08_trends_desc_64latbins.mat'];
  else
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2021_08_trends_asc_64latbins.mat'];
  end
elseif iNumYears == 20
  if iDorA > 0
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2022_08_trends_desc_64latbins.mat'];
  else
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2022_08_trends_asc_64latbins.mat'];
  end
elseif iNumYears == 22
  if iDorA > 0
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2024_08_trends_desc_64latbins.mat'];
  else
    fnamesave_anom_trend = ['AIRSL3_atm_data_2002_09_to_2024_08_trends_asc_64latbins.mat'];
  end
end
saver = ['save ' fnamesave_anom_trend ' comment trend* anom* dayOFtime StartY StopY StartYM StopYM'];
eval(saver);

