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
iaMax = 18*12; %% 18 year
iaMax = 19*12; %% 19 year

iNumYears = 18; 
iNumYears = 19; 
iaMax = iNumYears*12;

comment = 'see computeERA5_trends.m';
comment = 'see driver_computeERA5_monthly_trends_desc_or_asc.m';
if iNumYears == 18
  if iDorA > 0
    load ERA5_atm_data_2002_09_to_2020_08_desc.mat
  else
    load ERA5_atm_data_2002_09_to_2020_08_asc.mat
  end
elseif iNumYears == 19
  if iDorA > 0
    load ERA5_atm_data_2002_09_to_2021_08_desc.mat
  else
    load ERA5_atm_data_2002_09_to_2021_08_asc.mat
  end
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

computeERA5_surface_trends_64latbins
computeERA5_atmos_trends_64latbins

if iNumYears == 18
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


figure(1); pcolor(trend64_ptemp); colormap(llsmap5); caxis([-1 +1]*0.15); colorbar; shading interp; set(gca,'ydir','reverse'); title('dT(z)/dt')
figure(2); pcolor(trend64_ptemp_err); colormap(llsmap5); caxis([0 +1]*0.05); colorbar; shading interp; set(gca,'ydir','reverse'); title('dT(z)/dt unc'); colormap(jet)
figure(3); pcolor(trend64_ptemp_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); title('dT(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_ptemp_lag)./(1-trend64_ptemp_lag)));
figure(4); pcolor(trend64_ptemp_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*0.50); colorbar; shading interp; set(gca,'ydir','reverse'); title('dT(z)/dt unc.*corr'); colormap(jet)

figure(1); pcolor(trend64_gas_1); colormap(llsmap5); caxis([-1 +1]*0.01); colorbar; shading interp; set(gca,'ydir','reverse'); title('dWVfrac(z)/dt')
figure(2); pcolor(trend64_gas_1_err); colormap(llsmap5); caxis([0 +1]*0.001); colorbar; shading interp; set(gca,'ydir','reverse'); title('dWVfrac(z)/dt unc'); colormap(jet)
figure(3); pcolor(trend64_gas_1_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); title('dWVfrac(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_gas_1_lag)./(1-trend64_gas_1_lag)));
figure(4); pcolor(trend64_gas_1_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*0.010); colorbar; shading interp; set(gca,'ydir','reverse'); title('dWVfrac(z)/dt unc.*corr'); colormap(jet)

figure(1); pcolor(trend64_RH); colormap(llsmap5); caxis([-1 +1]*0.5); colorbar; shading interp; set(gca,'ydir','reverse'); title('dRH(z)/dt')
figure(2); pcolor(trend64_RH_err); colormap(llsmap5); caxis([0 +1]*0.1); colorbar; shading interp; set(gca,'ydir','reverse'); title('dRH(z)/dt unc'); colormap(jet)
figure(3); pcolor(trend64_RH_lag); colormap(llsmap5); caxis([0.75 +1]*1.1); colorbar; shading interp; set(gca,'ydir','reverse'); title('dRH(z)/dt lag1'); colormap(jet)
nc_cor = real(sqrt((1+trend64_RH_lag)./(1-trend64_RH_lag)));
figure(4); pcolor(trend64_RH_err.*nc_cor); colormap(llsmap5); caxis([0 +1]*1); colorbar; shading interp; set(gca,'ydir','reverse'); title('dRH(z)/dt unc.*corr'); colormap(jet)

