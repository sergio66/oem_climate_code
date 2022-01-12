addpath /asl/matlib/maps/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/PLOTTER

iFirstTime = -1;
if ~exist('pavgLAY')
  iFirstTime = +1;
  pjunk = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = pjunk(1:100)-pjunk(2:101);
  pjunkD = log(pjunk(1:100)./pjunk(2:101));
  pavgLAY = flipud(pjunkN./pjunkD)*ones(1,64);
end
[mmX,nnX] = size(pavgLAY);
if nnX == 4608
  iCntr = 3000;
elseif nnX == 64
  iCntr = 32;
end

if ~exist('iNorD')
  iNorD = input('Enter (+1,DEFAULT) night (-1) day  trends : ');
  if length(iNorD) == 0
    iNorD = 1;
  end
end

if ~exist('iAorOorL')
  iAorOorL = input('Enter (-1) land (0,default) both (+1) ocean trends : ');
  if length(iAorOorL) == 0 
    iAorOorL = 0;
  end

  clear maskLF
  maskLF = zeros(1,4608);
  if iAorOorL == 0
    maskLF = ones(1,4608);
  elseif iAorOorL == -1
    maskLF(landfrac == 1) = 1;
  elseif iAorOorL == +1
    maskLF(landfrac == 0) = 1;
  end
  maskLFmatr = reshape(maskLF,72,64)';

end

iNumYears = 18;
iNumYears = 19;
iNumYears = input('ERA-I always 17 years (2002/09 -2019/08) --- Enter Number of Years for ERA5/AIRSL3 (18 or 19) [19 = default] : ');
if length(iNumYears) == 0
  iNumYears = 19;
end

%% airsL3 : 'native' = 180 bins from L3, 'zonal' = 40 equal area latbins, [] = 64x72
if iNorD > 0
  strNorD = 'NIGHT';
  if iNumYears == 18
    airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Aug2020_18yr_desc.mat');
    airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Aug2020_18yr_desc.mat');
    airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2020_18yr_desc.mat');
    %era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2019_08_trends_desc.matWRONGYEAREND');
    era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2020_08_trends_desc.mat');
  elseif iNumYears == 19
    airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_desc.mat');
  end
  era   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA_atm_data_2002_09_to_2019_08_16day_trends_desc.mat');
  cmip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08_trends.mat');
else
  strNorD = 'DAY';
  if iNumYears == 18
    airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Aug2020_18yr_asc.mat');
    airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Aug2020_18yr_asc.mat');
    airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2020_18yr_asc.mat');
    %era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2019_08_trends_asc.matWRONGYEAREND');
    %era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2020_08_trends_asc.mat');
    era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_asc.mat');
  elseif iNumYears == 19  
    airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
    airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
    airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
    era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat');
    %era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_asc.mat');   %%% OOER may not have run this off yet
  end
  era    = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA_atm_data_2002_09_to_2019_08_16day_trends_asc.mat');
  cmip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08_trends.mat');
end

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

load llsmap5;
iFig = 40;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(airsL3.thestats64x72.RHSurfrate',1), [-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.5 +0.5]);   title([strNorD ' RHsurf d/dt AIRS L3 /yr']);   
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(era.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ERA K/yr']);   
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(era5.trend_stemp,72,64)',1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ERA5 K/yr']); 
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(airsL3.thestats64x72.stemprate',1), [-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt AIRSL3 K/yr']);
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(cmip6.trend_stemp,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.15 +0.15]); title(['D/N stemp d/dt CMIP6 K/yr']);

if ~exist('pavgLAY')
  boo = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = boo(1:100)-boo(2:101);
  pjunkD = log(boo(1:100)./boo(2:101));
  pavgLAY = pjunkN./pjunkD;
  pavgLAY = flipud(pavgLAY);
  pavgLAY = pavgLAY*ones(1,iCntr);
end

if ~exist('rlat')
  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
end

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ERA K/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ERA5 K/yr']);

boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.ptemprate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Tlevs,smoothn(junk(1:24,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt AIRS L3 K/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* reshape(cmip6.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title(['D/N dT/dt CMIP6  K/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ERA percent/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ERA5 percent/yr']);

boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.RHrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt AIRS L3 percent/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(cmip6.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title(['D/N dRH/dt CMIP6 percent/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ERA 1/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ERA5 1/yr']);

boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.waterrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt AIRS L3 1/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(cmip6.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title(['D/N d(fracWV)/dt CMIP6 1/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iFirstTime > 0

  figure(1);
  pcolor(rlat,airsL3.Qlevs,squeeze(nanmean(airsL3.thestats64x72.RHrate,1))'); shading interp; 
  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = boo.*airsL3.thestats64x72.RHrate; junk = squeeze(nanmean(junk,1));
  pcolor(rlat,airsL3.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dRH/dt AIRS L3 \newline from zonally averaged 64x72')
  
  figure(2);
  junk40 = equal_area_spherical_bands(20);
  junk40 = 0.5*(junk40(1:end-1)+junk40(2:end));
  pcolor(junk40,airsL3zonal.Qlevs,airsL3zonal.thestats.RHrate'); shading interp; 
  boo = zeros(40,12); for ijunk = 1 : 12; boo(:,ijunk) = 1; end
  junk = boo.*airsL3zonal.thestats.RHrate; 
  pcolor(junk40,airsL3zonal.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dRH/dt AIRS L3 \newline from zonally EQ AREA 40 bins')
  
  figure(3);
  junk180 = (-90:1:89)+0.5;
  pcolor(junk180,airsL3native.Qlevs,airsL3native.thestats.RHrate'); shading interp; 
  boo = zeros(180,12); for ijunk = 1 : 12; boo(:,ijunk) = 1; end
  junk = boo.*airsL3native.thestats.RHrate; 
  pcolor(junk180,airsL3native.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dRH/dt AIRS L3 \newline from zonally NATIVE 180 bins')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(4);
  pcolor(rlat,airsL3.Tlevs,squeeze(nanmean(airsL3.thestats64x72.ptemprate,1))'); shading interp; 
  boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
  junk = boo.*airsL3.thestats64x72.ptemprate; junk = squeeze(nanmean(junk,1));
  pcolor(rlat,airsL3.Tlevs,smoothn(junk',1)); shading interp; 
  ylim([10 1000]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dT/dt AIRS L3 \newline from zonally averaged 64x72')
  
  figure(5);
  junk40 = equal_area_spherical_bands(20);
  junk40 = 0.5*(junk40(1:end-1)+junk40(2:end));
  pcolor(junk40,airsL3zonal.Tlevs,airsL3zonal.thestats.ptemprate'); shading interp; 
  boo = zeros(40,24); for ijunk = 1 : 24; boo(:,ijunk) = 1; end
  junk = boo.*airsL3zonal.thestats.ptemprate; 
  pcolor(junk40,airsL3zonal.Tlevs,smoothn(junk',1)); shading interp; 
  ylim([10 1000]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dT/dt AIRS L3 \newline from zonally EQ AREA 40 bins')
  
  figure(6);
  junk180 = (-90:1:89)+0.5;
  pcolor(junk180,airsL3native.Tlevs,airsL3native.thestats.ptemprate'); shading interp; 
  boo = zeros(180,24); for ijunk = 1 : 24; boo(:,ijunk) = 1; end
  junk = boo.*airsL3native.thestats.ptemprate; 
  pcolor(junk180,airsL3native.Tlevs,smoothn(junk',1)); shading interp; 
  ylim([10 1000]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]); colorbar; title('dT/dt AIRS L3 \newline from zonally NATIVE 180 bins')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(7);
  pcolor(rlat,airsL3.Qlevs,squeeze(nanmean(airsL3.thestats64x72.waterrate,1))'); shading interp; 
  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = boo.*airsL3.thestats64x72.waterrate; junk = squeeze(nanmean(junk,1));
  pcolor(rlat,airsL3.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]/10); colorbar; title('dWVfrac/dt AIRS L3 \newline from zonally averaged 64x72')
  
  figure(8);
  junk40 = equal_area_spherical_bands(20);
  junk40 = 0.5*(junk40(1:end-1)+junk40(2:end));
  pcolor(junk40,airsL3zonal.Qlevs,airsL3zonal.thestats.waterrate'); shading interp; 
  boo = zeros(40,12); for ijunk = 1 : 12; boo(:,ijunk) = 1; end
  junk = boo.*airsL3zonal.thestats.waterrate; 
  pcolor(junk40,airsL3zonal.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]/10); colorbar; title('dWVfrac/dt AIRS L3 \newline from zonally EQ AREA 40 bins')
  
  figure(9);
  junk180 = (-90:1:89)+0.5;
  pcolor(junk180,airsL3native.Qlevs,airsL3native.thestats.waterrate'); shading interp; 
  boo = zeros(180,12); for ijunk = 1 : 12; boo(:,ijunk) = 1; end
  junk = boo.*airsL3native.thestats.waterrate; 
  pcolor(junk180,airsL3native.Qlevs,smoothn(junk',1)); shading interp; 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5); caxis([-0.15 +0.15]/10); colorbar; title('dWVfrac/dt AIRS L3 \newline from zonally NATIVE 180 bins')
end  
