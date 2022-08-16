function nwpChoice   = getdata_NWP(iNWP,iNorD,iAorOrL,iNumYears);

addpath /asl/matlib/maps/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/PLOTTER

%   iNWP = 1   ERA
%          2   MERRA2
%          5   ERA5 (default)

if nargin == 0
  iNWP = 5; %% ERA5
  iNorD = 1;
  iAorOrL = 0;
  iNumYears  = 19;
elseif nargin == 1
  iNorD = 1;
  iAorOrL = 0;
  iNumYears  = 19;
elseif nargin == 2
  iAorOrL = 0;
  iNumYears  = 19;
elseif nargin == 3
  iNumYears  = 19;
end

if length(intersect(iNWP,[1 2 5])) == 0
  error('Need iNWP = 1,2,5 (ERA-I, MERRA2, ERA5)');
end

%iFirstTime = -1;
%if ~exist('pavgLAY')
  iFirstTime = +1;
  pjunk = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = pjunk(1:100)-pjunk(2:101);
  pjunkD = log(pjunk(1:100)./pjunk(2:101));
  pavgLAY = flipud(pjunkN./pjunkD)*ones(1,64);
%end
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
  %iAorOorL = input('Enter (-1) land (0,default) both (+1) ocean trends : ');
  iAorOorL = 0;
  if length(iAorOorL) == 0 
    iAorOorL = 0;
  end
end

if ~exist('maskLF')
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

% iNumYears = 18;
% iNumYears = 19;
% iNumYears = input('ERA-I always 17 years (2002/09 -2019/08) --- Enter Number of Years for ERA5/AIRSL3 (18 or 19, or 14 for AMPI/CMIP overlap) [19 = default] : ');
% if length(iNumYears) == 0
%   iNumYears = 19;
% end
% iNumYears = 19;

%% airsL3 : 'native' = 180 bins from L3, 'zonal' = 40 equal area latbins, [] = 64x72
if iNorD > 0
  strNorD = 'NIGHT';
  if iNWP == 5
    %nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat');
    if iNumYears == 19
      nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_desc.mat');
    elseif iNumYears == 12
      nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2014_08_trends_desc.mat');
    end
    strChoice  = 'ERA5';
  elseif iNWP == 1
    nwpChoice   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA_atm_data_2002_09_to_2019_08_16day_trends_desc.mat');
    strChoice  = 'ERA-I';
  elseif iNWP == 2
    if iNumYears == 19
      nwpChoice = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/MERRA2_atm_data_2002_09_to_2021_08_trends_desc.mat');
    elseif iNumYears == 12
      nwpChoice = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/MERRA2_atm_data_2002_09_to_2014_08_trends_desc.mat');
    end
    strChoice  = 'MERRA2';
  else
    error('unknown NWP choice')
  end
  %cmip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08_trends.mat');
  %amip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/AMIP6_atm_data_2002_09_to_2014_08_trends.mat');
else
  strNorD = 'DAY';
  if iNWP == 5
    nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat');
    %if iNumYears == 19
    %  nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_asc.mat'); %%% OOER may not have run this off yet
    %elseif iNumYears == 12
    %  nwpChoice  = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2014_08_trends_asc.mat'); %%% OOER may not have run this off yet
    %end
    strChoice  = 'ERA5';
  elseif iNWP == 1
    nwpChoice   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA_atm_data_2002_09_to_2019_08_16day_trends_asc.mat');
    strChoice  = 'ERA-I';
  elseif iNWP == 2
    if iNumYears == 19
      nwpChoice = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/MERRA2_atm_data_2002_09_to_2021_08_trends_desc.mat');
    elseif iNumYears == 12
      nwpChoice = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/MERRA2_atm_data_2002_09_to_2014_08_trends_desc.mat');
    end
    strChoice  = 'MERRA2';
  else
    error('unknown NWP choice')
  end
  %cmip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08_trends.mat');
  %amip6 = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/AMIP6_atm_data_2002_09_to_2014_08_trends.mat');
end

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

load llsmap5;
iFig = 0;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(nwpChoice.trend_stemp,72,64)',1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ' strChoice ' K/yr']); 

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
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(nwpChoice.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ' strChoice ' K/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(nwpChoice.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ' strChoice ' percent/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(nwpChoice.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ' strChoice ' 1/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
