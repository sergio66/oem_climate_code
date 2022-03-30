function airsChoice  = getdata_AIRS_L3vsCLIMCAPS(iA,iNorD,iAorOrL);

addpath /asl/matlib/maps/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/PLOTTER

if nargin == 0
  iA = 1;
  iNorD = 1;
  iAorOrL = 0;
elseif nargin == 1
  iNorD = 1;
  iAorOrL = 0;
elseif nargin == 2
  iAorOrL = 0;
end

if iA == 1
  strXL3 = 'AIRS L3';
else
  strXL3 = 'CLIMCAPS';
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

%iNumYears = 18;
iNumYears = 19;
%iNumYears = input('ERA-I always 17 years (2002/09 -2019/08) --- Enter Number of Years for ERA5/AIRSL3 (18 or 19) [19 = default] : ');
%if length(iNumYears) == 0
%  iNumYears = 19;
%end

%% airsL3 : 'native' = 180 bins from L3, 'zonal' = 40 equal area latbins, [] = 64x72
if iNorD > 0
  strNorD = 'NIGHT';
  if iA == 1
    if iNumYears == 19
      %airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
      %airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
      %airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
      airsChoice    = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2021_19yr_desc.mat');
    end
  else
    airsChoice       = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2021_19yr_desc.mat');
  end
else
  strNorD = 'DAY';
  if iA == 1
    if iNumYears == 19  
      %airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
      %airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
      airsChoice       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
    end
  else
    airsChoice       = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2021_19yr_asc.mat');
  end
end

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

if iA == 1
  iT = 24;
  iW = 12;
else
  iT = 100;
  iW = 66;
  airsChoice.Tlevs = airsChoice.Tlevs/100;  %% Pa --> mb
  airsChoice.Qlevs = airsChoice.Qlevs/100;  %% Pa --> mb
  airsChoice.thestats64x72.RHrate = airsChoice.thestats64x72.RHrate * 100; 
  airsChoice.thestats64x72.RHSurfrate = airsChoice.thestats64x72.RHSurfrate * 100;
end

load llsmap5;
iFig = 0;
%iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(airsChoice.thestats64x72.RHSurfrate',1), [-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.5 +0.5]);   title([strNorD ' RHsurf d/dt ' strXL3 ' /yr']);   
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(airsChoice.thestats64x72.stemprate',1), [-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ' strXL3 ' K/yr']);

%if ~exist('pavgLAY')
%  boo = load('/home/sergio/MATLABCODE/airslevels.dat');
%  pjunkN = boo(1:100)-boo(2:101);
%  pjunkD = log(boo(1:100)./boo(2:101));
%  pavgLAY = pjunkN./pjunkD;
%  pavgLAY = flipud(pavgLAY);
%  pavgLAY = pavgLAY*ones(1,iCntr);
%end

if ~exist('rlat')
  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
end

maskLFmatr = ones(64,72);
boo = zeros(72,64,iT); for ijunk = 1 : iT; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsChoice.thestats64x72.ptemprate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsChoice.Tlevs,smoothn(junk(1:iT,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ' strXL3 ' K/yr']);

boo = zeros(72,64,iW); for ijunk = 1 : iW; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsChoice.thestats64x72.RHrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsChoice.Qlevs,smoothn(junk(1:iW,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ' strXL3 ' percent/yr']);

boo = zeros(72,64,iW); for ijunk = 1 : iW; boo(:,:,ijunk) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsChoice.thestats64x72.waterrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsChoice.Qlevs,smoothn(junk(1:iW,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ' strXL3 ' 1/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
