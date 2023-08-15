function airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL,iNumYears);

%% function airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL);   iA = +1 (AIRS L3) or -1 CLIMCAPS

addpath /asl/matlib/maps/
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/PLOTTER

if nargin == 0
  iA = 1;
  iNorD = 1;
  iAorOrL = 0;
  iNumYears = 19;
elseif nargin == 1
  iNorD = 1;
  iAorOrL = 0;
  iNumYears = 19;
elseif nargin == 2
  iAorOrL = 0;
  iNumYears = 19;
elseif nargin == 3
  iNumYears = 19;
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
%iNumYears = 19;
%iNumYears = input('ERA-I always 17 years (2002/09 -2019/08) --- Enter Number of Years for ERA5/AIRSL3 (18 or 19) [19 = default] : ');
%if length(iNumYears) == 0
%  iNumYears = 19;
%end

iNumYears0 = iNumYears;
zoo = [05 10 12 15 19 20];
moo = abs(iNumYears - zoo);
moo = find(moo == min(moo));
iNumYears = zoo(moo);
if iNumYears ~= iNumYears0
  fprintf(1,'getdata_AIRSL3vsCLIMCAPSL3.m : You want trends for %2i years but can only find closest = %2i years \n',iNumYears0,iNumYears);
end

%% airsL3 : 'native' = 180 bins from L3, 'zonal' = 40 equal area latbins, [] = 64x72
if iNorD > 0
  strNorD = 'NIGHT';
  if iA == 1
    if length(intersect(iNumYears,[05 10 12 15 18 19 20])) == 1
      fAIRS = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
    else
      iaJunk = [5 10 15 20 12 18 19];
      moo = abs(iaJunk-iNumYears);
      moo = find(moo == min(moo),1);    
      iNumYearsPowWow = 2002+iaJunk(moo);
      fAIRS = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug' num2str(2002+iNumYearsPowWow) '_' num2str(iNumYearsPowWow) 'yr_desc.mat'];
      disp('AIRS L3 : needs 12,18,19 or 05,10,15,20 years, subbing in closest year')
    end
    % if iNumYears == 19
    %   %airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    %   %airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    %   %airsL3       = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_desc.mat');
    %   fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_desc.mat';
    %   fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_desc.mat';
    %   fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2021_19yr_desc.mat';
    % elseif iNumYears == 20
    %   fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat';
    % elseif iNumYears == 12
    %   fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2014_12yr_desc.mat';
    % else
    %   iNumYears
    %   error('need 12, 19, 20 years')
    % end
  else
    if length(intersect(iNumYears,[05 10 12 15 19 20])) == 1
      fAIRS = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
    else
      iaJunk = [5 10 15 20 12 18 19];
      moo = abs(iaJunk-iNumYears);
      moo = find(moo == min(moo),1);    
      iNumYearsPowWow = 2002+iaJunk(moo);
      fAIRS = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug' num2str(2002+iNumYearsPowWow) '_' num2str(iNumYearsPowWow) 'yr_desc.mat'];
      disp('CLIMCAPS L3 : need 05 10 12 15 19 or 20 years, subbing in closest year')
    end    
  end
else
  strNorD = 'DAY';
  if iA == 1
    if iNumYears == 19  
      %airsL3native = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_native_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
      %airsL3zonal  = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_Sept2002_Jul2021_19yr_asc.mat');
      fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Jul2021_19yr_asc.mat';
    elseif iNumYears == 12
      fAIRS = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2014_12yr_asc.mat';
    else
      error('need 12 or 19 years')
    end
  else
    if iNumYears == 19  
      fAIRS = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2021_19yr_asc.mat';
    elseif iNumYears == 20
      fAIRS = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2022_20yr_asc.mat';
    elseif iNumYears == 12
     fAIRS = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_Sept2002_Aug2014_12yr_asc.mat';
    else
      error('need 12 or 19 years')
    end
  end
end

fprintf(1,' getdata_AIRSL3vsCLIMCAPSL3.m : iA = %2i iNumYears = %2i : fAIRS = %s \n',iA,iNumYears,fAIRS);
airsChoice       = load(fAIRS);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

if iA == 1
  iT = 24;
  iW = 12;
  %% pressure units already in mb, RH units already in percent
else
  iT = 100;
  iW = 66;
  airsChoice.Tlevs = airsChoice.Tlevs/100;  %% Pa --> mb
  airsChoice.Qlevs = airsChoice.Qlevs/100;  %% Pa --> mb
  airsChoice.thestats64x72.RHrate = airsChoice.thestats64x72.RHrate * 100;           %% fraction --> percent
  airsChoice.thestats64x72.RHSurfrate = airsChoice.thestats64x72.RHSurfrate * 100;   %% fraction --> percent
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
