%% driver_compute_AIRSL3_trends_desc_or_ascNOQuestioN --> do_the_AIRSL3_trends_OQuestioN --> do_the_fits_airsL3_ratesv7_tiles --> 

%% load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc.mat   has the computed trends in thestats64x72,thestats64x72_other
%% load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc.mat         has the data you WANT TO FIT YAYAYAYAYA

addpath /asl/matlib/aslutil/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/FIND_TRENDS/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

if ~exist('save64x72_olr')
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc.mat
end

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
drlon = 5; 
rlon = -180 : drlon : +180;      %% 73 edges, so 72 bins
rlat = latB2;                    %% 65 edges, so 64 bins

iNumYears = 22;
timespan = iNumYears;
fprintf(1,'timespan = %2i years \n',timespan)

savestr_version = ['Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(timespan) 'yr'];
StartY = 2002;             StartYM = 9;   %% start 09/2002
StopY  = StartY+iNumYears; StopYM  = 8;   %% stop  08/2021  
                           StopYM  = 6;   %% stop  06/2021  

junkcnt = 0;
for yyjunk = StartY : StopY
  mmS = 1; mmE = 12;
  if yyjunk == StartY
    mmS = 9;
  elseif yyjunk == StopY
    mmE = 8;
    mmE = 6;
  end
  for mmjunk = mmS : mmE
    junkcnt = junkcnt + 1;
    yy(junkcnt) = yyjunk;
    mm(junkcnt) = mmjunk;
  end
end

zonkA = find(yy == StartY & mm == StartYM);
if length(zonkA) == 0
  zonkA = 1;
end
zonkB = find(yy == StopY  & mm == StopYM);
if length(zonkB) == 0
  zonkB = length(mm);
end

%zonk = 1:zonk;
zonk = zonkA : zonkB;
fprintf(1,'going from %4i/%2i to %4i/%2i .. found %3i time points out of %3i\n',StartY,StartYM,StopY,StopYM,length(zonk),length(yy))

disp('forgot OLR')
iL3orCLIMCAPS = +1;
iNumCycles = 4;
if iL3orCLIMCAPS > 0
  for ii = 1 : 72
    fprintf(1,'   .... do_the_fits_airsL3_ratesv7_tiles.m : cloud and OLR stuff for AIRS L3 only : lonbin %2i of 72 \n\n',ii)
    xthestats64x72_other = do_profilerate_fit_other_fraction(squeeze(save64x72_CO(:,ii,:,zonk)),squeeze(save64x72_olr(:,ii,zonk)),squeeze(save64x72_clrolr(:,ii,zonk)),days(zonk),rlat,iNumCycles);

    thestats64x72_other.olrrate(ii,:) = xthestats64x72_other.olrrate;
    thestats64x72_other.olrratestd(ii,:) = xthestats64x72_other.olrratestd;
    thestats64x72_other.clrolrrate(ii,:) = xthestats64x72_other.clrolrrate;
    thestats64x72_other.clrolrratestd(ii,:) = xthestats64x72_other.clrolrratestd;
    thestats64x72_other.olranom(ii,:,:)  = xthestats64x72_other.olranom;
    thestats64x72_other.clr_olranom(ii,:,:)  = xthestats64x72_other.clrolranom;
  end
end

save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc_oopOLRanom.mat thestats64x72_other
