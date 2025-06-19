%% based on /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA/clustbatch_make_merra2_monthly_surfaceT.m

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/

disp('also see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/clustbatch_make_L3_airsV7_climcaps_surfaceT.m')

%{
%% check things done using 
iCnt = 0;
for yy = 2002 : 2022
  mS = 01; mE = 12;
  if yy == 2002
    mS = 09;
  elseif yy == 2022
    mE = 08;
  end
  for mm = mS : mE
    iCnt = iCnt + 1;
    fname = ['L3_SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(yy) '_' num2str(mm,'%02d') '.mat'];
    if ~exist(fname)
      fprintf(1,'%3i \n',iCnt);
    end
  end
end

%% then do eg    kleenslurm; sbatch -p 2021 --array=227 sergio_matlab_chip.sbatch 17
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_XX_YY_from_X_Y

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 120;
  JOB = 227;
end

yy = 2002;
mm = 08;
for ii = 1 : JOB
  mm = mm + 1;
  if mm > 12
    yy = yy + 1;
    mm = 1;
  end
  ysave(ii) = yy;
  msave(ii) = mm;
  fprintf(1,'%3i %4i/%2i \n',ii,ysave(ii),msave(ii))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%driver_compute_AIRS_CLIMCAPS_trends_desc_or_asc.m:61:maindir = '/asl/airs/AIRS3STM/v7/';
%driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m:74:maindir = '/asl/airs/AIRS3STM/v7/';
%driver_compute_AIRSL3_trends_desc_or_asc.m:88:maindir = '/asl/airs/AIRS3STM/v7/'; %% https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STM.006/2018/
%driver_compute_AIRSL3_trends_desc_or_ascNOQuestioN.m:88:maindir = '/asl/airs/AIRS3STM/v7/'; %% https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STM.006/2018/
fnameV7   = ['/home/sergio/asl/airs/AIRS3STM/v7/'            num2str(yy) '/AIRS.' num2str(yy) '.' num2str(mm,'%02d') '.01.L3.RetStd_IR031.v7.0.3.0.G20217190152.hdf'];
fnameCLIM = ['/home/sergio/asl/airs/CLIMCAPS_SNDR_AIRS_L3/'  num2str(yy) '/SNDR.AQUA.AIRS.' num2str(yy) num2str(mm,'%02d') '01.M01.L3_CLIMCAPS_QCC.std.v02_38.G.210410074410.nc'];

%%%%%%%%%%%%%%%%%%%%%%%%%

junk    = ['/home/sergio/asl/airs/AIRS3STM/v7/'            num2str(yy) '/AIRS.' num2str(yy) '.' num2str(mm,'%02d') '.01.L3.RetStd_IR*'];
thedir  = dir(junk);
fnameV7 = ['/home/sergio/asl/airs/AIRS3STM/v7/'            num2str(yy) '/' thedir(1).name];

Airs_Lat = hdfread(fnameV7, 'location', 'Fields', 'Latitude');
Airs_Lon = hdfread(fnameV7, 'location', 'Fields', 'Longitude');

Airs_STemp_A = hdfread(fnameV7, 'ascending', 'Fields','SurfSkinTemp_A');
Airs_STemp_D = hdfread(fnameV7, 'descending', 'Fields','SurfSkinTemp_D');

aV7.lon = Airs_Lon';
aV7.lat = Airs_Lat';
aV7.lon = unique(Airs_Lon');
aV7.lat = unique(Airs_Lat');
aV7.surf_temp(:,:,1) = fliplr(Airs_STemp_A');
aV7.surf_temp(:,:,2) = fliplr(Airs_STemp_D');
  bad = find(aV7.surf_temp < 100); aV7.surf_temp(bad) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%

junk      = ['/home/sergio/asl/airs/CLIMCAPS_SNDR_AIRS_L3/'  num2str(yy) '/SNDR.AQUA.AIRS.' num2str(yy) num2str(mm,'%02d') '01.M01.L3_CLIMCAPS_QCC*'];
thedir    = dir(junk);
fnameCLIM = ['/home/sergio/asl/airs/CLIMCAPS_SNDR_AIRS_L3/'  num2str(yy) '/' thedir(1).name];

aCLIM = read_netcdf_lls(fnameCLIM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1 : 64
  fprintf(1,'latbin %2i of 64 \n',jj);
  for ii = 1 : 72
    tile = (jj-1)*72 + ii;

    rlonIJ(ii,jj)  = rlon(ii);    
    rlatIJ(ii,jj)  = rlat(jj);    

    %%%%%%%%%%%%%%%%%%%%%%%%% climcaps

    a = aCLIM;
    bad = find(a.surf_temp < 100); a.surf_temp(bad) = NaN;

    xlim1 = rlon73(ii);
    xlim2 = rlon73(ii+1);
    ylim1 = rlat65(jj);
    ylim2 = rlat65(jj+1);
    booX = find(a.lon >= xlim1 & a.lon < xlim2);
    booY = find(a.lat >= ylim1 & a.lat < ylim2);
    wahA = squeeze(a.surf_temp(booX,booY,1));
    wahD = squeeze(a.surf_temp(booX,booY,2));
    wahA = wahA(:);
    wahD = wahD(:);
    meanA_stempIJ_clim(ii,jj) = nanmean(wahA);
    meanD_stempIJ_clim(ii,jj) = nanmean(wahD);
    mean_rlonIJ_clim(ii,jj)  = nanmean(a.lon(booX));
    mean_rlatIJ_clim(ii,jj)  = nanmean(a.lat(booY));    
    
    xctr1 = rlon(ii);
    yctr1 = rlat(jj);
    booX = find(a.lon >= xctr1,1);
    booY = find(a.lat >= yctr1,1);
    wahA = squeeze(a.surf_temp(booX,booY,1));
    wahD = squeeze(a.surf_temp(booX,booY,2));
    cntrA_stempIJ_clim(ii,jj) = wahA;
    cntrD_stempIJ_clim(ii,jj) = wahD;
    cntr_rlonIJ_clim(ii,jj)  = a.lon(booX);
    cntr_rlatIJ_clim(ii,jj)  = a.lat(booY);    

    %%%%%%%%%%%%%%%%%%%%%%%%% v7

    a = aV7;
    bad = find(a.surf_temp < 100); a.surf_temp(bad) = NaN;

    xlim1 = rlon73(ii);
    xlim2 = rlon73(ii+1);
    ylim1 = rlat65(jj);
    ylim2 = rlat65(jj+1);
    booX = find(a.lon >= xlim1 & a.lon < xlim2);
    booY = find(a.lat >= ylim1 & a.lat < ylim2);
    wahA = squeeze(a.surf_temp(booX,booY,1));
    wahD = squeeze(a.surf_temp(booX,booY,2));
    wahA = wahA(:);
    wahD = wahD(:);
    meanA_stempIJ_v7(ii,jj) = nanmean(wahA);
    meanD_stempIJ_v7(ii,jj) = nanmean(wahD);
    mean_rlonIJ_v7(ii,jj)  = nanmean(a.lon(booX));
    mean_rlatIJ_v7(ii,jj)  = nanmean(a.lat(booY));    
    
    xctr1 = rlon(ii);
    yctr1 = rlat(jj);
    booX = find(a.lon >= xctr1,1);
    booY = find(a.lat >= yctr1,1);
    wahA = squeeze(a.surf_temp(booX,booY,1));
    wahD = squeeze(a.surf_temp(booX,booY,2));
    cntrA_stempIJ_v7(ii,jj) = wahA;
    cntrD_stempIJ_v7(ii,jj) = wahD;
    cntr_rlonIJ_v7(ii,jj)  = a.lon(booX);
    cntr_rlatIJ_v7(ii,jj)  = a.lat(booY);    

    %%%%%%%%%%%%%%%%%%%%%%%%%

  end
end

addpath /home/sergio/MATLABCODE/COLORMAP
figure(1); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_clim); colormap jet; caxis([200 340])
figure(2); scatter_coast(rlonIJ,rlatIJ,50,cntrA_stempIJ_clim); colormap jet; caxis([200 340])
figure(3); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_clim - cntrA_stempIJ_clim); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - mean_rlatIJ_clim); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - cntr_rlatIJ_clim); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - mean_rlonIJ_clim); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - cntr_rlonIJ_clim); colormap(usa2); caxis([-1 +1])

figure(1); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_v7); colormap jet; caxis([200 340])
figure(2); scatter_coast(rlonIJ,rlatIJ,50,cntrA_stempIJ_v7); colormap jet; caxis([200 340])
figure(3); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_v7 - cntrA_stempIJ_v7); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - mean_rlatIJ_v7); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlatIJ - cntr_rlatIJ_v7); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - mean_rlonIJ_v7); colormap(usa2); caxis([-1 +1])
figure(4); scatter_coast(rlonIJ,rlatIJ,50,rlonIJ - cntr_rlonIJ_v7); colormap(usa2); caxis([-1 +1])

figure(1); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_clim); colormap jet; caxis([200 340])
figure(2); scatter_coast(rlonIJ,rlatIJ,50,cntrA_stempIJ_v7); colormap jet; caxis([200 340])
figure(3); scatter_coast(rlonIJ,rlatIJ,50,meanA_stempIJ_clim - cntrA_stempIJ_v7); colormap(usa2); caxis([-1 +1])

figure(1); scatter_coast(rlonIJ,rlatIJ,50,meanD_stempIJ_clim); colormap jet; caxis([200 340])
figure(2); scatter_coast(rlonIJ,rlatIJ,50,cntrD_stempIJ_v7); colormap jet; caxis([200 340])
figure(3); scatter_coast(rlonIJ,rlatIJ,50,meanD_stempIJ_clim - cntrD_stempIJ_v7); colormap(usa2); caxis([-1 +1])

comment = 'see clustbatch_make_merra2_monthly_surfaceT.m';
saver = ['save L3_SKT_TIMESERIES_2002_09_to_2022_08/surfaceT_' num2str(yy) '_' num2str(mm,'%02d') '.mat rlonIJ rlatIJ comment '];
saver = [saver ' cntrA_stempIJ_clim cntrD_stempIJ_clim cntr_rlonIJ_clim cntr_rlatIJ_clim meanA_stempIJ_clim meanD_stempIJ_clim mean_rlonIJ_clim mean_rlatIJ_clim '];
saver = [saver ' cntrA_stempIJ_v7   cntrD_stempIJ_v7   cntr_rlonIJ_v7   cntr_rlatIJ_v7   meanA_stempIJ_v7   meanD_stempIJ_v7   mean_rlonIJ_v7   mean_rlatIJ_v7 '];
eval(saver)
