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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('takes about 2-3 hours')

system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 20 for the iNumYears 
fprintf(1,'JOB = %2i \n',JOB)

if ~exist('iDorA')
  iDorA = -1; %% asc
  iDorA = +1; %% desc
end

if iDorA == -1
  error('errrrrrrr actually monthly MERRA2 was one time, so cannot really do ASC vs DESC')
end

clear iaFound
%iaMax = 18*12; %% 18 year
%iaMax = 19*12; %% 19 year

iNumYears = 070; %% 2012/05-2019/04 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iNumYears = 12; %% 2002/09 to 2014/08 
iNumYears = 18; %% 2002/09 to 2020/08 
iNumYears = 19; %% 2002/09 to 2021/08 
iNumYears = 20; %% 2002/09 to 2022/08 
%iaMax = iNumYears*12;

% iNumYears = input('Enter iNumYears : ');
if length(JOB) == 0
  JOB = 20;
end
iNumYears = JOB;

fprintf(1,'iNumYears = %2i \n',iNumYears);

if iNumYears <= 69
  iaMax = iNumYears*12;
elseif iNumYears == 70
  iaMax = 7*12;
end

iAllorSeasonal = +1;  %% all
iAllorSeasonal = -1;  %% DJF
iAllorSeasonal = -2;  %% MAM

iAllorSeasonal = -4;  %% SON
iAllorSeasonal = -3;  %% JJA

[iDorA iaMax]

find_computeMERRA2_monthly_trends_foutname
if exist(fout_trendjunk)
  fout_trendjunk
  error('fout_trendjunk exists')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m
for ii = 1 : iaMax
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/DESC/merra2_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/ASC/merra2_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  end
  if exist(fin)
    iaFound(ii) = 1;
  else
    iaFound(ii) = 0;
  end
end
fprintf(1,'found %3i of expected %3i files \n', [sum(iaFound) length(iaFound)])
plot(1:iaMax,iaFound,'+-')
if sum(iaFound) < length(iaFound)
  find(iaFound == 0)
  fprintf(1,'last file looked for = %s \n',fin)
  error('not enough')
end

disp('mat files made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA2/clust_loop_make_monthly_tile_center.m')
disp('mat files made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA2/clust_loop_make_monthly_tile_center.m')
disp('mat files made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA2/clust_loop_make_monthly_tile_center.m')

disp('reading in monthly MERRA2 data in /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/ made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MERRA2/clust_loop_make_monthly_tile_center_asc_or_desc.m')
for ii = 1 : iaMax
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/DESC/merra2_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/ASC/merra2_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  end
  if exist(fin)
    if mod(ii,100) == 0
      fprintf(1,'+ \n')
    elseif mod(ii,10) == 0
      fprintf(1,'x')
    else
      fprintf(1,'.')
    end
    iaFound(ii) = 1;

    a = load(fin);
    %% already done this in /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_CMIP6/clust_compute_cmip6_profile_rtpfiles.m
    %a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'mixing ratio','relative humidity');
    if ~isfield(a.pnew_ip,'rh')
      a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'specific humidity','relative humidity');
    end

    all.yy(ii) = a.thedateS(1);
    all.mm(ii) = a.thedateS(2);
    all.dd(ii) = a.thedateS(3);

    all.nwp_ptemp(ii,:,:) = a.pnew_ip.ptemp;
    all.nwp_gas_1(ii,:,:) = a.pnew_ip.gas_1;
    all.nwp_gas_3(ii,:,:) = a.pnew_ip.gas_3;
    all.nwp_rh(ii,:,:)    = a.pnew_ip.rh;
    all.nwp_plevs(ii,:,:) = a.pnew_ip.plevs;

    all.gas_1(ii,:,:) = a.pnew_op.gas_1;
    all.gas_3(ii,:,:) = a.pnew_op.gas_3;
    all.ptemp(ii,:,:) = a.pnew_op.ptemp;
    all.stemp(ii,:)   = a.pnew_op.stemp;
    all.mmw(ii,:)     = a.pnew_op.mmw;
    all.nlays(ii,:)   = a.pnew_op.nlevs-1;
    all.RH(ii,:,:)    = a.pnew_op.RH;
    all.TwSurf(ii,:)  = a.pnew_op.TwSurf;
    all.RHSurf(ii,:)  = a.pnew_op.RHSurf;
  else
    iaFound(ii) = 0;
  end
end
fprintf(1,'\n');
all.rlon = a.pnew_op.rlon;
all.rlat = a.pnew_op.rlat;

monitor_memory_whos

comment = 'see computeMERRA2_trends.m';
comment = 'see driver_computeMERRA2_monthly_trends_desc_or_asc.m';

find_computeMERRA2_monthly_foutname

if iNumYears <= 100
  iSave = -1;
else
  iSave = input('Save huge file (-1/+1) : ');
end

iSave = -1;  %%% unless now you do a whole new set of files Steve brings down from ERA5

if iSave > 0
  %%% foutjunk = ['MERRA2_atm_data_2002_09_to_*.mat'];
  fprintf(1,'saving huge file : can type in a separate window         watch "ls -lt %s " \n',foutjunk)
  saver = ['save ' foutjunk ' comment all'];
  eval(saver);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(all.stemp,1)); colormap(jet); title('MERRA2 mean stemp')
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(all.RHSurf,1)); colormap(jet); title('MERRA2 mean RHsurf DO NOT BELIEVE')
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(all.TwSurf,1)); colormap(jet); title('MERRA2 mean TWSurf DO NOT BELIEVE')
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(all.mmw,1)); colormap(jet); title('MERRA2 mean mmw')

figure(5); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.stemp); colormap(jet); title('MERRA2 mean stemp')
figure(6); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.RHSurf); colormap(jet); title('MERRA2 mean RHsurf DO NOT BELIEVE')
figure(7); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.TwSurf); colormap(jet); title('MERRA2 mean TWSurf DN NOT BELIEVE')
figure(8); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.mmw); colormap(jet); title('MERRA2 mean mmw')

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = flipud(pN./pD);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

figure(9); junk = reshape(a.pnew_op.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar; 
  title('Mean T')
figure(10); junk = reshape(a.pnew_op.RH,100,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([00 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar; 
  title('Mean RH')
figure(11); junk = reshape(log10(a.pnew_op.gas_1),101,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([13 23]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar; 
  title('Mean log10(SH)');

pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

computeERA5_surface_trends

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('MERRA2 trend  stemp K/yr');    caxis([-0.2 +0.2]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('MERRA2 trend  RHsurf pc/yr UGH'); caxis([-0.4 +0.4]); colormap(usa2);
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('MERRA2 trend  TWSurf K/yr UGH');  caxis([-0.2 +0.2]); colormap(usa2);
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('MERRA2 trend  colwater mm/yr');  caxis([-0.2 +0.2]); colormap(usa2);

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('MERRA2 trend  stemp K/yr');    caxis([-0.1 +0.1]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('MERRA2 trend  RHsurf pc/yr UGH'); caxis([-0.4 +0.4]); colormap(usa2);
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('MERRA2 trend  TWSurf K/yr UGH');  caxis([-0.1 +0.1]); colormap(usa2);
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('MERRA2 trend  colwater mm/yr');  caxis([-0.2 +0.2]); colormap(usa2);
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeERA5_atmos_trends

if isfield(all,'nwp_plevs')
  trend_nwp_plevs_mean = mean(squeeze(nanmean(all.nwp_plevs,1)),2);
end

trend_plays = flipud(pN./pD);

trend_rlat = all.rlat;
trend_rlon = all.rlon;
trend_rlat64 = rlat; trend_rlon72 = rlon;
%trend_plevs37 = permute(all.nwp_plevs,[2 1 3]); trend_plevs37 = reshape(trend_plevs37,37,227*4608); trend_plevs37 = mean(trend_plevs37,2);

%find_computeMERRA2_monthly_trends_foutname
if ~exist(fout_trendjunk)
  fprintf(1,'saving trend file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk)
  saver = ['save ' fout_trendjunk ' comment trend*'];
  eval(saver);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('DO NOT BELIEVE RHSURF stuff since I forgot to use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('DO NOT BELIEVE RHSURF stuff since I forgot to use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('DO NOT BELIEVE RHSURF stuff since I forgot to use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('MERRA2 trend  stemp K/yr');    caxis([-0.1 +0.1]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('MERRA2 trend  RHsurf DO NOT BELIEVE pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
addpath /asl/matlib/maps/
aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('MERRA2 dST/dt');      caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('MERRA2 UGH dRHSurf/dt'); caxis([-1 +1]*0.25); colormap(llsmap5)

figure(3); junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 trend ptemp K/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 trend RH percent/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(5); junk = reshape(trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 100 layer trend WVfrac /yr');  caxis([-1 +1]*0.01); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(6); junk = squeeze(nanmean(all.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 mean ptemp K');  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(7); junk = squeeze(nanmean(all.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 mean RH percent');  caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%{
figure(7); junk = squeeze(nanstd(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 std ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(8); junk = squeeze(nanstd(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 std RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(7); junk = squeeze(max(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 std ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(8); junk = squeeze(max(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MERRA2 std RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
%}

figure(3); ylim([1 1000])
figure(5); set(gca,'yscale','linear'); caxis([-1 +1]*0.015);
figure(4); set(gca,'yscale','linear'); caxis([-1 +1]*0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%

iNlev = 42;
trend_nwp_rh = real(trend_nwp_rh);
figure(8); junk = reshape(trend_nwp_ptemp,iNlev,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('MERRA2 trend ptemp K/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(9); junk = reshape(trend_nwp_rh,iNlev,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('MERRA2 trend RH percent/yr');  caxis([-0.25 +0.25]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(10); junk = reshape(trend_nwp_gg,iNlev,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('MERRA2 trend SH g/g/yr');  caxis([-2.5 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(11); junk = reshape(trend_nwp_ppmv,iNlev,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('MERRA2 trend PPMV /yr');  caxis([-40 40]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(12); junk = reshape(trend_nwp_frac,iNlev,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('MERRA2 frac /yr');  caxis([-10 +10]*1e-3); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = squeeze(nanmean(all.nwp_gas_1,1)); junk1 = reshape(junk1,iNlev,72,64); junk1 = squeeze(nanmean(junk1,2)); 
junk2 = reshape(trend_nwp_frac,iNlev,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(13); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('MERRA2 SH g/g/yr VERS2');  caxis([-2.5 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = toppmv(all.nwp_plevs,all.nwp_ptemp,all.nwp_gas_1,18,21); junk1 = squeeze(nanmean(junk1,1));
junk1 = reshape(junk1,iNlev,72,64); junk1 = squeeze(nanmean(junk1,2));
junk2 = reshape(trend_nwp_frac,iNlev,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1);
figure(14); loglog(nanmean(junk1,2),trend_nwp_plevs_mean);  set(gca,'ydir','reverse'); xlim([1 1e4]); grid
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('MERRA2 PPMV/yr VERS2');  caxis([0 40]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% this is MLS comparisons, Frank Werner JPL
figure(11); ylim([1 300]); colormap usa2
figure(07); ylim([1 300]);
figure(09); ylim([1 300]); caxis([-1 1]*1e-7); colormap(usa2)
figure(12); ylim([1 300]); caxis([-1 1]*1e-7); colormap(usa2)
figure(10); ylim([1 300]); caxis([-1 1]*1e-1); colormap(usa2)
figure(13); ylim([1 300]); caxis([-1 1]*1e-1); colormap(usa2)
%}
