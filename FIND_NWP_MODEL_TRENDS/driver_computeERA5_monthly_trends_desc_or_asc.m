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

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m
for ii = 1 : iaMax
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
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
  error('not enough')
end

disp('reading in monthly ERA5 data in /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m')
for ii = 1 : iaMax
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
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

comment = 'see computeERA5_trends.m';
comment = 'see driver_computeERA5_monthly_trends_desc_or_asc.m';
if iNumYears == 18
  if iDorA > 0
    %save -v7.3 ERA5_atm_data_2002_09_to_2019_08_desc.mat comment all
    save -v7.3 ERA5_atm_data_2002_09_to_2020_08_desc.mat comment all
  else
    %save -v7.3 ERA5_atm_data_2002_09_to_2019_08_asc.mat comment all
    save -v7.3 ERA5_atm_data_2002_09_to_2020_08_asc.mat comment all
  end
elseif iNumYears == 19
  if iDorA > 0
%    save -v7.3 ERA5_atm_data_2002_09_to_2021_07_desc.mat comment all
    save -v7.3 ERA5_atm_data_2002_09_to_2021_08_desc.mat comment all
  else
%    save -v7.3 ERA5_atm_data_2002_09_to_2021_07_asc.mat comment all
    save -v7.3 ERA5_atm_data_2002_09_to_2021_08_asc.mat comment all
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(all.stemp,1)); colormap(jet); title('ERA5 mean stemp')
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(all.RHSurf,1)); colormap(jet); title('ERA5 mean RHsurf')
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(all.TwSurf,1)); colormap(jet); title('ERA5 mean TWSurf')
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(all.mmw,1)); colormap(jet); title('ERA5 mean mmw')

figure(5); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.stemp); colormap(jet); title('ERA5 mean stemp')
figure(6); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.RHSurf); colormap(jet); title('ERA5 mean RHsurf')
figure(7); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.TwSurf); colormap(jet); title('ERA5 mean TWSurf')
figure(8); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.mmw); colormap(jet); title('ERA5 mean mmw')

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
  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(10); junk = reshape(a.pnew_op.RH,100,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([00 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

computeERA5_surface_trends

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');    caxis([-0.2 +0.2]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf K/yr');  caxis([-0.2 +0.2]); colormap(usa2);
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');  caxis([-0.2 +0.2]); colormap(usa2);

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');    caxis([-0.1 +0.1]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
figure(3); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf K/yr');  caxis([-0.1 +0.1]); colormap(usa2);
figure(4); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');  caxis([-0.2 +0.2]); colormap(usa2);
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

if iNumYears == 18
  if iDorA > 0
    %save ERA5_atm_data_2002_09_to_2019_08_trends_desc.mat comment trend*
    save ERA5_atm_data_2002_09_to_2020_08_trends_desc.mat comment trend*
  else
    %save ERA5_atm_data_2002_09_to_2019_08_trends_asc.mat comment trend*
    save ERA5_atm_data_2002_09_to_2020_08_trends_asc.mat comment trend*
  end
elseif iNumYears == 19
  if iDorA > 0
    %save ERA5_atm_data_2002_09_to_2021_07_trends_desc.mat comment trend*
    save ERA5_atm_data_2002_09_to_2021_08_trends_desc.mat comment trend*
  else
    %save ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat comment trend*
    save ERA5_atm_data_2002_09_to_2021_08_trends_asc.mat comment trend*
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');    caxis([-0.1 +0.1]); colormap(usa2);
figure(2); scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);

figure(3); junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 trend ptemp K/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 trend RH percent/yr');  caxis([-0.25 +0.25]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

figure(5); junk = squeeze(nanmean(all.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 mean ptemp K');  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(6); junk = squeeze(nanmean(all.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 mean RH percent');  caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%{
figure(7); junk = squeeze(nanstd(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 std ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(8); junk = squeeze(nanstd(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 std RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(7); junk = squeeze(max(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 std ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(8); junk = squeeze(max(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 std RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7); junk = reshape(trend_nwp_ptemp,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 trend ptemp K/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(8); junk = reshape(trend_nwp_rh,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 trend RH percent/yr');  caxis([-0.25 +0.25]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(9); junk = reshape(trend_nwp_gg,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 trend SH g/g/yr');  caxis([0 +2.5]*1e-5); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(10); junk = reshape(trend_nwp_ppmv,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 trend PPMV /yr');  caxis([0 40]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(11); junk = reshape(trend_nwp_frac,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 frac /yr');  caxis([-10 +10]*1e-3); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

pause(0.1);

junk1 = squeeze(nanmean(all.nwp_gas_1,1)); junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2)); 
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(12); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 SH g/g/yr VERS2');  caxis([0 +2.5]*1e-5); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = toppmv(all.nwp_plevs,all.nwp_ptemp,all.nwp_gas_1,18,21); junk1 = squeeze(nanmean(junk1,1));
junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2));
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(13); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1);
figure(13); loglog(nanmean(junk1,2),trend_nwp_plevs_mean);  set(gca,'ydir','reverse'); xlim([1 1e4]); grid
figure(13); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 PPMV/yr VERS2');  caxis([0 40]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar


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
