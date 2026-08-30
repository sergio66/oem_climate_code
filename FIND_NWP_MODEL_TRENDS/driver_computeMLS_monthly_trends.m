%% monthly, 18 years x 12 months/year = 216
%% monthly, 19 years x 12 months/year = 228

addpath0

%if ~exist('iDorA')
%  iDorA = +1; %% desc
%  iDorA = -1; %% asc
%end

clear iaFound

iNumYears = 16;  %% 2004-2020
iNumYears = 17;  %% 2004-2021
iNumYears = 21;  %% 2004-2025
iNumYears = 18;  %% 2004-2022

iaMax = iNumYears*12;

%% see /home/sergio/git/rtpmake/CLUST_RTPMAKE/CLUSTMAKE_MLS/clust_compute_mls_profile_rtpfiles.m
dirMLS = '/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon//TimeSeries/MLS/Tile_Center/';

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MLS/clust_compute_mls_profile_rtpfiles.m
for ii = 1 : iaMax
  fin = [dirMLS 'mls_tile_center_monthly_timestep_' num2str(ii,'%03d') '.mat'];

  if exist(fin)
    iaFound(ii) = 1;
  else
    iaFound(ii) = 0;
  end
end
[sum(iaFound) length(iaFound)]
plot(1:iaMax,iaFound,'+-')

for ii = 1 : iaMax
  fin = [dirMLS 'mls_tile_center_monthly_timestep_' num2str(ii,'%03d') '.mat'];
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

    %% already done this in /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_MLS/clust_compute_mls_profile_rtpfiles.m
    %%a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'mixing ratio','relative humidity');
    if ~isfield(a.pnew_ip,'rh')
      a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'specific humidity','relative humidity');
    end

    pall.yy(ii) = a.yyuseII;
    pall.mm(ii) = a.mmuseII;
    pall.dd(ii) = 15;

    pall.nwp_ptemp(ii,:,:) = a.pnew_ip.ptemp;
    pall.nwp_gas_1(ii,:,:) = a.pnew_ip.gas_1;
    %pall.nwp_gas_3(ii,:,:) = a.pnew_ip.gas_3;
    pall.nwp_rh(ii,:,:)    = a.pnew_ip.rh;
    pall.nwp_plevs(ii,:,:) = a.pnew_ip.plevs;

    pall.gas_1(ii,:,:) = a.pnew_op.gas_1;
    pall.gas_3(ii,:,:) = a.pnew_op.gas_3;
    pall.ptemp(ii,:,:) = a.pnew_op.ptemp;
    pall.stemp(ii,:)   = a.pnew_op.stemp;
    pall.mmw(ii,:)     = a.pnew_op.mmw;
    pall.nlays(ii,:)   = a.pnew_op.nlevs-1;
    pall.RH(ii,:,:)    = a.pnew_op.RH;
    pall.TwSurf(ii,:)  = a.pnew_op.TwSurf;
    pall.RHSurf(ii,:)  = a.pnew_op.RHSurf;
  else
    iaFound(ii) = 0;

    pall.yy(ii) = NaN;
    pall.mm(ii) = NaN;
    pall.dd(ii) = NaN;;

    pall.nwp_ptemp(ii,:,:)  = NaN;
    pall.nwp_gas_1(ii,:,:)  = NaN;
    %pall.nwp_gas_3(ii,:,:) = NaN;
    pall.nwp_rh(ii,:,:)     = NaN;
    pall.nwp_plevs(ii,:,:)  = NaN;

    pall.gas_1(ii,:,:) = NaN;
    pall.gas_3(ii,:,:) = NaN;
    pall.ptemp(ii,:,:) = NaN;
    pall.stemp(ii,:)   = NaN;
    pall.mmw(ii,:)     = NaN;
    pall.nlays(ii,:)   = NaN;
    pall.RH(ii,:,:)    = NaN;
    pall.TwSurf(ii,:)  = NaN;
    pall.RHSurf(ii,:)  = NaN;

  end
end

fprintf(1,'\n');
pall.rlon = a.pnew_op.rlon;
pall.rlat = a.pnew_op.rlat;

monitor_memory_whos

comment = 'see driver_computeMLS_monthly_trends.m';
if iNumYears == 16
  save -v7.3 MLS_atm_data_2004_09_to_2020_08.mat comment pall
elseif iNumYears == 17
  save -v7.3 MLS_atm_data_2004_09_to_2021_08.mat comment pall
elseif iNumYears == 18
  save -v7.3 MLS_atm_data_2004_09_to_2022_08.mat comment pall
elseif iNumYears == 21
  save -v7.3 MLS_atm_data_2004_09_to_2025_08.mat comment pall
end

figure(1); scatter_coast(pall.rlon,pall.rlat,40,nanmean(pall.stemp,1)); colormap(jet); title('MLS mean stemp')
figure(2); scatter_coast(pall.rlon,pall.rlat,40,nanmean(pall.RHSurf,1)); colormap(jet); title('MLS mean RHsurf')
figure(3); scatter_coast(pall.rlon,pall.rlat,40,nanmean(pall.TwSurf,1)); colormap(jet); title('MLS mean TWSurf')
figure(4); scatter_coast(pall.rlon,pall.rlat,40,nanmean(pall.mmw,1)); colormap(jet); title('MLS mean mmw')

figure(5); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.stemp); colormap(jet); title('MLS mean stemp')
figure(6); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.RHSurf); colormap(jet); title('MLS mean RHsurf')
figure(7); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.TwSurf); colormap(jet); title('MLS mean TWSurf')
figure(8); scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.mmw); colormap(jet); title('MLS mean mmw')

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = flipud(pN./pD);

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

figure(9); junk = reshape(a.pnew_op.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
  title('Mean T')
figure(10); junk = reshape(a.pnew_op.RH,100,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([00 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
  title('Mean RH')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(pall.yy,pall.mm,pall.dd,2002);

computeERA5_surface_trends

figure(1); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_stemp,1)); title('MLS trend  stemp K/yr');    caxis([-0.2 +0.2]); colormap(usa2);
figure(2); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_RHSurf,1)); title('MLS trend  RHsurf pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
figure(3); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_TwSurf,1)); title('MLS trend  TWSurf K/yr');  caxis([-0.2 +0.2]); colormap(usa2);
figure(4); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_mmw,1)); title('MLS trend  colwater mm/yr');  caxis([-0.2 +0.2]); colormap(usa2);

figure(1); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_stemp,1)); title('MLS trend  stemp K/yr');    caxis([-0.1 +0.1]/10); colormap(usa2);
figure(2); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_RHSurf,1)); title('MLS trend  RHsurf pc/yr'); caxis([-0.4 +0.4]/10); colormap(usa2);
figure(3); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_TwSurf,1)); title('MLS trend  TWSurf K/yr');  caxis([-0.1 +0.1]/10); colormap(usa2);
figure(4); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_mmw,1)); title('MLS trend  colwater mm/yr');  caxis([-0.2 +0.2]/10); colormap(usa2);
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeERA5_atmos_trends

if isfield(pall,'nwp_plevs')
  trend_nwp_plevs_mean = squeeze(nanmean(pall.nwp_plevs,1));
end

trend_plays = flipud(pN./pD);

trend_rlat = pall.rlat;
trend_rlon = pall.rlon;
trend_rlat64 = rlat; trend_rlon72 = rlon;

if iNumYears == 16
  save MLS_atm_data_2004_09_to_2020_08_trends.mat comment trend*
elseif iNumYears == 17
  save MLS_atm_data_2004_09_to_2021_08_trends.mat comment trend*
elseif iNumYears == 18
  save MLS_atm_data_2004_09_to_2022_08_trends.mat comment trend*
elseif iNumYears == 21
  save MLS_atm_data_2004_09_to_2025_08_trends.mat comment trend*
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('llsmap5.mat');
figure(1); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_stemp,1)); title('MLS trend  stemp K/yr');    caxis([-0.1 +0.1]); colormap(usa2);
figure(2); scatter_coast(pall.rlon,pall.rlat,40,nanmean(trend_RHSurf,1)); title('MLS trend  RHsurf pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);

figure(3); junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MLS trend ptemp K/yr');  caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MLS trend RH percent/yr');  caxis([-0.25 +0.25]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MLS trend WVfrac/yr');  caxis([-0.10 +0.10]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

figure(5); junk = squeeze(nanmean(pall.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MLS mean ptemp K');  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(6); junk = squeeze(nanmean(pall.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('MLS mean RH percent');  caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(7);
wah = pall.gas_1; whos wah
plot(2004 + dayOFtime/365,wah(:,[10 50 90],2000))
