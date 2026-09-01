load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/h2645structure.mat

iIndia = 2788;
iIndia = 2786;  %% Arabian Sea
iIndia = 2787;  %% Arabian Sea off India, and Maharastra
iIndia = 2788;  %% central India, next to MH

iIndia = input('Enter tile (default 2788 = central India) : ');
if length(iIndia) == 0
  iIndia = 2788;
end

figure(1); clf; plot_72x64_tiles(iIndia);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf; 

era5spectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat','trend');
plot(h.vchan,era5spectra.trend(:,iIndia),'b'); hold on
era5spectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat','trend');
plot(h.vchan,era5spectra.trend(:,iIndia),'c'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,obsspectra.b_asc(:,iIndia),'r',h.vchan,obsspectra.b_desc(:,iIndia),'m'); hold off
plotaxis2; 

xlim([640 1640])
hl = legend('ERA5 D','ERA5 N','OBS D','OBS N','location','best');
title('ERA5 + AIRS spectral rate for India')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_put_together_QuantileChoose_anomalies.m

iCnt = 0;
for yy = 2002 : 2022
  mmS = 1; mmE = 12;
  if yy == 2002
    mmS = 09;
  elseif yy == 2022
    mmS = 08;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;
    daysSince2002(iCnt) = change2days(yy,mm,15,12);
  end
end

moo = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin38/LonBin51/iQAX_3_summarystats_LatBin38_LonBin51_timesetps_001_457_V1.mat');

%plot(2002 + daysSince2002/365,btanomA(iIndia,:),'b',2002 + daysSince2002/365,btanomD(iIndia,:),'r','linewidth',2); %% LatBin = 39, LonBin = 51 >> (39-1)*72 + 51 = iIndia
%plotaxis2;
%hl = legend('asc','desc','location','best'); title(['BTXXXX anom, tile iIndia over India/Arabian Sea Q' num2str(iQuant)]); xlim([2002.75 2022.75])

figure(3); clf
mootime = moo.year_desc + (moo.month_desc-1)/12 + (moo.day_desc-1)/12/30;

plot(mootime,moo.meanBT1231_asc+10,'b',mootime,moo.meanBT1231_desc-10,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile iIndia')
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.meanBT1231_asc,4);  fprintf(1,'BT1231 trend for asc node = %8.6f K/yr \n',B(2))
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.meanBT1231_desc,4); fprintf(1,'BT1231 trend for desc node = %8.6f K/yr \n',B(2))

figure(4); clf
moo.rad1231_asc  = squeeze(moo.rad_quantile_asc(:,1520,5));
moo.rad1231_desc = squeeze(moo.rad_quantile_desc(:,1520,5));
plot(mootime,moo.rad1231_asc+0,'b',mootime,moo.rad1231_desc-0,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean RAD1231 for Tile iIndia, Q05')
plot(mootime,moo.rad1231_asc+5,'b',mootime,moo.rad1231_desc-5,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean RAD1231 for Tile iIndia, Q05')
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.rad1231_asc,4);  fprintf(1,'rad1231 trend Q05 for asc node = %8.6f K/yr \n',B(2))
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.rad1231_desc,4); fprintf(1,'rad1231 trend Q05 for desc node = %8.6f K/yr \n',B(2))

figure(5); clf
moo.bt1231_asc  = rad2bt(1231,squeeze(moo.rad_quantile_asc(:,1520,5)));
moo.bt1231_desc = rad2bt(1231,squeeze(moo.rad_quantile_desc(:,1520,5)));
plot(mootime,moo.bt1231_asc+0,'b',mootime,moo.bt1231_desc-0,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile iIndia, Q05')
%plot(mootime,moo.bt1231_asc+5,'b',mootime,moo.bt1231_desc-5,'r','linewidth',2); hl = legend('asc','desc','fontsize',10); xlim([2002.75 2022.75]); title('Mean BT1231 for Tile iIndia, Q05')
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.bt1231_asc,4);  fprintf(1,'bt1231 trend Q05 for asc node = %8.6f K/yr \n',B(2))
B = Math_tsfit_lin_robust((mootime-2002)*365,moo.bt1231_desc,4); fprintf(1,'bt1231 trend Q05 for desc node = %8.6f K/yr \n',B(2))

figure(6); clf
sktjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/surface_jac_new.mat');
tjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac_new.mat');
g1jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');
g101jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g101_jac_new.mat');
g102jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');
plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2))
plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2)*0.1); xlim([640 1640]); plotaxis2;

figure(7); clf
sktjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/surface_jac.mat');
tjac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/temp_jac.mat');
g1jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
g101jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g101_jac.mat');
g102jac = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2))
plot(sktjac.fout,sktjac.jsurface(:,1),sktjac.fout,sum(g1jac.jout + g101jac.jout + g102jac.jout,2)*0.1); xlim([640 1640]); plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8); clf; 

afrac = 1;
nfrac = 0;
obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'b--',h.vchan,nfrac*obsspectra.b_desc(:,iIndia),'m--'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'r--',h.vchan,nfrac*obsspectra.b_desc(:,iIndia),'m--'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'k',h.vchan,nfrac*obsspectra.b_desc(:,iIndia),'g'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'b',h.vchan,nfrac*obsspectra.b_desc(:,iIndia),'c'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'r',h.vchan,nfrac*obsspectra.b_desc(:,iIndia),'m'); hold on

hold off
plotaxis2; 

xlim([640 1640])
hl = legend('D 01','N 01','D 02','N 02','D 03','N 03','D 04','N 04','D 05','N 05','location','best');
title('AIRS spectral rate for India')

%%%%%%%%%%%%%%%%%%%%%%%%%

afrac = 1;
nfrac = 0;
obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'c'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'g'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'k'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'b'); hold on

obsspectra = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat','b_asc','b_desc');
obsspectra.b_asc  = (reshape(obsspectra.b_asc,4608,2645))';
obsspectra.b_desc = (reshape(obsspectra.b_desc,4608,2645))';
plot(h.vchan,afrac*obsspectra.b_asc(:,iIndia),'r'); hold on

hold off
plotaxis2; 

xlim([640 1640])
hl = legend('01','02','03','04','05','location','best');
title('AIRS spectral rate for India')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
