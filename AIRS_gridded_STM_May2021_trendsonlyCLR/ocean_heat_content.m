function ohc = ocean_heat_content()

%% https://psl.noaa.gov/data/gridded/tables/ocean.html
%% https://www.ncei.noaa.gov/access/global-ocean-heat-content/monthly_analysis.html
%% look at Globally analyzed fields https://www.ncei.noaa.gov/access/global-ocean-heat-content/heat_monthly_global.html and get the netcds files
%% Howver this paper : https://essd.copernicus.org/articles/16/3517/2024/ : says the OHC (ocean heat contents) is in ZK = 10^21 Joules ie 1000 more
%{
Analyzed fields:
0 - 700 meters layer (2005 to present)
File naming conventions:
[p]_[ud]-[ld]_[by][ey][bm]-[em].dat where:
[p] - parameter (HC = heat content)
[ud] - upper depth of layer
[ld] - lower depth of layer
[by] - last two digits of beginning year
[ey] - last two digits of ending year
[bm] - two digits of beginning month
[em] - two digits of ending month
Note: year 2000 is presented as A0, 2001 is A1, etc.
year 2010 is presented as B0, 2011 is B1, etc.
year 2020 is presented as C0, 2021 is C1, etc.

Monthly heat content (10^18 joules) at 1 m, and everyhing else at 10^22 J = 10 ZJ (1ZJ = 10^21 J)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
from NCDUMP

     float h18_hc(time, depth, lat, lon) ;
     	   h18_hc:long_name = "Ocean heat content anomaly calculated from objectively analyzed temperature anomaly fields over given depth limits" ;
           h18_hc:coordinates = "time lat lon depth" ;
           h18_hc:cell_methods = "area: mean depth: mean time: mean" ;
	   h18_hc:grid_mapping = "crs" ;
           h18_hc:units = "10^18_joules" ;

	   float month_h22_WO(time) ;
	   	 month_h22_WO:long_name = "global_heat_content_integral_0-700 " ;
                 month_h22_WO:units = "10^22_joules" ;
    	   float month_h22_se_WO(time) ;
	         month_h22_se_WO:long_name = "global_standarderror_heat_content_integral_0-700" ;
	  	 month_h22_se_WO:units = "10^22_joules" ;
           float month_h22_NH(time) ;
	         month_h22_NH:long_name = "northern_hemisphere_heat_content_integral_0-700 " ;
		 month_h22_NH:units = "10^22_joules" ;
          float month_h22_se_NH(time) ;
       	         month_h22_se_NH:long_name = "northern_hemisphere_standarderror_heat_content_integral_0-700" ;
		 month_h22_se_NH:units = "10^22_joules" ;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('YY')

  %% run this as stand alone to check anonaly code
  %% run this as stand alone to check anonaly code
  %% run this as stand alone to check anonaly code

  do_XX_YY_from_X_Y
  addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
  addpath /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE_Git/COLORMAP
  addpath /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE_Git/PLOTTER
  addpath /home/sergio/MATLABCODE/TIME
  addpath /asl/matlib/aslutil
end

ohc.ocean_0700 = read_netcdf_lls('heat_content_anomaly_0-700_monthly.nc');
ohc.ocean_2000 = read_netcdf_lls('heat_content_anomaly_0-2000_monthly.nc');

ohc.oceantime = 2005 + ((1:231)-1)/12;
pcolor(squeeze(ohc.ocean_0700.h18_hc(:,:,1,1))'); colorbar; shading interp; colorbar; caxis([-1 +1]*80); colormap(usa2);
  title('1 m heat content 10^18 Joules')
plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO,'r.-',ohc.oceantime,ohc.ocean_0700.month_h22_NH,'b',ohc.oceantime,ohc.ocean_0700.month_h22_SH,'b--',ohc.oceantime,ohc.ocean_0700.month_h22_NH+ohc.ocean_0700.month_h22_SH,'m')
  plotaxis2; hl = legend('All Ocean','NH','SH','NH+SH','location','best','fontsize',10); title('0-0700m'); ylabel('10^{22} Joules')
plot(ohc.oceantime,ohc.ocean_2000.month_h22_WO,'r.-',ohc.oceantime,ohc.ocean_2000.month_h22_NH,'b',ohc.oceantime,ohc.ocean_2000.month_h22_SH,'b--',ohc.oceantime,ohc.ocean_2000.month_h22_NH+ohc.ocean_2000.month_h22_SH,'m')
  plotaxis2; hl = legend('All Ocean','NH','SH','NH+SH','location','best','fontsize',10); title('0-2000m'); ylabel('10^{22} Joules')

plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO,'r',ohc.oceantime,ohc.ocean_2000.month_h22_WO,'b'); title('All Ocean');
  plotaxis2; hl = legend('0-0700','0-2000','location','best','fontsize',10); title('Heat Content'); ylabel('10^{22} Joules')
plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO*10,'r',ohc.oceantime,ohc.ocean_2000.month_h22_WO*10,'b'); title('All Ocean');
  plotaxis2; hl = legend('0-0700','0-2000','location','best','fontsize',10); title('Heat Content'); ylabel('ZettaJoules (1 ZJ = 10^{21} J)')
%% see Fig 12 of https://essd.copernicus.org/articles/16/3517/2024/
plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO*10-200,'r',ohc.oceantime,ohc.ocean_2000.month_h22_WO*10-200,'b'); title('All Ocean');
  plotaxis2; hl = legend('0-0700','0-2000','location','best','fontsize',10); title('Heat Content'); ylabel('ZettaJoules (1 ZJ = 10^{21} J)')

addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
boo = 1 : length(ohc.oceantime);
[B, stats] = Math_tsfit_lin_robust((ohc.oceantime-ohc.oceantime(1))*365.25,ohc.ocean_0700.month_h22_WO,4);
ohc.ocean_0700.month_h22_WO_anom = compute_anomaly(boo,(ohc.oceantime-ohc.oceantime(1))*365.25,B,[],ohc.ocean_0700.month_h22_WO,-1);
  plot(1:231,ohc.ocean_0700.month_h22_WO_anom+B(1),1:231,ohc.ocean_0700.month_h22_WO)
[B, stats] = Math_tsfit_lin_robust((ohc.oceantime-ohc.oceantime(1))*365.25,ohc.ocean_2000.month_h22_WO,4);
ohc.ocean_2000.month_h22_WO_anom = compute_anomaly(boo,(ohc.oceantime-ohc.oceantime(1))*365.25,B,[],ohc.ocean_2000.month_h22_WO,-1);
  plot(1:231,ohc.ocean_2000.month_h22_WO_anom+B(1),1:231,ohc.ocean_2000.month_h22_WO)

%% see /home/sergio/MATLABCODE/MODIS_CLOUD/driver_read_modis_aerosol_monthly_L3.m
do_XX_YY_from_X_Y

for ii = 1 : 72
  for jj = 1 : 64
    booX = find(ohc.ocean_0700.lon_bnds(1,:) >= rlon73(ii) & ohc.ocean_0700.lon_bnds(2,:) < rlon73(ii+1));
    booY = find(ohc.ocean_0700.lat_bnds(2,:) >= rlat65(jj) & ohc.ocean_0700.lat_bnds(1,:) < rlat65(jj+1));
    for iCnt = 1 : 231
      junk = squeeze(ohc.ocean_0700.h18_hc(booX,booY,1,iCnt));
      ohc.ocean_0700.h18_hc_72x64(ii,jj,iCnt) = nanmean(junk(:));
      junk = squeeze(ohc.ocean_2000.h18_hc(booX,booY,1,iCnt));
      ohc.ocean_2000.h18_hc_72x64(ii,jj,iCnt) = nanmean(junk(:));
    end
  end
end

iPlotOceanAnom = -1;
if iPlotOceanAnom > 0
  for iCnt = 1 : 231
    pcolor(reshape(XX,72,64),reshape(YY,72,64),squeeze(ohc.ocean_0700.h18_hc_72x64(:,:,iCnt))); 
    simplemap(reshape(YY,72,64),reshape(XX,72,64),squeeze(ohc.ocean_0700.h18_hc_72x64(:,:,iCnt)),5);
    title(num2str(iCnt,'%03d'));    title([num2str(ohc.oceantime(iCnt)) ' anomaly at 1 m in 10^{18} J']); 
    shading interp;  colorbar; caxis([-1 +1]*50); colormap(usa2); pause(0.1);
  end
end

ohc.iNumOceantimeSteps = length(ohc.oceantime);
coslatOcean = cos(YY'*pi/180)*ones(1,ohc.iNumOceantimeSteps);
coslatOcean = cos(YY'*pi/180);
iCnt = 1;
for ii = 1 : ohc.iNumOceantimeSteps
  moo = squeeze(ohc.ocean_0700.h18_hc_72x64(:,:,ii));
  ohc.ocean_0700.h18_hc_72x64_cos(ii) = nansum(moo(:).*coslatOcean)./nansum(coslatOcean);
  moo = squeeze(ohc.ocean_2000.h18_hc_72x64(:,:,ii));
  ohc.ocean_2000.h18_hc_72x64_cos(ii) = nansum(moo(:).*coslatOcean)./nansum(coslatOcean);
end

plot(ohc.oceantime,ohc.ocean_0700.h18_hc_72x64_cos*5,'b',ohc.oceantime,ohc.ocean_0700.month_h22_WO,'r',ohc.oceantime,ohc.ocean_2000.h18_hc_72x64_cos*20,'c',ohc.oceantime,ohc.ocean_2000.month_h22_WO,'k'); plotaxis2; 
  hl = legend('1m energy x 5','integrated 0-0700 m','1m energy x 20','integrated 0-2000 m','location','best','fontsize',10);
plot(ohc.oceantime,ohc.ocean_0700.h18_hc_72x64_cos*5,'b',ohc.oceantime,ohc.ocean_0700.month_h22_WO,'r',ohc.oceantime,ohc.ocean_2000.h18_hc_72x64_cos*20/4,'c',ohc.oceantime,ohc.ocean_2000.month_h22_WO,'k'); plotaxis2; 
  hl = legend('1m energy x 5','integrated 0-0700 m','1m energy x 20/4','integrated 0-2000 m','location','best','fontsize',10);

P = polyfit(ohc.oceantime,ohc.ocean_0700.month_h22_WO,1); boohoo = polyval(P,ohc.oceantime);
plot(ohc.oceantime,ohc.ocean_0700.h18_hc_72x64_cos*5,ohc.oceantime,ohc.ocean_0700.month_h22_WO - boohoo'); plotaxis2;

P = polyfit(ohc.oceantime,ohc.ocean_2000.month_h22_WO,1); boohoo = polyval(P,ohc.oceantime);
plot(ohc.oceantime,ohc.ocean_2000.h18_hc_72x64_cos*5,ohc.oceantime,ohc.ocean_2000.month_h22_WO - boohoo'); plotaxis2;

[B07, stats07, junkanom07] = compute_anomaly_wrapper(1:length(ohc.oceantime),ohc.oceantime-ohc.oceantime(1),ohc.ocean_0700.month_h22_WO,4,[],-1,-1);
[B20, stats20, junkanom20] = compute_anomaly_wrapper(1:length(ohc.oceantime),ohc.oceantime-ohc.oceantime(1),ohc.ocean_2000.month_h22_WO,4,[],-1,-1);
plot(ohc.oceantime,junkanom07,'b',ohc.oceantime,junkanom20,'r','linewidth',2); plotaxis2; hl = legend('0700m','2000m','location','best'); title('Ocean Heat Content Anomaly'); ylabel('10^{22} J');
