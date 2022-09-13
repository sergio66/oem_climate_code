xumbc.maskLF = maskLF';

xumbc.trend_ST = results(:,6);
xumbc.trend_T  = deltaTlat(:,1:97);
xumbc.trend_RH  = deltaRHlat(:,1:97);
xumbc.trend_WV  = fracWVlat(:,1:97);

xumbc.pavgLAY = pavgLAY(1:97,3000);
xumbc.rlat65   = rlat65;
xumbc.rlat     = rlat;
xumbc.rlon73   = rlon73;

xcomment = 'see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/find_T_RH_trends.m and quickprint_find_T_RH_trends.m';
saver = ['save ' diroutQuick '/umbc_trends_iQuantile' num2str(iQuantile) '.mat xumbc xcomment iQuantile'];
eval(saver);

figure(06); clf;
figure(28); clf;
figure(29); clf;
figure(30); clf;

aslmap(6,xumbc.rlat65,xumbc.rlon73,smoothn((reshape(xumbc.maskLF.*xumbc.trend_ST,72,64)') ,1), [-90 +90],[-180 +180]); title('UMBC dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(29); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_T',1)); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([10 1000]); caxis([-1 +1]*0.15); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt T UMBC'])
colormap(cmap); ylim([1 1000])

figure(28); 
pcolor(xumbc.rlat,xumbc.pavgLAY,xumbc.trend_RH'); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_RH',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([10 1000]); caxis([-1 +1]*0.5); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH UMBC'])
colormap(cmap); ylim([100 1000])

figure(30); 
pcolor(xumbc.rlat,xumbc.pavgLAY,xumbc.trend_WV'); 
pcolor(xumbc.rlat,xumbc.pavgLAY,smoothn(xumbc.trend_WV',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
ylim([10 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac UMBC'])
colormap(cmap); ylim([100 1000])

%{
figure(06); aslprint([diroutQuick '/smoothed_stemp_rates.pdf'])
figure(28); aslprint([diroutQuick '/smoothed_rh_rates.pdf'])
figure(29); aslprint([diroutQuick '/smoothed_tz_rates.pdf'])
figure(30); aslprint([diroutQuick '/smoothed_wv_rates.pdf'])
%}
