figure(4); clf; semilogy(nanmean(resultsT,1),pavg,'bx-',nanmean(resultsTunc,1)/sqrt(4608),pavg,'c.-');   set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('T and \sigma T')
figure(5); clf; semilogy(nanmean(resultsWV,1),pavg,'bx-',nanmean(resultsWVunc,1)/sqrt(4608),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('WV and \sigma WV')
figure(6); clf; semilogy(nanmean(resultsO3,1),pavg,'bx-',nanmean(resultsO3unc,1)/sqrt(4608),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('O3 and \sigma O3')

figure(4); clf; semilogy(nanmean(resultsT,1),pavg,'bx-',nanmean(resultsTunc,1)/sqrt(72),pavg,'c.-');   set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('T and \sigma T')
figure(5); clf; semilogy(nanmean(resultsWV,1),pavg,'bx-',nanmean(resultsWVunc,1)/sqrt(72),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('WV and \sigma WV')
figure(6); clf; semilogy(nanmean(resultsO3,1),pavg,'bx-',nanmean(resultsO3unc,1)/sqrt(72),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('O3 and \sigma O3')

if exist('plays') & exist('era5')
  figure(4); clf; semilogy(nanmean(resultsT,1),pavg,'bx-',nanmean(resultsTunc,1)/sqrt(72),pavg,'b',nanmean(era5.trend_ptemp,2),plays,'rx-',nanmean(era5.trend_ptemp_err,2)/sqrt(72),plays,'r');   
    set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('T and \sigma T')
    hl = legend('mean UMBC','mean \sigma UMBC)/sqrt(72)','mean ERA5','mean \sigma ERA5)/sqrt(72)','location','best','fontsize',10);
  figure(5); clf; semilogy(nanmean(resultsWV,1),pavg,'bx-',nanmean(resultsWVunc,1)/sqrt(72),pavg,'b',nanmean(era5.trend_gas_1,2),plays,'rx-',nanmean(era5.trend_gas_1_err,2)/sqrt(72),plays,'r');   
    set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('WV and \sigma WV')
    hl = legend('mean UMBC','mean \sigma UMBC)/sqrt(72)','mean ERA5','mean \sigma ERA5)/sqrt(72)','location','best','fontsize',10);
  figure(6); clf; semilogy(nanmean(resultsO3,1),pavg,'bx-',nanmean(resultsO3unc,1)/sqrt(72),pavg,'b',nanmean(era5.trend_gas_3,2),plays,'rx-',nanmean(era5.trend_gas_3_err,2)/sqrt(72),plays,'r');   
    set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('O3 and \sigma O3')
    hl = legend('mean UMBC','mean \sigma UMBC)/sqrt(72)','mean ERA5','mean \sigma ERA5)/sqrt(72)','location','best','fontsize',10);
end

figure(1); clf;
  wah = abs(resultsT'./resultsTunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); title('S/N : T/\sigma T'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp; 
  set(gca,'yscale','log'); ylim([1 1000])
figure(2); clf;
  wah = abs(resultsWV'./resultsWVunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); title('S/N : WV/\sigma WV'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp
figure(3); clf; 
  wah = abs(resultsO3'./resultsO3unc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); title('S/N : O3/\sigma O3'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); colorbar; colormap jet; shading interp

