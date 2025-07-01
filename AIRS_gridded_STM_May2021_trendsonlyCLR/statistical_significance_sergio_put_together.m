figure(69); clf;
  wah = abs(resultsTunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); set(gca,'ydir','reverse'); ylim([1 1000]); colormap jet; shading interp; % title('\sigma T'); ; plotaxis2;
  set(gca,'yscale','log'); ylim([10 1000])
  ylim([10 1000]); caxis([0 0.02]);  xlabel('Latitude [deg]'); ylabel('Pressure [mb]'); 
  % colorbar('vertical'); text(90,08,'K/yr'); set(gca,'fontsize',12)
  c = colorbar('horizontal'); text(-105,6000,'K/yr','fontsize',12); set(gca,'fontsize',sizefont)
  c.Ruler.Exponent = 0; 
  c.Ruler.TickLabelFormat = '%0.3f';

figure(70); clf;
  wah = abs(resultsWVunc'); wah = squeeze(nanmean(reshape(wah,length(pavg),72,64),2)); 
  pcolor(rlat,pavg,wah); set(gca,'ydir','reverse'); ylim([100 1000]); colormap jet; shading interp; % title('\sigma WV'); plotaxis2; 
  set(gca,'yscale','linear'); ylim([100 1000])
  ylim([100 1000]); caxis([0 0.003]); xlabel('Latitude [deg]'); ylabel('Pressure [mb]'); 
  % colorbar('vertical'); text(90,50,'1/yr'); set(gca,'fontsize',12)
  c = colorbar('horizontal'); text(-105,1350,'1/yr'); set(gca,'fontsize',sizefont)
  c.Ruler.Exponent = 0; 
  c.Ruler.TickLabelFormat = '%0.4f';

iOffset = 80;
statistical_significance_sergio
figure(69); set(gca,'fontsize',sizefont); text(-80,20, '(A)','fontsize',12); caxis([0 0.015]);
figure(70); set(gca,'fontsize',sizefont); text(-80,200,'(C)','fontsize',12); caxis([0 2e-3]);
figure(86); text(-105,6000,'K/yr','fontsize',12);; set(gca,'fontsize',sizefont); text(-80,20, '(B)','fontsize',12);  
figure(90); text(-105,1350,'1/yr');set(gca,'fontsize',sizefont); text(-80,200,'(D)','fontsize',12);  
