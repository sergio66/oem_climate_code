save Output/strow.mat latx plays water_airs* temp_airs* coljacs_airs* 

addpath /asl/matlib/plotutils
figure(1); aslprint('Output/allsky_bias_latbins.pdf')
figure(4); aslprint('Output/allsky_tracegas_latbins.pdf')
figure(5); aslprint('Output/allsky_cloud_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf
  pcolor(latx,(plays),water_airs'); shading flat; title('d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(3); clf
  pcolor(latx,(plays),temp_airs'); shading flat; title('d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/allsky_WV_z_latbins.pdf')
figure(3); aslprint('Output/allsky_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
  pcolor(latx,Qlevs,double(l3waterrate')); shading flat; title('L3 d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(3); clf
  pcolor(latx,(Tlevs),double(l3temprate')); shading flat; title('L3 d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/L3_WV_z_latbins.pdf')
figure(3); aslprint('Output/L3_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
  pcolor(latx,plays,double(gas1_rate_cloud')); shading flat; title('ERA d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(3); clf
  pcolor(latx,plays,double(ptemp_rate_cloud')); shading flat; title('ERA d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/era_WV_z_latbins.pdf')
figure(3); aslprint('Output/era_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%

waterrate_cloud_ak(36,:) = NaN;
figure(2); clf
  pcolor(latx,plays,double(waterrate_cloud_ak')); shading flat; title('ERA*AK d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

temprate_cloud_ak(36,:) = NaN;
figure(3); clf
  pcolor(latx,plays,double(temprate_cloud_ak')); shading flat; title('ERA*AK d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/era_ak_WV_z_latbins.pdf')
figure(3); aslprint('Output/era_ak_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
  pcolor(latx,plays,squeeze(merrarates.water_allpars(:,:,2)')); shading flat; title('MERRA d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(3); clf
  pcolor(latx,plays,squeeze(merrarates.ptemp_allpars(:,:,2)')); shading flat; title('MERRA d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/merra_WV_z_latbins.pdf')
figure(3); aslprint('Output/merra_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%

waterrate_merra_ak(36,:) = NaN;
figure(2); clf
  pcolor(latx,plays,double(waterrate_merra_ak')); shading flat; title('MERRA*AK d(W)/dt frac/yr')
  colormap(color5.llsmap5); caxis([-0.05 +0.05]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

temprate_merra_ak(36,:) = NaN;
figure(3); clf
  pcolor(latx,plays,double(temprate_merra_ak')); shading flat; title('MERRA*AK d(T)/dt K/yr')
  colormap(color5.llsmap5); caxis([-0.15 +0.15]); colorbar
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  axis([-90 +90 9 1000])                                                   
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})

figure(2); aslprint('Output/merra_ak_WV_z_latbins.pdf')
figure(3); aslprint('Output/merra_ak_T_z_latbins.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

