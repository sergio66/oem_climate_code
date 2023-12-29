clear *bias_* 

allbias_o = nansum(obsrates.rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(obsrates.rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(obsrates.rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(obsrates.rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(obsrates.rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(obsrates.rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(obsrates.rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(obsrates.rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(obsrates.rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(obsrates.rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(obsrates.rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(obsrates.rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('Obervations');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbias_o = nansum(umbcL3.umbcL3_spectral_rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(umbcL3.umbcL3_spectral_rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(umbcL3.umbcL3_spectral_rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(umbcL3.umbcL3_spectral_rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(umbcL3.umbcL3_spectral_rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(umbcL3.umbcL3_spectral_rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(umbcL3.umbcL3_spectral_rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(umbcL3.umbcL3_spectral_rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(umbcL3.umbcL3_spectral_rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(umbcL3.umbcL3_spectral_rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(umbcL3.umbcL3_spectral_rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(umbcL3.umbcL3_spectral_rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('CHIRP\_A');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbias_o = nansum(airsL3.airsL3_spectral_rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(airsL3.airsL3_spectral_rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(airsL3.airsL3_spectral_rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(airsL3.airsL3_spectral_rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(airsL3.airsL3_spectral_rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(airsL3.airsL3_spectral_rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(airsL3.airsL3_spectral_rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(airsL3.airsL3_spectral_rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(airsL3.airsL3_spectral_rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(airsL3.airsL3_spectral_rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(airsL3.airsL3_spectral_rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(airsL3.airsL3_spectral_rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('AIRS L3');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbias_o = nansum(climcapsL3.climcapsL3_spectral_rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(climcapsL3.climcapsL3_spectral_rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('CLIMCAPS L3');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbias_o = nansum(era5.era5_spectral_rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(era5.era5_spectral_rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(era5.era5_spectral_rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(era5.era5_spectral_rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(era5.era5_spectral_rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(era5.era5_spectral_rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(era5.era5_spectral_rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(era5.era5_spectral_rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(era5.era5_spectral_rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(era5.era5_spectral_rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(era5.era5_spectral_rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(era5.era5_spectral_rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('ERA5');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allbias_o = nansum(merra2.merra2_spectral_rates(:,oceanX) .* (ones(2645,1) * cos(YY(oceanX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(oceanX)'*pi/180),2);
tropicalbias_o = nanmean(merra2.merra2_spectral_rates(:,find(tropics_o == 1)),2);
smidlatbias_o = nanmean(merra2.merra2_spectral_rates(:,find(smidlats_o == 1)),2);
nmidlatbias_o = nanmean(merra2.merra2_spectral_rates(:,find(nmidlats_o == 1)),2);
spolarbias_o = nanmean(merra2.merra2_spectral_rates(:,find(spoles_o == 1)),2);
npolarbias_o = nanmean(merra2.merra2_spectral_rates(:,find(npoles_o == 1)),2);

allbias_l = nansum(merra2.merra2_spectral_rates(:,landX) .* (ones(2645,1) * cos(YY(landX)'*pi/180)),2) ./ nansum(ones(2645,1) * cos(YY(landX)'*pi/180),2);
tropicalbias_l = nanmean(merra2.merra2_spectral_rates(:,find(tropics_l == 1)),2);
smidlatbias_l = nanmean(merra2.merra2_spectral_rates(:,find(smidlats_l == 1)),2);
nmidlatbias_l = nanmean(merra2.merra2_spectral_rates(:,find(nmidlats_l == 1)),2);
spolarbias_l = nanmean(merra2.merra2_spectral_rates(:,find(spoles_l == 1)),2);
npolarbias_l = nanmean(merra2.merra2_spectral_rates(:,find(npoles_l == 1)),2);

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,allbias_o,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('OCEAN')
  title('MERRA2');

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; 
  plot(f,allbias_l,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','south','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.075 +0.075]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for iii = [1]
   tafov(iii).XTickLabel = '';
   tafov(iii).XLabel.String = [];
  end
  for iii = [1 2]
   tafov(iii).FontSize = 8;
  end

