clear *bias_* 

if iFrac == +1
  %% Night
  ii = 2;
elseif iFrac == -1
  %% Day
  ii = 3;
elseif iFrac == 0
  %% Avg
  ii = 4;
elseif iFrac == -10
  %% Diff
  ii = 5;
end
ii = ii+3;

ii = ii + 1; figure(ii); clf; 
  ta = tiledlayout(2,3);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tafov(1) = nexttile; 
allbias_11      = nansum(obsrates.rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_11 = nanmean(obsrates.rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_11   = nanmean(obsrates.rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_11    = nanmean(obsrates.rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_11,'k.-',f,polarbias_11,'b',f,midlatbias_11,'r',f,tropicalbias_11,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('AIRS L1C');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(2) = nexttile; 
allbias_12      = nansum(era5.era5_spectral_rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_12 = nanmean(era5.era5_spectral_rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_12   = nanmean(era5.era5_spectral_rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_12    = nanmean(era5.era5_spectral_rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_12,'k.-',f,polarbias_12,'b',f,midlatbias_12,'r',f,tropicalbias_12,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('ERA5');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(3) = nexttile; 
allbias_13      = nansum(airsL3.airsL3_spectral_rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_13 = nanmean(airsL3.airsL3_spectral_rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_13   = nanmean(airsL3.airsL3_spectral_rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_13    = nanmean(airsL3.airsL3_spectral_rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_13,'k.-',f,polarbias_13,'b',f,midlatbias_13,'r',f,tropicalbias_13,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('AIRS L3');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(4) = nexttile; 
allbias_13      = nansum(umbcL3.umbcL3_spectral_rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_13 = nanmean(umbcL3.umbcL3_spectral_rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_13   = nanmean(umbcL3.umbcL3_spectral_rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_13    = nanmean(umbcL3.umbcL3_spectral_rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_13,'k.-',f,polarbias_13,'b',f,midlatbias_13,'r',f,tropicalbias_13,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('AIRS\_RT');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(5) = nexttile; 
allbias_22      = nansum(merra2.merra2_spectral_rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_22 = nanmean(merra2.merra2_spectral_rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_22   = nanmean(merra2.merra2_spectral_rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_22    = nanmean(merra2.merra2_spectral_rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_22,'k.-',f,polarbias_22,'b',f,midlatbias_22,'r',f,tropicalbias_22,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('MERRA2');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(6) = nexttile; 
allbias_23      = nansum(climcapsL3.climcapsL3_spectral_rates .* muYY2645,2) ./ nansum(muYY2645,2);
tropicalbias_23 = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(abs(p.rlat) <= 30)),2);
midlatbias_23   = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(abs(p.rlat) > 30 & abs(p.rlat) < 60)),2);
polarbias_23    = nanmean(climcapsL3.climcapsL3_spectral_rates(:,find(abs(p.rlat) > 60)),2);
  
plot(f,allbias_23,'k.-',f,polarbias_23,'b',f,midlatbias_23,'r',f,tropicalbias_23,'g');
 xlim([645 1608]); xticks([750:250:1500]); ylim([-0.01 +0.005]*10); 
 if iFrac == -10
   ylim([-0.01 +0.01]);
 end
 ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); t = title('CLIMCAPS L3');
 t.FontSize = 10; t.FontWeight = 'normal'; 
 xlabel('Wavenumber cm-1','Fontsize',10); ylabel('dBT/dt K/yr','Fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 3 5 6]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end
% Remove all xtick labels except for 2nd row
for ii = [1 2 3]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end
