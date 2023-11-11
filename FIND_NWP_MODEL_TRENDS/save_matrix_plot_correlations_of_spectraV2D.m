iOffSet = 38;
fprintf(1,'Figs %2i : %2i showing Biases : (Obs-Cal) cosine averaged : Land/Ocean : All_Lats, Tropical, MidLats, Polar \n',iOffSet+1,iOffSet+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncYYA = nansum(cosYY,2);

cosYYPolar = cos(YY*pi/180)';
cosYYPolar = ones(size(YY))';
cosYYPolar(abs(YY) <= 60) = 0;
cosYYPolar = ones(2645,1) * cosYYPolar;
ncYYP = nansum(cosYYPolar,2);

cosYYMidLats = cos(YY*pi/180)';
cosYYMidLats = ones(size(YY))';
cosYYMidLats(abs(YY) > 60 | abs(YY) <= 30) = 0;
cosYYMidLats = ones(2645,1) * cosYYMidLats;
ncYYM = nansum(cosYYMidLats,2);

cosYYTropics = cos(YY*pi/180)';
cosYYTropics = ones(size(YY))';
cosYYTropics(abs(YY) > 30) = 0;
cosYYTropics = ones(2645,1) * cosYYTropics;
ncYYT = nansum(cosYYTropics,2);

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iOffSet+1); clf;

%%want a 2x3 tiled layout
ta = tiledlayout(2,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = jet; 
dxyzaspect = [1.5 1.75 1]; %% not bad
dxyzaspect = [2 1.75 1]; %% not bad

smoothPts = 1;
smoothPts = 5;
smoothPts = 3;

tafov(1) = nexttile;
plot(fchanx,nansum(cosYY.*obsrates.rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*obsrates.rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*obsrates.rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*obsrates.rates,2)./ncYYT,'g',...
     'linewidth',0.5);
plot(fchanx,smooth(nansum(cosYY.*obsrates.rates,2)./ncYYA,smoothPts),'k.-',fchanx,smooth(nansum(cosYYPolar.*obsrates.rates,2)./ncYYP,smoothPts),'b',...
     fchanx,smooth(nansum(cosYYMidLats.*obsrates.rates,2)./ncYYM,smoothPts),'r',fchanx,smooth(nansum(cosYYTropics.*obsrates.rates,2)./ncYYT,smoothPts),'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('Obs')
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

tafov(2) = nexttile;
plot(fchanx,nansum(cosYY.*era5.era5_spectral_rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*era5.era5_spectral_rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*era5.era5_spectral_rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*era5.era5_spectral_rates,2)./ncYYT,'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('ERA5')
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

tafov(3) = nexttile;
plot(fchanx,nansum(cosYY.*airsL3.airsL3_spectral_rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*airsL3.airsL3_spectral_rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*airsL3.airsL3_spectral_rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*airsL3.airsL3_spectral_rates,2)./ncYYT,'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('AIRS L3')
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

tafov(4) = nexttile;
plot(fchanx,nansum(cosYY.*umbcL3.umbcL3_spectral_rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*umbcL3.umbcL3_spectral_rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*umbcL3.umbcL3_spectral_rates,2)./ncYYT,'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('UMBC')
 ylabel('dBT/dt K/yr')

tafov(5) = nexttile;
plot(fchanx,nansum(cosYY.*merra2.merra2_spectral_rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*merra2.merra2_spectral_rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*merra2.merra2_spectral_rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*merra2.merra2_spectral_rates,2)./ncYYT,'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('MERRA2')
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

tafov(6) = nexttile;
plot(fchanx,nansum(cosYY.*climcapsL3.climcapsL3_spectral_rates,2)./ncYYA,'k.-',fchanx,nansum(cosYYPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYYP,'b',...
     fchanx,nansum(cosYYMidLats.*climcapsL3.climcapsL3_spectral_rates,2)./ncYYM,'r',fchanx,nansum(cosYYTropics.*climcapsL3.climcapsL3_spectral_rates,2)./ncYYT,'g',...
     'linewidth',0.5);
 xlim([645 1608]); ylim([-0.01 +0.005]*10)
 plotaxis2;
 hl = legend('All Lats','Polar','MidLats','Tropics','location','best','fontsize',8); title('CLIMCAPS L3')
 ylabel('dBT/dt K/yr')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strAll  = 'All (600-1640 cm^{-1})'; 
str15um = '15um (640-800 cm^{-1})';
strWin  = 'Window (800-960 cm^{-1})';
strWV   = 'WV 6.7 um (1370-1620 cm^{-1})';

strAll  = 'All (LW/MW)';
str15um = 'CO2 15um';
strWin  = 'Window';
strWV   = 'WV 6.7 um';

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(70); clf; 

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = jet; 
dxyzaspect = [1.5 1.75 1]; %% not bad
dxyzaspect = [2 1.75 1]; %% not bad

colormap(cmap);

tafov(1) = nexttile; imagesc(Rallchans(2:6,:)'); title(strAll)
cx = caxis; cx = [0 1]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(2) = nexttile; imagesc(R15umchans(2:6,:)'); title(str15um)
cx = caxis; cx = [0 1]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(3) = nexttile;imagesc(Rwinchans(2:6,:)'); title(strWin)
cx = caxis; cx = [0 1]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(4) = nexttile; imagesc(RWVchans(2:6,:)'); title(strWV)
cx = caxis; cx = [0 1]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

set(tafov,'colormap',cmap,'CLim',cx);
%set(tafov,'Xlim',xlimits);
%set(tafov,'Ylim',ylimits);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

%% main title
%title(ta,'R^2')

%% assign colorbar to one tile
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(71); clf; 

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5; 
dxyzaspect = [1.5 1.75 1]; %% not bad
dxyzaspect = [2 1.75 1]; %% not bad

colormap(cmap);

tafov(1) = nexttile; imagesc(Mallchans(2:6,:)');
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.02 +0.02]; caxis(cx);
title(strAll)
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(2) = nexttile; imagesc(M15umchans(2:6,:)');
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.02 +0.02]; caxis(cx);
title(str15um)
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(3) = nexttile; imagesc(Mwinchans(2:6,:)');
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.02 +0.02]; caxis(cx);
title(strWin)
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(4) = nexttile; imagesc(MWVchans(2:6,:)');
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.02 +0.02]; caxis(cx);
title(strWV)
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

set(tafov,'colormap',cmap,'CLim',cx);
%set(tafov,'Xlim',xlimits);
%set(tafov,'Ylim',ylimits);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

%% main title
%title(ta,'\mu (K)')

%% assign colorbar to one tile
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(72); clf; 

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = jet; 
cmap = llsmap5; 
cmap = llsmap5(32:63,:);

dxyzaspect = [1.5 1.75 1]; %% not bad
dxyzaspect = [2 1.75 1]; %% not bad

colormap(cmap);

tafov(1) = nexttile; imagesc(Sallchans(2:6,:)'); title(strAll)
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(2) = nexttile; imagesc(S15umchans(2:6,:)'); title(str15um)
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(3) = nexttile; imagesc(Swinchans(2:6,:)'); title(strWin)
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

tafov(4) = nexttile; imagesc(SWVchans(2:6,:)'); title(strWV)
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

set(tafov,'colormap',cmap,'CLim',cx);
%set(tafov,'Xlim',xlimits);
%set(tafov,'Ylim',ylimits);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

%% main title
%title(ta,'\sigma (K)')

%% assign colorbar to one tile
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%
