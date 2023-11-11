ncYY = nansum(cosYY,2);
iOffSet = 39;
fprintf(1,'Figs %2i : %2i showing Biases : All lats (Obs-Cal) cosine averaged, Tropical, MidLats, Polar \n',iOffSet+1,iOffSet+4)
figure(iOffSet+1); clf;
plot(fchanx,nansum(cosYY.*obsrates.rates - cosYY.*obsrates.rates,2)./ncYY,'k.-',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYY.*obsrates.rates - cosYY.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',......
     fchanx,nansum(cosYY.*obsrates.rates - cosYY.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*7.5)
 plotaxis2;
 title('All lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYPolar = cos(YY*pi/180)';
cosYYPolar = ones(size(YY))';
cosYYPolar(abs(YY) <= 60) = 0;
cosYYPolar = ones(2645,1) * cosYYPolar;
ncYY = nansum(cosYYPolar,2);
figure(iOffSet+2); clf;
plot(fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*obsrates.rates,2)./ncYY,'k.-',fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYPolar.*obsrates.rates - cosYYPolar.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*7.5)
 plotaxis2;
 title('Polar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYMidLats = cos(YY*pi/180)';
cosYYMidLats = ones(size(YY))';
cosYYMidLats(abs(YY) > 60 | abs(YY) <= 30) = 0;
cosYYMidLats = ones(2645,1) * cosYYMidLats;
ncYY = nansum(cosYYMidLats,2);
figure(iOffSet+3); clf;
plot(fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*obsrates.rates,2)./ncYY,'k.-',fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYMidLats.*obsrates.rates - cosYYMidLats.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*7.5)
 plotaxis2;
 title('MidLatitudes : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYTropics = cos(YY*pi/180)';
cosYYTropics = ones(size(YY))';
cosYYTropics(abs(YY) > 30) = 0;
cosYYTropics = ones(2645,1) * cosYYTropics;
ncYY = nansum(cosYYTropics,2);
figure(iOffSet+4); clf;
plot(fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*obsrates.rates,2)./ncYY,'k.-',fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYTropics.*obsrates.rates - cosYYTropics.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*7.5)
 plotaxis2;
 title('Tropics lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

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
