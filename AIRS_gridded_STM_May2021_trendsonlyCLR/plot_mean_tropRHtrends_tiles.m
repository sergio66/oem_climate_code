iFig = 42;

clear ta taHandle tafov;

figure(42); clf
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 1/5;

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(1) = nexttile;
aslmapSergio(rlat65,rlon73,smoothn((reshape(meanAIRSrhtrend_trop,72,64)'),1), [-90 +90],[-180 +180]);   title('UMBC'); caxis([-1 +1]/5); colormap(cmap)

tafov(2) = nexttile;
aslmapSergio(rlat65,rlon73,smoothn((reshape(meanAIRSL3rhtrend_trop,72,64)'),1), [-90 +90],[-180 +180]); title('AIRS L3'); caxis([-1 +1]/5); colormap(cmap)

tafov(3) = nexttile;
aslmapSergio(rlat65,rlon73,smoothn((reshape(meanERA5rhtrend_trop,72,64)'),1), [-90 +90],[-180 +180]);   title('ERA5'); caxis([-1 +1]/5); colormap(cmap)

tafov(4) = nexttile;
aslmapSergio(rlat65,rlon73,smoothn((reshape(meanCMIP6rhtrend_trop,72,64)'),1), [-90 +90],[-180 +180]);  title('CMIP6'); caxis([-1 +1]/5); colormap(cmap)

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

set(tafov,'colormap',cmap,'CLim',[-cmax*1.01 +cmax*1.01]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';
