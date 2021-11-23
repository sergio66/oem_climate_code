figure(32)
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.15;

tafov(1) = nexttile;
pcolor(rlat,pavgLAY(1:97,3000),smoothn(xdeltaTlat(:,1:97)',1)); %title('UMBC','Fontsize',12);
shading interp;

tafov(2) = nexttile;
boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = xmaskLFmatr'; end
junk = boo.* airsL3.thestats64x72.ptemprate; 
junk = squeeze(nanmean(junk,1))'; 
pcolor(rlat,airsL3.Tlevs,smoothn(junk(1:24,:),1)); %title('AIRSL3','Fontsize',12);
shading interp;

tafov(3) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.*reshape(era5.trend_ptemp,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); %title('ERA5','Fontsize',12); 
shading interp;

tafov(4) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.* reshape(cmip6.trend_ptemp,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1));  %title('MERRA','Fontsize',12); 
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

set(tafov,'ydir','reverse');   
set(tafov,'yscale','log');  
set(tafov,'colormap',cmap,'CLim',[-cmax*1.01 +cmax*1.01]);
set(tafov,'Ylim',[T_Ylim 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = 'P(mb)'; tafov(3).YLabel.FontSize = 10;
tafov(3).XLabel.String = 'Latitude'; tafov(3).XLabel.FontSize = 10;
tafov(4).XLabel.String = 'Latitude'; tafov(4).XLabel.FontSize = 10;
%xlabel(tafov(2),'Latitude','fontsize',10)
%ylabel(tafov(1),'P(mb)','fontsize',10)

%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(3),'ERA5', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(4),'CMIP6', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(33)
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.0125;

tafov(1) = nexttile;
pcolor(rlat,pavgLAY(1:97,3000),smoothn(xfracWVlat(:,1:97)',1)); %title('UMBC','Fontsize',12);
shading interp;

tafov(2) = nexttile;
boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
junk = boo.* airsL3.thestats64x72.waterrate; 
junk = squeeze(nanmean(junk,1))'; 
pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); %title('AIRSL3','Fontsize',12);
shading interp;

tafov(3) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.*reshape(era5.trend_gas_1,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); %title('ERA5','Fontsize',12); 
shading interp;

tafov(4) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.* reshape(cmip6.trend_gas_1,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1));  %title('MERRA','Fontsize',12); 
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

set(tafov,'ydir','reverse');   
set(tafov,'yscale','log');  
set(tafov,'colormap',cmap,'CLim',[-cmax*1.01 +cmax*1.01]);
set(tafov,'Ylim',[WV_Ylim 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = 'P(mb)'; tafov(3).YLabel.FontSize = 10;
tafov(3).XLabel.String = 'Latitude'; tafov(3).XLabel.FontSize = 10;
tafov(4).XLabel.String = 'Latitude'; tafov(4).XLabel.FontSize = 10;
%xlabel(tafov(2),'Latitude','fontsize',10)
%ylabel(tafov(1),'P(mb)','fontsize',10)

%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(3),'ERA5', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(4),'CMIP6', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(34)
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.15;
cmax = 0.25;

tafov(1) = nexttile;
pcolor(rlat,pavgLAY(1:97,3000),smoothn(xdeltaRHlat(:,1:97)',1)); %title('UMBC','Fontsize',12);
shading interp;

tafov(2) = nexttile;
boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
junk = boo.* airsL3.thestats64x72.RHrate; 
junk = squeeze(nanmean(junk,1))'; 
pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); %title('AIRSL3','Fontsize',12);
shading interp;

tafov(3) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.*reshape(era5.trend_RH,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); %title('ERA5','Fontsize',12); 
shading interp;

tafov(4) = nexttile;
boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
junk = boo.* reshape(cmip6.trend_RH,100,72,64); 
junk = squeeze(nanmean(junk,2)); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1));  %title('MERRA','Fontsize',12); 
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

set(tafov,'ydir','reverse');   
set(tafov,'yscale','log');  
set(tafov,'colormap',cmap,'CLim',[-cmax*1.01 +cmax*1.01]);
set(tafov,'Ylim',[WV_Ylim 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = 'P(mb)'; tafov(3).YLabel.FontSize = 10;
tafov(3).XLabel.String = 'Latitude'; tafov(3).XLabel.FontSize = 10;
tafov(4).XLabel.String = 'Latitude'; tafov(4).XLabel.FontSize = 10;
%xlabel(tafov(2),'Latitude','fontsize',10)
%ylabel(tafov(1),'P(mb)','fontsize',10)

%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(3),'ERA5', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(4),'CMIP6', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

