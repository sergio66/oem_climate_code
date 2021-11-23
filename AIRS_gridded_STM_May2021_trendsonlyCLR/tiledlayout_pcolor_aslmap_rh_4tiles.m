addpath /asl/matlib/maps
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

clear all

close all
figure(1); scr_siz = get(gcf); a0 = scr_siz.Position;

figure(2); set(gcf, 'Position',  [100, 100,  560*1.25, 420]);
figure(3); set(gcf, 'Position',  [850, 100,  560*1.25, 420]);

figure(2);
umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_RH_zonal_trends.mat');
airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/airsL3_rh_zonal_trends.mat');
era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_rh_zonal_trends.mat');
merra  = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/merra_rh_zonal_trends.mat');

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);  
end

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
cmax = 0.1;
cmax = 0.5;

tafov(1) = nexttile;
pcolor(umbc.rlat,umbc.p97,umbc.data);  %title('UMBC','Fontsize',12); 
shading interp;

tafov(2) = nexttile;
pcolor(airsL3.rlat,airsL3.hlevs,airsL3.data); %title('AIRSL3','Fontsize',12);
shading interp;

tafov(3) = nexttile;
pcolor(era.rlat,era.plays,era.data); %title('ERA','Fontsize',12); 
shading interp;

tafov(4) = nexttile;
pcolor(merra.rlat,merra.plevsM,merra.data'); %title('MERRA','Fontsize',12); 
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

set(tafov,'ydir','reverse');   
set(tafov,'yscale','log');  
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[100 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

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
%title(tafov(1),'UMBC','Fontsize',12); 
%title(tafov(2),'AIRS L3','Fontsize',12); 
%title(tafov(3),'ERA','Fontsize',12); 
%title(tafov(4),'ERA','Fontsize',12); 
%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(4),'MERRA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%{
boo = get(cbh,'position')
%where the position arguments are [xposition yposition width height].
%}

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_zonal_rh_trends.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

figure(3); clf
umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_RH_zonal_trends.mat');
airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/airsL3_rh_zonal_trends.mat');
era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_rh_zonal_trends.mat');
merra  = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/merra_rh_zonal_trends.mat');

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);  
end

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
cmax = 0.1;
cmax = 0.5;

tafov(1) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,umbc.dataMap, [-90 +90],[-180 +180]); 

tafov(2) = nexttile;
aslmapSergio(airsL3.rlat65,airsL3.rlon73,airsL3.dataMap, [-90 +90],[-180 +180]); 

tafov(3) = nexttile;
aslmapSergio(era.rlat65,era.rlon73,era.dataMap, [-90 +90],[-180 +180]); 

tafov(4) = nexttile;
aslmapSergio(merra.rlat65,merra.rlon73,merra.dataMap, [-90 +90],[-180 +180]); 

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(4),'MERRA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%title(tafov(1),'UMBC','Fontsize',12); 
%title(tafov(2),'AIRS L3','Fontsize',12); 
%title(tafov(3),'ERA','Fontsize',12); 
%title(tafov(4),'MERRA','Fontsize',12); 

%% assign colorbar to one tile
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_global_rh_500mb.pdf');
