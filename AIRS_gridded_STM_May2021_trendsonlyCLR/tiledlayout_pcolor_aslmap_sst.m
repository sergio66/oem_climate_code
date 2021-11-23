addpath /asl/matlib/maps
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

clear all

%{
-rw-rw-r-- 1 sergio pi_strow     73837 May 12 21:06 ../AIRS_gridded_May2021_trendsonlyCLR/umbc_trends.mat
-rw-rw-r-- 1 sergio pi_strow    255231 May 12 20:27 ../AIRS_gridded_May2021_trendsonlyCLR/giss_trends.mat
-rw-rw-r-- 1 sergio pi_strow   1099034 May 12 20:26 ../AIRS_gridded_May2021_trendsonlyCLR/airsL3_trends.mat
-rw-rw-r-- 1 sergio pi_strow 20604849 May 15 08:40 ../AIRS_gridded_May2021_anomalyonlyCLR/era_trends_StartSept2002_v2.mat
%}

figure(1); clf;
figure(2); set(gcf, 'Position',  [100, 100,  560*1.25, 420]);
figure(3); set(gcf, 'Position',  [100, 100,  560*1.25, 420]);

umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_trends.mat');
airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/airsL3_trends.mat');
era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_trends_StartSept2002_v2.mat');
giss   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/giss_trends.mat');

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);

%%want a 2x2 tiled layout
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
cmax = 0.15;

tafov(1) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(umbc.umbc_st_trend4608',1), [-90 +90],[-180 +180]); 

tafov(2) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(airsL3.airsL3_trend4608',1), [-90 +90],[-180 +180]); 

tafov(3) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(reshape(era.skt_trends,72,64)',1), [-90 +90],[-180 +180]); 

tafov(4) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(giss.giss_trend4608',1), [-90 +90],[-180 +180]); 

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

%% put titles
title(tafov(1),'  UMBC ', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(3),'  ERA  ', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(4),'  GISS ', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(1),'UMBC','Fontsize',12); 
%title(tafov(2),'AIRS L3','Fontsize',12); 
%title(tafov(3),'ERA','Fontsize',12); 
%title(tafov(4),'GISS','Fontsize',12); 

%% assign colorbar to one tile
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_global_sst.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);

%%want a 1x2 tiled layout
ta = tiledlayout(1,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
cmax = 0.15;

tafov(1) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(umbc.airs_quantile16_bt1231_trend4608',1), [-90 +90],[-180 +180]); 

tafov(2) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,smoothn(umbc.umbc_st_trend4608',1), [-90 +90],[-180 +180]); 

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

%% put titles
title(tafov(1),'BT1231', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(2),' UMBC ', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(1),'BT1231','Fontsize',12); 
%title(tafov(2),'UMBC','Fontsize',12); 

%% assign colorbar to one tile
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
cbh = colorbar(tafov(1));
cbh.Layout.Tile = 'south';

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_global_bt1231_sst.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
