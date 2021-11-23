addpath /asl/matlib/maps
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDo = -1;
if iDo > 0
  clear all
  
  figure(1);
  % Want a 3x3 plot
  ta = tiledlayout(3,3);
  % nexttile says go to the obvious, save it’s handle
  
  xx = 1 : 100;
  phase = randn(1:9)*180;
  odd_ndesc(1:9,1:100) = 0;
  even_ndesc(1:9,1:100) = 0;
  for ii = 1 : 9
    odd_ndesc(ii,1:100)  = cos(xx*pi/180 + phase(ii));
    even_ndesc(ii,1:100) = sin(xx*pi/180 + phase(ii));
    even_ndesc(ii,1:100) = even_ndesc(ii,1:100) + odd_ndesc(ii,1:100);
  end
  
  for ifov = 1:9
     tafov(ifov) = nexttile;
     plot(xx,odd_ndesc(ifov,:)./max(odd_ndesc(ifov,:)),'linewidth',1);
     hold on;
     plot(xx,even_ndesc(ifov,:)./max(even_ndesc(ifov,:)),'linewidth',1);
     grid;
  end
  % Get rid of all extra space I can
  ta.Padding = 'none';
  ta.TileSpacing = 'none';
  % Remove all xtick labels except for bottom row
  for i=1:6
     tafov(i).XTickLabel = '';
  end
  % Remove all ytick labels except for 1st column
  for i = [2 3 5 6 8 9]
     tafov(i).YTickLabel = '';
  end
  % Put in xlabel and ylable in the “middle”
  tafov(4).YLabel.String = 'PDF';
  tafov(8).XLabel.String = '\Delta \nu';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = -1;
if iDo > 0
  clear all
  
  figure(2);
  umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_o3_zonal_trends.mat');
  airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/airsL3_o3_zonal_trends.mat');
  era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_o3_zonal_trends.mat');
  
  load('llsmap5.mat');
  if length(llsmap5) == 64
    %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
    llsmap5 = llsmap5(2:64,:);  
  end

  %%want a 1x3 tiled layout
  ta = tiledlayout(1,3);
  
  cmap = llsmap5;
  cmax = 0.01;
  %cmax = 0.1;
  
  tafov(1) = nexttile;
  pcolor(umbc.rlat,umbc.p97,umbc.data);  title('UMBC'); 
    set(gca,'ydir','reverse');   set(gca,'yscale','log');  ylim([1 100]);  colormap(cmap); caxis([-cmax +cmax]); shading interp
  
  tafov(2) = nexttile;
  pcolor(airsL3.rlat,airsL3.plevs,airsL3.data); title('AIRSL3'); 
    set(gca,'ydir','reverse');   set(gca,'yscale','log');  ylim([1 100]);  colormap(cmap); caxis([-cmax +cmax]); shading interp
  
  tafov(3) = nexttile;
  pcolor(era.rlat,era.plays,era.data/10); title('ERA/10'); 
    set(gca,'ydir','reverse');   set(gca,'yscale','log');  ylim([1 100]);  colormap(cmap); caxis([-cmax +cmax]); shading interp
  
  % Get rid of all extra space I can
  ta.Padding = 'none';
  ta.TileSpacing = 'none';
  
  % Remove all ytick labels except for 1st column
  for ii = [2 3]
     tafov(ii).YTickLabel = '';
  end
  
  % Put in xlabel and ylable in the “middle”
  tafov(1).YLabel.String = 'P(mb)';
  tafov(2).XLabel.String = 'Latitude';
  
  %for ii = 1 : 3
    %set(gca,'ydir','reverse');   set(gca,'yscale','log');  ylim([1 100]);  colormap(cmap); caxis([-cmax +cmax]);
  %  tafov(ii).Ydir = 'reverse';
  %end
  
  cb = colorbar;
  cb.Layout.Tile = 'south';
  %cb.Layout.Tile = 'east';
  %boo = get(gcf,'position')
  %boo(1) = boo(1)+100; boo(3) = boo(3)-2*100;
  %set(cb,'position',boo)
  %boo = get(gcf,'position')
  %where the position arguments are [xposition yposition width height].
  
  %% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_zonal_o3_trends.pdf');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

figure(3);
umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_o3_zonal_trends.mat');
airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/airsL3_o3_zonal_trends.mat');
era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_o3_zonal_trends.mat');

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);  
end

%%want a 1x3 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
%cmax = 0.1;

tafov(1) = nexttile;
pcolor(umbc.rlat,umbc.p97,umbc.data);  %title('UMBC','Fontsize',12); 
shading interp;

tafov(2) = nexttile;
pcolor(airsL3.rlat,airsL3.plevs,airsL3.data); %title('AIRSL3','Fontsize',12);
shading interp;

tafov(3) = nexttile;
pcolor(era.rlat,era.plays,era.data/10); %title('ERA/10','Fontsize',12); 
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');   
set(tafov,'yscale','log');  
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[1 100]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)';
tafov(2).XLabel.String = 'Latitude';
tafov(1).YLabel.FontSize = 10;
tafov(2).XLabel.FontSize = 10;
%xlabel(tafov(2),'Latitude','fontsize',10)
%ylabel(tafov(1),'P(mb)','fontsize',10)

%% put titles
title(tafov(1),'UMBC','Fontsize',12); 
title(tafov(2),'AIRS L3','Fontsize',12); 
title(tafov(3),'ERA','Fontsize',12); 

%% assign colorbar to one tile
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';

%{
boo = get(cbh,'position')
%where the position arguments are [xposition yposition width height].
%}

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_zonal_o3_trends.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

figure(4); clf
umbc   = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_trendsonlyCLR/umbc_o3_zonal_trends.mat');
airsL3 = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/airsL3_o3_zonal_trends.mat');
era    = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/era_o3_zonal_trends.mat');

load('llsmap5.mat');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);  
end

%%want a 1x3 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.01;
%cmax = 0.1;

tafov(1) = nexttile;
aslmapSergio(umbc.rlat65,umbc.rlon73,umbc.dataMap, [-90 +90],[-180 +180]); 

tafov(2) = nexttile;
aslmapSergio(airsL3.rlat65,airsL3.rlon73,airsL3.dataMap, [-90 +90],[-180 +180]); 

tafov(3) = nexttile;
aslmapSergio(era.rlat65,era.rlon73,era.dataMap/10, [-90 +90],[-180 +180]); 

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

%% put titles
title(tafov(1),'UMBC','Fontsize',12); 
title(tafov(2),'AIRS L3','Fontsize',12); 
title(tafov(3),'ERA','Fontsize',12); 

%% assign colorbar to one tile
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';

%% aslprint
%% aslprint_asis('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/all_global_o3_025mb.pdf');
