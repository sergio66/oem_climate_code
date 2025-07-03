function [] = profile_plots_2x1x2tiledlayout_tall(x,y,z11,z12,z13,z21,z22,iFig,plotoptions);

%% zXY = rowX, column Y   X=1-3,Y=1-2
%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/tiledlayout_pcolor_aslmap_rh_4tiles_eradonebyCLH.m

%{
EXAMPLE /home/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_6panel_zonaltrends.m

[h,ha,p,pa] = rtpread('../summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp');
z11 = p.ptemp;
z11 = reshape(p.ptemp,101,72,64);
z11 = squeeze(nanmean(z11,2));
z12 = z11 * 2;
zz = z11 * 3;
z21 = z11 * 4;
z22 = z11 * 5;
y = p.plevs(:,2300);
x = p.rlat; x = reshape(x,72,64); x = nanmean(x,1);
profile_plots_2x1x2tiledlayout_wide(x,y,z11,z12,zz,z21,z22,2);
profile_plots_2x1x2tiledlayout_tall(x,y,z11,z12,zz,z21,z22,3);
%}

% addpath /asl/matlib/maps
% addpath /asl/matlib/plotutils
% addpath /home/sergio/MATLABCODE
% addpath /home/sergio/MATLABCODE/PLOTTER
% addpath /home/sergio/MATLABCODE/COLORMAP
% addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% plotoptions can have following fields
%%   cx                    = [caxis1 caxis2]
%%   cmap                  = colormap
%%   yLinearOrLog          = +1/-1
%%   yReverseDir           = +1/-1
%%   x/yLimits             = [-90 +90] [1 1000];
%%   str11,str12,str3,str4 = give the subplots some flair
%%   xstr,ystr             = xlabel, ylabel
%%   xaxis_metric          = 'linear' or 'sine'
%%   maintitle             = maintitle
%%   smooth                  = +1 for smothn, +0.5 for smoothdata(movmean,5),0 for no smooth

%{
%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/tiledlayout_pcolor_aslmap_sst.m
load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end
%}

jet64 = jet(64);
jett = jet64; jett(1,:) = 1;

cx1 = min([nanmin(z11(:)) nanmin(z12(:)) nanmin(z13(:)) nanmin(z21(:)) nanmin(z22(:))]);
cx2 = max([nanmax(z11(:)) nanmax(z12(:)) nanmax(z13(:)) nanmax(z21(:)) nanmax(z22(:))]);

xlimits = [min(x) max(x)];  
ylimits = [min(y) max(y)];  
plotoptions0.cx           = [cx1 cx2];
plotoptions0.cmap   = jett;
plotoptions0.yLinearOrLog = +1;
plotoptions0.yReverseDir  = -1;
plotoptions0.xLimits      = xlimits;
plotoptions0.yLimits      = ylimits;
plotoptions0.str11 = '(A)';
plotoptions0.str12 = '(B)';
plotoptions0.str13 = '(C)';
plotoptions0.str21 = '(D)';
plotoptions0.str22 = '(E)';
plotoptions0.xstr = 'Xstr';
plotoptions0.ystr = 'Ystr';
plotoptions0.cstr = 'Cstr';
plotoptions0.barstr = 'barstr';
plotoptions0.xaxis_metric = 'linear';
plotoptions0.maintitle = 'Main Title';
plotoptions0.smooth    = +0; %% no smoothing
plotoptions0.smooth    = +1; %% smooth Matlab, 1 pt window
plotoptions0.smooth    = +1; %% smoothn <default>

if nargin == 7
  iFig = 1;
  plotoptions = plotoptions0;
elseif nargin == 8
  plotoptions = plotoptions0;
end

if ~isfield(plotoptions,'smooth')
  plotoptions.smooth = plotoptions0.smooth;
end
if ~isfield(plotoptions,'xaxis_metric')
  plotoptions.xaxis_metric = plotoptions0.xaxis_metric;
else
  if strfind(plotoptions.xaxis_metric,'sine')
    xlimits = [-1 +1];  
    plotoptions0.xLimits = xlimits;
    plotoptions.xLimits  = xlimits;
  end
end
if ~isfield(plotoptions,'cx')
  plotoptions.cx = plotoptions0.cx;
end
if ~isfield(plotoptions,'cstr')
  plotoptions.cstr = plotoptions0.cstr;
end
if ~isfield(plotoptions,'cmap')
  plotoptions.cmap = plotoptions0.cmap;
end
if ~isfield(plotoptions,'yLinearOrLog')
  plotoptions.yLinearOrLog = plotoptions0.yLinearOrLog;
end
if ~isfield(plotoptions,'yReverseDir')
  plotoptions.yReverseDir = plotoptions0.yReverseDir;
end
if ~isfield(plotoptions,'xLimits')
  plotoptions.xLimits = plotoptions0.xLimits;
end
if ~isfield(plotoptions,'yLimits')
  plotoptions.yLimits = plotoptions0.yLimits;
end
if ~isfield(plotoptions,'str11')
  plotoptions.str11 = plotoptions0.str11;
end
if ~isfield(plotoptions,'str12')
  plotoptions.str12 = plotoptions0.str12;
end
if ~isfield(plotoptions,'str21')
  plotoptions.str21 = plotoptions0.str21;
end
if ~isfield(plotoptions,'str22')
  plotoptions.str22 = plotoptions0.str22;
end
if ~isfield(plotoptions,'str13')
  plotoptions.str13 = plotoptions0.str13;
end
if ~isfield(plotoptions,'xstr')
  plotoptions.xstr = plotoptions0.xstr;
end
if ~isfield(plotoptions,'ystr')
  plotoptions.ystr = plotoptions0.ystr;
end
if ~isfield(plotoptions,'barstr')
  plotoptions.barstr = plotoptions0.barstr;
end
if ~isfield(plotoptions,'maintitle')
  plotoptions.maintitle = plotoptions0.maintitle;
end

fig = figure(iFig); clf; 
figure(iFig); scr_siz = get(gcf); a0 = scr_siz.Position;
figure(iFig); %set(gcf, 'Position',  [100, 100,  560*1.25, 420]);

xlimits = plotoptions.xLimits;
ylimits = plotoptions.yLimits;

%%want a 2x1x2 tiled layout
ta = tiledlayout(3,4,'TileSpacing','None', 'Padding','None');
ta = tiledlayout(6,6);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cx   = plotoptions.cx;
cmap = plotoptions.cmap;
figure(iFig);

if strfind(plotoptions.xaxis_metric,'linear')
  disp('  xaxis will be linear in spacing (rlat)')
  iMetric = +1; %% default
elseif strfind(plotoptions.xaxis_metric,'sine')
  disp('  xaxis will be spaced by sin(rlat) so narrowest at poles')
  iMetric = +2;
else
  error('  xaxis needs sine or linear spacing (xaxis_metric)')
end

if plotoptions.smooth >= 1
  disp('using smoothn(z,1)')
elseif plotoptions.smooth == 0.5
  disp('using smooth(z,''movmean'',5)')
elseif plotoptions.smooth <= 0
  disp('no smoothing')
end

if plotoptions.smooth >= 1
  wah11 = smoothn(z11,1);
  wah12 = smoothn(z12,1);
  wah13 = smoothn(z13,1);
  wah21 = smoothn(z21,1);
  wah22 = smoothn(z22,1);
elseif plotoptions.smooth == 0.5
  wah11 = smoothdata(z11,"movmean",5);
  wah12 = smoothdata(z12,"movmean",5);
  wah13 = smoothdata(z13,"movmean",5);
  wah21 = smoothdata(z21,"movmean",5);
  wah22 = smoothdata(z22,"movmean",5);
elseif plotoptions.smooth <= 0
  wah11 = z11;
  wah12 = z12;
  wah13 = z13;
  wah21 = z21;
  wah22 = z22;
end

if iMetric == 1
  %% 2 + 2 = 4
  tafov(1) = nexttile([3,2]);
  pcolor(x,y,wah11); shading interp; colormap(plotoptions.cmap);
  
  tafov(2) = nexttile([3,2]);
  pcolor(x,y,wah12); shading interp; colormap(plotoptions.cmap);
  
  tafov(3) = nexttile([3,2]);
  pcolor(x,y,wah13); shading interp; colormap(plotoptions.cmap);

  nexttile([1,1])
  axis off
  
  tafov(4) = nexttile([3,2]);
  pcolor(x,y,wah21); shading interp; colormap(plotoptions.cmap);
  
  tafov(5) = nexttile([3,2]);
  pcolor(x,y,wah22); shading interp; colormap(plotoptions.cmap);
  
elseif iMetric == 2
  %% see MATLABCODE/PLOTTER/pcolor_sin
  xsin = sin(x*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

  tafov(1) = nexttile([3,2]);;
  pcolor(x,y,wah11); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  
  tafov(2) = nexttile([3,2]);
  pcolor(x,y,wah12); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  
  tafov(3) = nexttile([3,2]);
  pcolor(x,y,wah13); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])

  nexttile([1,1])
  axis off
  
  tafov(4) = nexttile([3,2]);
  pcolor(x,y,wah21); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  
  tafov(5) = nexttile([3,2]);
  pcolor(x,y,wah22); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iASLMAP = -1;
common_2x1x2_tall_tiled_layout
