function [] = profile_figure6_wiley_2x1tiledlayout(x,y,T11,T21,W11,W21,iFig,plotoptionsT,plotoptionsW);

%% see /umbc/xfs2/strow/asl/s1/sergio/home/git/MATLABCODE_Git/PLOTTER/TILEDPLOTS/profile_plots_2x1tiledlayout.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/maps
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% plotoptions can have following fields
%%   cx                      = [caxis1 caxis2]
%%   cmap                    = colormap
%%   yLinearOrLog            = +1/-1
%%   yReverseDir             = +1/-1
%%   x/yLimits               = [-90 +90] [1 1000];
%%   str11,str21             = give the subplots some flair
%%   xstr,ystr               = xlabel, ylabel
%%   xaxis_metric            = 'linear' or 'sine'
%%   maintitle               = maintitle
%%   smooth                  = +1 for smothn, +0.5 for smoothdata(movmean,5),0 for no smooth
%%   clf                     = +1 to clear figure, -1 to leave

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

cT1 = min([nanmin(T11(:)) nanmin(T21(:))]);
cT2 = max([nanmax(T11(:)) nanmax(T21(:))]);

cW1 = min([nanmin(W11(:)) nanmin(W21(:))]);
cW2 = max([nanmax(W11(:)) nanmax(W21(:))]);

xlimits = [min(x) max(x)];  
ylimits = [min(y) max(y)];

plotoptions0.cx           = [cT1 cT1];
plotoptions0.cmap         = jett;
plotoptions0.yLinearOrLog = +1;
plotoptions0.yReverseDir  = -1;
plotoptions0.xLimits      = xlimits;
plotoptions0.yLimits      = ylimits;
plotoptions0.str11        = '(A)'; 
plotoptions0.str21        = '(B)'; 
plotoptions0.str11        = ' '; 
plotoptions0.str21        = ' '; 
plotoptions0.xstr         = 'Xstr';
plotoptions0.ystr         = 'Ystr';
plotoptions0.barstr       = 'barstr';
plotoptions0.xaxis_metric = 'linear';
plotoptions0.maintitle    = 'Main Title';
plotoptions0.clf          = +1;

%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 6
  iFig = 1;
  plotoptionsT = plotoptions0;
  plotoptionsA = plotoptions0;  
elseif nargin == 7
  plotoptionsT = plotoptions0;
  plotoptionsW = plotoptions0;
elseif nargin == 8
  plotoptionsW = plotoptions0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(plotoptionsT,'xaxis_metric')
  plotoptionsT.xaxis_metric = plotoptions0.xaxis_metric;
else
  if strfind(plotoptionsT.xaxis_metric,'sine')
    xlimits = [-1 +1];  
    plotoptions0.xLimits = xlimits;
    plotoptionsT.xLimits  = xlimits;
  end
end
if ~isfield(plotoptionsT,'smooth')
  plotoptionsT.smooth = plotoptions0.smooth;
end
if ~isfield(plotoptionsT,'clf')
  plotoptionsT.clf = plotoptions0.clf;
end
if ~isfield(plotoptionsT,'cx')
  plotoptionsT.cx = plotoptions0.cx;
end
if ~isfield(plotoptionsT,'cmap')
  plotoptionsT.cmap = plotoptions0.cmap;
end
if ~isfield(plotoptionsT,'yLinearOrLog')
  plotoptionsT.yLinearOrLog = plotoptions0.yLinearOrLog;
end
if ~isfield(plotoptionsT,'yReverseDir')
  plotoptionsT.yReverseDir = plotoptions0.yReverseDir;
end
if ~isfield(plotoptionsT,'xLimits')
  plotoptionsT.xLimits = plotoptions0.xLimits;
end
if ~isfield(plotoptionsT,'yLimits')
  plotoptionsT.yLimits = plotoptions0.yLimits;
end
if ~isfield(plotoptionsT,'str11')
  plotoptionsT.str11 = plotoptions0.str11;
end
if ~isfield(plotoptionsT,'str21')
  plotoptionsT.str21 = plotoptions0.str21;
end
if ~isfield(plotoptionsT,'xstr')
  plotoptionsT.xstr = plotoptions0.xstr;
end
if ~isfield(plotoptionsT,'ystr')
  plotoptionsT.ystr = plotoptions0.ystr;
end
if ~isfield(plotoptionsT,'barstr')
  plotoptionsT.barstr = plotoptions0.barstr;
end
if ~isfield(plotoptionsT,'maintitle')
  plotoptionsT.maintitle = plotoptions0.maintitle;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(plotoptionsW,'xaxis_metric')
  plotoptionsW.xaxis_metric = plotoptions0.xaxis_metric;
else
  if strfind(plotoptionsz.xaxis_metric,'sine')
    xlimits = [-1 +1];  
    plotoptions0.xLimits = xlimits;
    plotoptionsW.xLimits  = xlimits;
  end
end
if ~isfield(plotoptionsW,'smooth')
  plotoptionsz.smooth = plotoptions0.smooth;
end
if ~isfield(plotoptionsW,'clf')
  plotoptionsz.clf = plotoptions0.clf;
end
if ~isfield(plotoptionsW,'cx')
  plotoptionsz.cx = plotoptions0.cx;
end
if ~isfield(plotoptionsW,'cmap')
  plotoptionsz.cmap = plotoptions0.cmap;
end
if ~isfield(plotoptionsW,'yLinearOrLog')
  plotoptionsz.yLinearOrLog = plotoptions0.yLinearOrLog;
end
if ~isfield(plotoptionsW,'yReverseDir')
  plotoptionsz.yReverseDir = plotoptions0.yReverseDir;
end
if ~isfield(plotoptionsW,'xLimits')
  plotoptionsz.xLimits = plotoptions0.xLimits;
end
if ~isfield(plotoptionsW,'yLimits')
  plotoptionsz.yLimits = plotoptions0.yLimits;
end
if ~isfield(plotoptionsW,'str11')
  plotoptionsz.str11 = plotoptions0.str11;
end
if ~isfield(plotoptionsW,'str21')
  plotoptionsz.str21 = plotoptions0.str21;
end
if ~isfield(plotoptionsW,'xstr')
  plotoptionsz.xstr = plotoptions0.xstr;
end
if ~isfield(plotoptionsW,'ystr')
  plotoptionsz.ystr = plotoptions0.ystr;
end
if ~isfield(plotoptionsW,'barstr')
  plotoptionsz.barstr = plotoptions0.barstr;
end
if ~isfield(plotoptionsW,'maintitle')
  plotoptionsz.maintitle = plotoptions0.maintitle;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iFig); scr_siz = get(gcf); a0 = scr_siz.Position;
figure(iFig);
set(gcf,'resize','off')
set(gcf, 'Position',  [100, 100,  560*2, 420]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlimitsT = plotoptionsT.xLimits;
ylimitsT = plotoptionsT.yLimits;

xlimitsW = plotoptionsW.xLimits;
ylimitsW = plotoptionsW.yLimits;

%%want a 2x1 tiled layout
ta = tiledlayout(2,2,'TileSpacing','None', 'Padding','None');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

%if plotoptionsT.yLinearOrLog == +1
%  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
%else
%  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
%  %ta.OuterPosition = [0.0375 0.0375 0.875 0.875];
%  %ta.OuterPosition = [0.0375 0.0375 0.825 0.825];
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iFig);
% oo = find(isfinite(z11));

if strfind(plotoptionsT.xaxis_metric,'linear')
  disp('  xaxis will be linear in spacing (rlat)')
  iMetric = +1; %% default
elseif strfind(plotoptionsT.xaxis_metric,'sine')
  disp('  xaxis will be spaced by sin(rlat) so narrowest at poles')
  iMetric = +2;
else
  error('  xaxis needs sine or linear spacing (xaxis_metric)')
end

if plotoptionsT.smooth >= 1
  disp('using smoothn(z,1)')
elseif plotoptionsT.smooth == 0.5
  disp('using smooth(z,''movmean'',5)')
elseif plotoptionsT.smooth <= 0
  disp('no smoothing')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotoptionsT.smooth >= 1
  tah11 = smoothn(T11,1);
  tah21 = smoothn(T21,1);
elseif plotoptionsT.smooth == 0.5
  tah11 = smoothdata(T11,"movmean",5);
  tah21 = smoothdata(T21,"movmean",5);
elseif plotoptionsT.smooth <= 0
  tah11 = T11;
  tah21 = T21;
end

if plotoptionsW.smooth >= 1
  wah11 = smoothn(W11,1);
  wah21 = smoothn(W21,1);
elseif plotoptionsW.smooth == 0.5
  wah11 = smoothdata(W11,"movmean",5);
  wah21 = smoothdata(W21,"movmean",5);
elseif plotoptionsW.smooth <= 0
  wah11 = W11;
  wah21 = W21;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = plotoptionsT.cmap;
cxT   = plotoptionsT.cx;
cxW   = plotoptionsW.cx;

%% fontsize
fs = 14;
fs = 12;

if iMetric == 1
  tafov(1) = nexttile;
  pcolor(x,y,tah11); shading interp; colormap(plotoptionsT.cmap);
  box on
  ax = gca;
  ax.LineWidth = 4;
  caxis(cxT);
  ylim(plotoptionsT.yLimits)
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  ylabel('Pressure [mb]'); set(gca,'fontsize',fs);
  
  tafov(2) = nexttile;
  pcolor(x,y,wah11); shading interp; colormap(plotoptionsW.cmap);
  box on
  ax = gca;
  ax.LineWidth = 4;
  caxis(cxW);
  ylim(plotoptionsW.yLimits)
  set(gca,'ydir','reverse')
  ylabel('Pressure [mb]'); set(gca,'fontsize',fs);
  
  tafov(3) = nexttile;
  pcolor(x,y,tah21); shading interp; colormap(plotoptionsT.cmap);
  box on
  ax = gca;
  ax.LineWidth = 4;
  caxis(cxT);
  ylim(plotoptionsT.yLimits)
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')  
  colorbar('horizontal'); caxis(cxT);
  text(-110,15000,'K yr^{-1}','fontsize',fs);  
  ylabel('Pressure [mb]'); xlabel('Latitude [deg]'); set(gca,'fontsize',fs);
  
  tafov(4) = nexttile;
  pcolor(x,y,wah21); shading interp; colormap(plotoptionsW.cmap);
  box on
  ax = gca;
  ax.LineWidth = 4;
  caxis(cxW);
  ylim(plotoptionsW.yLimits)
  set(gca,'ydir','reverse')
  colorbar('horizontal'); caxis(cxW);
  text(-100,1500,'yr^{-1}','fontsize',fs);  
  ylabel('Pressure [mb]'); xlabel('Latitude [deg]'); set(gca,'fontsize',fs);
  
elseif iMetric == 2
  %% see MATLABCODE/PLOTTER/pcolor_sin
  xsin = sin(x*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

  tafov(1) = nexttile;
  pcolor(xsin,y,wah11); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  box on
  ax = gca;
  ax.LineWidth = 4;

  error('yuk')
  tafov(2) = nexttile;
  pcolor(xsin,y,wah21); shading interp; colormap(plotoptions.cmap);
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])
  box on
  ax = gca;
  ax.LineWidth = 4;
end

for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

strspace = 'none';
strspace = 'tight';
strspace = 'compact';

% Get rid of all extra space I can
ta.Padding = strspace;
ta.TileSpacing = strspace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iASLMAP = -1;
%common_2x1tiled_layout
