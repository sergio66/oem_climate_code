function [] = aslmap_2x3tiledlayout(z11,z12,z13,z21,z22,z23,iFig,plotoptions,x,y);

%% zXY = rowX, column Y   X=1-3,Y=1-2
%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/tiledlayout_pcolor_aslmap_rh_4tiles_eradonebyCLH.m

%{
EXAMPLE ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_gather_spectralrates_AIRSL3_NWP_XMIP6.m
plotoptions.str11 = 'AMIP6'; plotoptions.str12 = 'AMIP6';       plotoptions.str21 = 'CMIP6'; plotoptions.str22 = 'CMIP6';   plotoptions.str31 = 'UMBC'; plotoptions.str32 = 'UMBC';
  figure(7); clf; aslmap_6tiledlayout(amip6(1520,:),amip6_geo_rates,cmip6(1520,:),cmip6_geo_rates,x.rates(1520,:),umbcst,7,plotoptions);
%}

% addpath /asl/matlib/maps
% addpath /asl/matlib/plotutils
% addpath /home/sergio/MATLABCODE
% addpath /home/sergio/MATLABCODE/PLOTTER
% addpath /home/sergio/MATLABCODE/COLORMAP
% addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% plotoptions can have following fields
%%   cx                = [caxis1 caxis2]
%%   cmap              = colormap
%%   yLinearOrLog      = +1/-1
%%   yReverseDir       = +1/-1
%%   str11,str12,str3,str4 = give the subplots some flair
%%   xstr,ystr           = xlabel, ylabel
%%   maintitle           = maintitle
%%   shift180            = -1 for no, +1 for yes
%%   smooth              = +1 for smothn, +0.5 for smoothdata(movmean,5),0 for no smooth

%{
%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/tiledlayout_pcolor_aslmap_sst.m
load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end
%}

%[numtimesteps,~] = size(airsl3_64x72.all.mmw);
rlat = load('latB64.mat'); rlat65 = rlat.latB2; rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon73 = (1:73); rlon73 = -180 + (rlon73-1)*5;  rlon = (1:72); rlon = -177.5 + (rlon-1)*5;
[Y,X] = meshgrid(rlat,rlon); X = X(:); Y = Y(:);
%figure(6); clf; aslmap(6,rlat65,rlon73,smoothn(reshape(xtrend_st,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.15 +0.15]);

jet64 = jet(64);
jett = jet64; jett(1,:) = 1;

%if length(llsmap5) == 64
%  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
%  llsmap5 = llsmap5(2:64,:);  
%end

cx1 = min([nanmin(z11(:)) nanmin(z12(:)) nanmin(z13(:)) nanmin(z21(:)) nanmin(z22(:)) nanmin(z23(:))]);
cx2 = max([nanmax(z11(:)) nanmax(z12(:)) nanmax(z13(:)) nanmax(z21(:)) nanmax(z22(:)) nanmax(z23(:))]);
plotoptions0.cx           = [cx1 cx2];
plotoptions0.cmap   = jett;
plotoptions0.yLinearOrLog = +1;
plotoptions0.yReverseDir  = -1;
plotoptions0.str11 = '(A)'; plotoptions0.str12 = '(B)'; plotoptions0.str13 = '(C)';
plotoptions0.str21 = '(D)'; plotoptions0.str22 = '(E)'; plotoptions0.str23 = '(F)';
plotoptions0.xstr = 'Xstr';
plotoptions0.ystr = 'Ystr';
plotoptions0.barstr = 'barstr';
plotoptions0.maintitle = 'Main Title';
plotoptions0.shift180 = -1;
plotoptions0.smooth    = +0; %% no smoothing
plotoptions0.smooth    = +1; %% smooth Matlab, 1 pt window
plotoptions0.smooth    = +1; %% smoothn <default>

if nargin == 6
  iFig = 1;
  plotoptions = plotoptions0;
  x = rlon;
  y = rlat;
elseif nargin == 7
  plotoptions = plotoptions0;
  x = rlon;
  y = rlat;
elseif nargin == 8
  x = rlon;
  y = rlat;
elseif nargin == 9
  error('huh, you gave x but not y?????')
end
if nargin == 10
  [Y,X] = meshgrid(y,x); X = X(:); Y = Y(:);
  rlat65 = [y y(end)+mean(diff(y))]; 
  rlon73 = [x x(end)+mean(diff(x))];   
end
xlimits = [min(x) max(x)];  
ylimits = [min(y) max(y)];  

if ~isfield(plotoptions,'smooth')
  plotoptions.smooth = plotoptions0.smooth;
end
if ~isfield(plotoptions,'cx')
  plotoptions.cx = plotoptions0.cx;
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
if ~isfield(plotoptions,'str11')
  plotoptions.str11 = plotoptions0.str11;
end
if ~isfield(plotoptions,'str12')
  plotoptions.str12 = plotoptions0.str12;
end
if ~isfield(plotoptions,'str13')
  plotoptions.str13 = plotoptions0.str13;
end
if ~isfield(plotoptions,'str21')
  plotoptions.str21 = plotoptions0.str21;
end
if ~isfield(plotoptions,'str22')
  plotoptions.str22 = plotoptions0.str22;
end
if ~isfield(plotoptions,'str23')
  plotoptions.str23 = plotoptions0.str23;
end
if ~isfield(plotoptions,'barstr')
  plotoptions.barstr = plotoptions0.barstr;
end
if ~isfield(plotoptions,'xstr')
  plotoptions.xstr = plotoptions0.xstr;
end
if ~isfield(plotoptions,'ystr')
  plotoptions.ystr = plotoptions0.ystr;
end
if ~isfield(plotoptions,'maintitle')
  plotoptions.maintitle = plotoptions0.maintitle;
end
if ~isfield(plotoptions,'shift180')
  plotoptions.shift180 = plotoptions0.shift180;
end

if plotoptions.shift180 == +1
  shift180 = 180;
  xlimits = [min(wrapTo360(x)) max(wrapTo360(x))];
  X = wrapTo360(X);
  rlon73 = wrapTo360(rlon73);
else
  shift180 = 0;
end

fig = figure(iFig); clf; 
figure(iFig); scr_siz = get(gcf); a0 = scr_siz.Position;
figure(iFig); %set(gcf, 'Position',  [100, 100,  560*1.25, 420]);

%%want a 2x3 tiled layout
ta = tiledlayout(2,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cx   = plotoptions.cx;
cmap = plotoptions.cmap;

dxyzaspect = [1 1 1];
dxyzaspect = [1.5 0.75 1];
dxyzaspect = [1.25 0.75 1]; %% not bad
dxyzaspect = [1.8 0.8 1];   %% too tall and not very wide
dxyzaspect = [1.125 0.875 1]; %% not bad
dxyzaspect = [1.375 0.925 1]; %% not bad
dxyzaspect = [1.6 1.0 1]; %% not bad, from aslmap_2x5tiledlayout.m

%%%%%%%%%%%%%%%%%%%%%%%%%

if plotoptions.smooth >= 1
  disp('using smoothn(z,1)')
elseif plotoptions.smooth == 0.5
  disp('using smooth(z,''movmean'',5)')
elseif plotoptions.smooth <= 0
  disp('no smoothing')
end

oo = find(isfinite(z11));
tafov(1) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z11,72,64)',1), [-90 +90],[-180 +180]); 
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z11,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z11,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z11,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(2) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z12,72,64)',1), [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z12,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z12,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z12,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(3) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z13,72,64)',1), [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z13,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z13,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z13,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(4) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z21,72,64)',1), [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z21,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z21,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z21,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(5) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z22,72,64)',1), [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z22,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z22,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z22,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(6) = nexttile;
%aslmapSergio(rlat65,rlon73,smoothn(reshape(z23,72,64)',1), [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
if plotoptions.smooth >= 1
  aslmapSergio(rlat65,rlon73,smoothn(reshape(z23,72,64)',1), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth == 0.5
  aslmapSergio(rlat65,rlon73,smoothdata(reshape(z23,72,64)',"movmean",5), [-90 +90],[-180 +180]); 
elseif plotoptions.smooth <= 0
  aslmapSergio(rlat65,rlon73,reshape(z23,72,64)',[-90 +90],[-180 +180]);
end
colormap(plotoptions.cmap);
daspect(dxyzaspect)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iASLMAP = +1;
common_2x3tiled_layout

