function [] = aslmap_3x4tiledlayout(z11,z12,z13,z14,z21,z22,z23,z24,z31,z32,z33,z34,iFig,plotoptions,x,y);

%% zXY = rowX, column Y   X=1-2,Y=1-4
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

cx1 = min([nanmin(z11(:)) nanmin(z12(:)) nanmin(z13(:)) nanmin(z14(:)) nanmin(z21(:)) nanmin(z22(:)) nanmin(z23(:)) nanmin(z24(:)) nanmin(z31(:)) nanmin(z32(:)) nanmin(z33(:)) nanmin(z34(:)) ]);
cx2 = max([nanmax(z11(:)) nanmax(z12(:)) nanmax(z13(:)) nanmax(z14(:)) nanmax(z21(:)) nanmax(z22(:)) nanmax(z23(:)) nanmax(z24(:)) nanmax(z31(:)) nanmax(z32(:)) nanmax(z33(:)) nanmax(z34(:)) ]);

plotoptions0.cx           = [cx1 cx2];
plotoptions0.cmap   = jett;
plotoptions0.yLinearOrLog = +1;
plotoptions0.yReverseDir  = -1;
plotoptions0.str11 = '(A)'; plotoptions0.str12 = '(B)'; plotoptions0.str13 = '(C)'; plotoptions0.str14 = '(D)';
plotoptions0.str21 = '(E)'; plotoptions0.str22 = '(F)'; plotoptions0.str23 = '(G)'; plotoptions0.str24 = '(H)';
plotoptions0.str31 = '(I)'; plotoptions0.str32 = '(J)'; plotoptions0.str33 = '(K)'; plotoptions0.str34 = '(L)';
plotoptions0.xstr = 'Xstr';
plotoptions0.ystr = 'Ystr';
plotoptions0.barstr = 'barstr';
plotoptions0.maintitle = 'Main Title';
plotoptions0.shift180 = -1;
plotoptions0.smooth    = +0; %% no smoothing
plotoptions0.smooth    = +1; %% smooth Matlab, 1 pt window
plotoptions0.smooth    = +1; %% smoothn <default>

if nargin == 12
  iFig = 1;
  plotoptions = plotoptions0;
  x = rlon;
  y = rlat;
elseif nargin == 13
  plotoptions = plotoptions0;
  x = rlon;
  y = rlat;
elseif nargin == 14
  x = rlon;
  y = rlat;
elseif nargin == 15
  error('huh, you gave x but not y?????')
end
if nargin == 16
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
if ~isfield(plotoptions,'str14')
  plotoptions.str14 = plotoptions0.str14;
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
if ~isfield(plotoptions,'str24')
  plotoptions.str24 = plotoptions0.str24;
end
if ~isfield(plotoptions,'str31')
  plotoptions.str31 = plotoptions0.str31;
end
if ~isfield(plotoptions,'str32')
  plotoptions.str32 = plotoptions0.str32;
end
if ~isfield(plotoptions,'str33')
  plotoptions.str33 = plotoptions0.str33;
end
if ~isfield(plotoptions,'str34')
  plotoptions.str34 = plotoptions0.str34;
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

if plotoptions.smooth >= 1
  disp('using smoothn(z,1)')
elseif plotoptions.smooth == 0.5
  disp('using smooth(z,''movmean'',5)')
elseif plotoptions.smooth <= 0
  disp('no smoothing')
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

%%want a 2x4 tiled layout
ta = tiledlayout(3,4);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cx   = plotoptions.cx;
cmap = plotoptions.cmap;

dxyzaspect = [1 1 1];
dxyzaspect = [1.5 0.75 1];
dxyzaspect = [1.25 0.75 1]; %% not bad
dxyzaspect = [1.8 0.8 1];   %% too tall and not very wide
dxyzaspect = [1.125 0.875 1]; %% not bad
dxyzaspect = [1.375 0.925 1]; %% not bad

if plotoptions.smooth >= 1
  wah11 = smoothn(reshape(z11,72,64)',1);
  wah12 = smoothn(reshape(z12,72,64)',1);
  wah13 = smoothn(reshape(z13,72,64)',1);
  wah14 = smoothn(reshape(z14,72,64)',1);
  wah21 = smoothn(reshape(z21,72,64)',1);
  wah22 = smoothn(reshape(z22,72,64)',1);
  wah23 = smoothn(reshape(z23,72,64)',1);
  wah24 = smoothn(reshape(z24,72,64)',1);
  wah31 = smoothn(reshape(z31,72,64)',1);
  wah32 = smoothn(reshape(z32,72,64)',1);
  wah33 = smoothn(reshape(z33,72,64)',1);
  wah34 = smoothn(reshape(z34,72,64)',1);
elseif plotoptions.smooth == 0.5
  wah11 = smoothdata(reshape(z11,72,64)',"movmean",5);
  wah12 = smoothdata(reshape(z12,72,64)',"movmean",5);
  wah13 = smoothdata(reshape(z13,72,64)',"movmean",5);
  wah14 = smoothdata(reshape(z14,72,64)',"movmean",5);
  wah21 = smoothdata(reshape(z21,72,64)',"movmean",5);
  wah22 = smoothdata(reshape(z22,72,64)',"movmean",5);
  wah23 = smoothdata(reshape(z23,72,64)',"movmean",5);
  wah24 = smoothdata(reshape(z24,72,64)',"movmean",5);
  wah31 = smoothdata(reshape(z31,72,64)',"movmean",5);
  wah32 = smoothdata(reshape(z32,72,64)',"movmean",5);
  wah33 = smoothdata(reshape(z33,72,64)',"movmean",5);
  wah34 = smoothdata(reshape(z34,72,64)',"movmean",5);
elseif plotoptions.smooth <= 0
  wah11 = reshape(z11,72,64)';
  wah12 = reshape(z12,72,64)';
  wah13 = reshape(z13,72,64)';
  wah14 = reshape(z14,72,64)';
  wah21 = reshape(z21,72,64)';
  wah22 = reshape(z22,72,64)';
  wah23 = reshape(z23,72,64)';
  wah24 = reshape(z24,72,64)';
  wah31 = reshape(z31,72,64)';
  wah32 = reshape(z32,72,64)';
  wah33 = reshape(z33,72,64)';
  wah34 = reshape(z34,72,64)';
end

tafov(1) = nexttile;
oo = find(isfinite(z11));
figure(iFig);
aslmapSergio(rlat65,rlon73,wah11, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(2) = nexttile;
aslmapSergio(rlat65,rlon73,wah12, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(3) = nexttile;
aslmapSergio(rlat65,rlon73,wah13, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(4) = nexttile;
aslmapSergio(rlat65,rlon73,wah14, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(5) = nexttile;
aslmapSergio(rlat65,rlon73,wah21, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(6) = nexttile;
aslmapSergio(rlat65,rlon73,wah22, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(7) = nexttile;
aslmapSergio(rlat65,rlon73,wah23, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(8) = nexttile;
aslmapSergio(rlat65,rlon73,wah24, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(9) = nexttile;
aslmapSergio(rlat65,rlon73,wah31, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(10) = nexttile;
aslmapSergio(rlat65,rlon73,wah32, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(11) = nexttile;
aslmapSergio(rlat65,rlon73,wah33, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

tafov(12) = nexttile;
aslmapSergio(rlat65,rlon73,wah34, [-90 +90],[-180 +180]); colormap(plotoptions.cmap);
daspect(dxyzaspect)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iASLMAP = +1;
common_3x4tiled_layout

