%addpath Tiledplots   % LLS path to Sergio's codes
addpath /asl/matlib/maps  % LLS path to this
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

if ~exist('plotoptions','var')

  figure(1); close

  hgload /home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN/skt_Day_Night_avg_6panel.fig       what was used for Ryan Feb 2024 version
  % Get Sergio's plotoptions
  load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/skt_trends_strow.mat','plotoptions');
end

% Set aspect ratio (by hand first time)
set(gcf,'Position',[1103 512 1007 382])

hf = get(gcf,'Children');
htc = get(hf(2),'Children')   % I think TiledChartLayout will always be hf(2)??

% Tile Titles and ind holds members of htc that are maps (axes), not colorbars
atype = 'axes';
ind = [];
for i = 1:length(htc);
   ttype = get(htc(i),'type');
   if all(ismember(ttype,atype)) ind = [ind i];
   end
   ttitle(i).text = htc(i).Title.String;
end

% Get the data
k = 1;
for i = ind
   md = findobj(htc(i),'type','Surface')
   x(i,:,:) = md.XData;
   y(i,:,:) = md.YData;
   c(i,:,:) = md.CData;
end

% Just Sergio's rlat65 rlat73
%load figlatlon
do_XX_YY_from_X_Y

% tricks to find axes title and change one of them!
for i=ind
   if contains('CHIRP\_A',ttitle(i).text)
      i
      ttitle(i).text = 'AIRS\_RT';
   end
end

% Matlab mapping toolbox adds NaN's to end of mapping arrays that it returns
cc = squeeze(c(:,1:end-1,1:end-1));

% Finally
ta = tiledlayout(2,3);

for i=flip(ind)
  tafov(i) = nexttile;
  aslmapSergio(rlat65,rlon73,squeeze(cc(i,:,:)), [-90 +90],[-180 +180]);colormap(plotoptions.cmap);
  title(ttitle(i).text);caxis([-0.15 0.1501]);%caxis([-0.12 0.12]);
end

cb = colorbar;
cb.Layout.Tile = 'south'
xl = xlabel(tafov(ind(2)),'T_{surface} (K/year)')

ta.Padding = 'compact';
ta.TileSpacing = 'compact';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
BAD IDEA aslprint_asis('fixed_6panel_SKTtrend.pdf');
BAD IDEA hf = gcf;
BAD IDEA exportgraphics(hf,'fixed_6panel_SKTtrend.png')
%}

%{
The first one isn’t that good in .pdf since it is an embedded image,
and for a paper Backgroundcolor none isn’t needed since white is the
default for paper.  The second one is probably best for a paper, and
again Backgroundcolor none is not needed for a paper, but it IS needed
for a beamer presentation since we use a colored background.  It is
high quality pdf text too.  So, this is best for paper.

The third one is probably best for presentations since most people
don’t zoom them and the .png is quite readable.  It appears to have a
transparent background anyway OR it has the identical background color
that beamer uses (which I may be setting somewhere in the beamer
code).

So, #2 for papers, #3 for presentations.  A presentation using all
three is attached in next message.

hf = gcf;
exportgraphics(hf,'LLS_Figs/fig8x_defaultpdf.pdf', 'Backgroundcolor','none');
exportgraphics(hf,'LLS_Figs/fig8x_defaultpdf2.pdf','Backgroundcolor','none','ContentType','vector');
exportgraphics(hf,'LLS_Figs/fig8x_defaultpng.png'); 
%}

