sizefont = 10;
sizefont = 16;
sizefont = 12;
sizefont = 14;
sizefont = 13;
fprintf(1,'common_2x1tiled_layout.m : sizefont = %2i \n',sizefont)
tafov(1).FontSize = sizefont;
tafov(2).FontSize = sizefont;

if plotoptions.yReverseDir == 1
  set(tafov,'ydir','reverse');   
end
if plotoptions.yLinearOrLog == -1
  set(tafov,'yscale','log');  
end
set(tafov,'colormap',cmap,'CLim',cx);
%% text(-100,1500,'yr^{-1}','fontsize',sizefont)     %% figure 65 of ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_trendspaper_figure15_uncertainty_ztest.m
%% text(-110,10000,'K yr^{-1}','fontsize',sizefont)  %% figure 68 of ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_trendspaper_figure15_uncertainty_ztest.m

if iASLMAP == -1 | iASLMAP == +2
  set(tafov,'Xlim',xlimits);
  set(tafov,'Ylim',ylimits);
  %set(tafov,'shading','interp');
end

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

ta.Padding = 'compact';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
%for ii = [2 4]
%   tafov(ii).YTickLabel = '';
%   tafov(ii).YLabel.String = [];
%end
% Remove all xtick labels except for 3rd row
for ii = [1]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

xstr = plotoptions.xstr;
ystr = plotoptions.ystr;

if iASLMAP == +2
  % Put in xlabel and ylabel in the “middle”
  tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = sizefont;
  tafov(2).YLabel.String = ystr; tafov(2).YLabel.FontSize = sizefont;
else
  % Put in xlabel and ylabel in the “middle”
  tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = sizefont;
  tafov(2).YLabel.String = ystr; tafov(2).YLabel.FontSize = sizefont;
  tafov(2).XLabel.String = xstr; tafov(2).XLabel.FontSize = sizefont;
  %tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = sizefont;
  %tafov(3).YLabel.String = ystr; tafov(3).YLabel.FontSize = sizefont;
  %tafov(5).YLabel.String = ystr; tafov(5).YLabel.FontSize = sizefont;
  %tafov(5).XLabel.String = xstr; tafov(5).XLabel.FontSize = sizefont;
  %tafov(6).XLabel.String = xstr; tafov(6).XLabel.FontSize = sizefont;
end  

%% put titles
title(tafov(1), plotoptions.str11, 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
title(tafov(2), plotoptions.str21, 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = axes(fig);
han = gca;
han.Visible = 'off';
% X label
%%xlabel(xstr,'Fontsize',16)
%%han.XLabel.Visible = 'on';

% Left label
%yyaxis(ax, 'left');
%% ylabel(ystr,'Fontsize',16)
%% han.YLabel.Visible = 'on';
% Right label
%yyaxis(ax, 'right');
%ylabel('Y2')
%han.YLabel.Visible = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main title
%title(ta,plotoptions.maintitle);

%keyboard_nowindow
%if plotoptions.yLinearOrLog == -1
%  title(ta,plotoptions.maintitle);
%end

%% assign colorbar to one tile
cbh = colorbar(tafov(2));
cbh.Layout.Tile = 'south';
%cbh.Title.String = plotoptions.maintitle;
%cbh.Title.String = plotoptions.barstr;
%set(cbh.Title,{'String','Rotation','Position','FontSize'},{cbh.Title.String,0,[200 -35],15});


%%%cbh.Title.String = xstr;
%%%set(cbh.Title,{'String','Rotation','Position','FontSize'},{xstr,0,[200 +15],15})

%{
boo = get(cbh,'position')
%where the position arguments are [xposition yposition width height].
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
