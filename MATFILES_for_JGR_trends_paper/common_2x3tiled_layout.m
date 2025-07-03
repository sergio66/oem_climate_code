tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;
tafov(5).FontSize = 10;
tafov(6).FontSize = 10;

if plotoptions.yReverseDir == 1
  set(tafov,'ydir','reverse');   
end
if plotoptions.yLinearOrLog == -1
  set(tafov,'yscale','log');  
end
set(tafov,'colormap',cmap,'CLim',cx);

if iASLMAP == -1
  set(tafov,'Xlim',xlimits);
  set(tafov,'Ylim',ylimits);
  %set(tafov,'shading','interp');
end

if iASLMAP == +1
  % Get rid of all extra space I can
  ta.Padding = 'none';
  ta.TileSpacing = 'none';

  ta.Padding = 'tight';
  ta.TileSpacing = 'tight';

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';

else
  ta.Padding = 'loose';
  ta.TileSpacing = 'loose';

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
end

% Remove all ytick labels except for 1st column
for ii = [2 3 5 6]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end
% Remove all xtick labels except for 3rd row
for ii = [1 2 3 5]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

xstr = plotoptions.xstr;
ystr = plotoptions.ystr;
% Put in xlabel and ylabel in the “middle”
%tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = 10;
%tafov(3).YLabel.String = ystr; tafov(3).YLabel.FontSize = 10;
%tafov(5).YLabel.String = ystr; tafov(5).YLabel.FontSize = 10;
%tafov(5).XLabel.String = xstr; tafov(5).XLabel.FontSize = 10;
%tafov(6).XLabel.String = xstr; tafov(6).XLabel.FontSize = 10;

%% put titles
if iASLMAP == +1
  ytitle = 1.01;
  title(tafov(1), plotoptions.str11, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(2), plotoptions.str12, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(3), plotoptions.str13, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(4), plotoptions.str21, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(5), plotoptions.str22, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(6), plotoptions.str23, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
elseif iASLMAP == -1
  ytitle = 1.01;
  title(tafov(1), plotoptions.str11, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(2), plotoptions.str12, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(3), plotoptions.str13, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(4), plotoptions.str21, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(5), plotoptions.str22, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
  title(tafov(6), plotoptions.str23, 'Units', 'normalized', 'Position', [0.5, ytitle, 0]);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = axes(fig);
han = gca;
han.Visible = 'off';
% X label
%%%% xlabel(xstr,'Fontsize',16)
han.XLabel.Visible = 'on';

% Left label
%yyaxis(ax, 'left');
%%%%% ylabel(ystr,'Fontsize',16)
han.YLabel.Visible = 'on';

% Right label
%yyaxis(ax, 'right');
%ylabel('Y2')
%han.YLabel.Visible = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main title
%title(ta,plotoptions.maintitle)

%% assign colorbar to one tile
cbh = colorbar(tafov(5));
cbh.Layout.Tile = 'south';
cbh.Title.String = plotoptions.maintitle;
cbh.Title.String = plotoptions.barstr;
set(cbh.Title,{'String','Rotation','Position','FontSize'},{cbh.Title.String,0,[200 -35],15});

%{
% https://www.mathworks.com/matlabcentral/answers/861340-colorbar-spacing-issue-using-tiledlayout
placeholder1 = tiledlayout('flow','Parent',ta);
placeholder1.Layout.Tile = 'south';
cbh = colorbar(tafov(5),'FontName','Arial Narrow','FontSize',8); 
%cbh.Label.String = plotoptions.maintitle;
cbh.Layout.Tile  = 'south';
set(cbh.Title,{'String','Rotation','Position','FontSize'},{plotoptions.maintitle,0,[50 +20],8})
%set(cbh.Title,{'String','Rotation','Position','FontSize'},{plotoptions.maintitle,0,[180 -1],8})
placeholder2 = tiledlayout('flow','Parent',ta);
placeholder2.Layout.Tile = 'south';
%placeholder3=tiledlayout('flow','Parent',ta);
%placeholder3.Layout.Tile='south';
%}

%{
boo = get(cbh,'position')
%where the position arguments are [xposition yposition width height].
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
