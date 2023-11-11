ta = tiledlayout(3,3,'TileSpacing','None', 'Padding','None');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(1,regions,:))'); colormap(tafov(1),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('RH (FRAC -0+)'); t.FontSize = 10; t.FontWeight = 'normal';

tafov(2) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(2,regions,:))'); colormap(tafov(2),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('WV (FRAC -0+)'); t.FontSize = 10; t.FontWeight = 'normal';

tafov(3) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(3,regions,:))'); colormap(tafov(3),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T(FRAC -0+)'); t.FontSize = 10; t.FontWeight = 'normal';

t = text(4.75,4,'200 mb','Rotation',90');

%%% 
tafov(4) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(2,regions,:))'); colormap(tafov(4),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('RH'); t.FontSize = 10; t.FontWeight = 'normal';

tafov(5) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(5,regions,:))'); colormap(tafov(5),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('WV'); t.FontSize = 10; t.FontWeight = 'normal';

tafov(6) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(8,regions,:))'); colormap(tafov(6),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('T'); t.FontSize = 10; t.FontWeight = 'normal';

t = text(4.75,4,'500 mb','Rotation',90');

%%% 
tafov(7) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(3,regions,:))'); colormap(tafov(7),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('RH'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');

tafov(8) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(6,regions,:))'); colormap(tafov(8),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('WV'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');

tafov(9) = nexttile;  imagesc(squeeze(allX_frac_neg0pos(7,regions,:))'); colormap(tafov(9),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames(1:4))
  set(gca,'xtick',[1:4],'xticklabel',{'GLOBAL','TROPICAL','MIDLAT','POLAR'});
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('T'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');

t = text(4.75,4,'800 mb','Rotation',90');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'tight';
ta.Padding = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 3 5 6 8 9]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end

% Remove all xtick labels except for 2nd row
for ii = [1 2 3 4 5 6]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end
