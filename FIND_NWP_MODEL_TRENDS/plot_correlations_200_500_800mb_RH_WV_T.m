%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% allX_BLAH(iY,iX,:) will have size 9 x 5 x 4 
  %%     == [200/500/800 RH/WV/T] X ['ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'] x [GLOBAL TRP MIDLAT POLAR]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total   [RH,WVfrac,T][RH,WVfrac,T][RH,WVfrac,T]
  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% whos allXchi = 9 x 5 x 4                                       correlate ERA5 with [AIRS L3, CLIMCAPS, MERRA2, UMBC]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remember IX = iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total but I don't care about midlats+tropics = region4
%% so regions should be 1,2,3,5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ta = tiledlayout(3,3,'TileSpacing','None', 'Padding','None');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;  imagesc(squeeze(allXchi(1,regions,:))'); colormap(tafov(1),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('RH (CORRELATION)'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(1,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(1,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(2) = nexttile;  imagesc(squeeze(allXchi(2,regions,:))'); colormap(tafov(2),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('WV (CORRELATION)'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(2,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(2,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(3) = nexttile;  imagesc(squeeze(allXchi(3,regions,:))'); colormap(tafov(3),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  t = title('T(CORRELATION)'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(3,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(3,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

t = text(4.75,4,'200 mb','Rotation',90');

%%% 
tafov(4) = nexttile;  imagesc(squeeze(allXchi(4,regions,:))'); colormap(tafov(4),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('RH'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(4,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(4,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(5) = nexttile;  imagesc(squeeze(allXchi(5,regions,:))'); colormap(tafov(5),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('WV'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(5,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(5,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(6) = nexttile;  imagesc(squeeze(allXchi(6,regions,:))'); colormap(tafov(6),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('T'); t.FontSize = 10; t.FontWeight = 'normal';
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(6,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(6,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

t = text(4.75,4,'500 mb','Rotation',90');

%%% 
tafov(7) = nexttile;  imagesc(squeeze(allXchi(7,regions,:))'); colormap(tafov(7),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('RH'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(7,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(7,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(8) = nexttile;  imagesc(squeeze(allXchi(8,regions,:))'); colormap(tafov(8),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('WV'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(8,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(8,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

tafov(9) = nexttile;  imagesc(squeeze(allXchi(9,regions,:))'); colormap(tafov(9),jet); caxis([0 1]); 
  set(gca,'ytick',[1:4],'yticklabel',mnames);
  set(gca,'xtick',[1:4],'xticklabel',znames);
  ax = gca; ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;
  %t = title('T'); t.FontSize = 10; t.FontWeight = 'normal';
  colorbar('southoutside');
if iWriteCorrelNumbers > 0
  wah = squeeze(allXchi(9,regions,:))';
  for iii = 1 : 4
    for jjj = 1 : 4
      junk = allXchi(9,regions(iii),jjj); junkstr = num2str(junk,'%4.3f');
      text(iii-0.30,jjj-0.10,junkstr,'Fontsize',7);
    end
  end
end

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
