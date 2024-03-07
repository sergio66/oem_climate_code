figure(1); clf; plot(1:len,thesummary.umbc_mean(:,1),'bo-',1:len,thesummary.merra2_mean(:,1),'g',1:len,thesummary.era5_mean(:,1),'r','linewidth',2); title('global dSKT/dt'); hl = legend('UMBC','MERRA2','ERA5','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

figure(2); clf; plot(1:len,thesummary.umbc_mean(:,2),'bo-',1:len,thesummary.merra2_mean(:,2),'g',1:len,thesummary.era5_mean(:,2),'r','linewidth',2); title('global dmmw/dt'); hl = legend('UMBC','MERRA2','ERA5','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

figure(3); clf; plot(1:len,thesummary.umbc_mean(:,3),'bo-',1:len,thesummary.merra2_mean(:,3),'g',1:len,thesummary.era5_mean(:,3),'r','linewidth',2); title('global dmmw percent/dt'); hl = legend('UMBC','MERRA2','ERA5','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf; 
ta = tiledlayout(3,1);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tafov(1) = nexttile; plot(1:len,thesummary.umbc_mean(:,4),'bo-',1:len,thesummary.merra2_mean(:,4),'g',1:len,thesummary.era5_mean(:,4),'r','linewidth',2);
    xlim([1 length(comment)]); hl = legend('UMBC','MERRA2','ERA5','location','best','fontsize',8);
    t = title('T/ML/P dSKT/dt'); t.FontSize = 12; t.FontWeight = 'normal';
  tafov(2) = nexttile; plot(1:len,thesummary.umbc_mean(:,7),'bo-',1:len,thesummary.merra2_mean(:,7),'g',1:len,thesummary.era5_mean(:,7),'r','linewidth',2);    xlim([1 length(comment)]); 
  tafov(3) = nexttile; plot(1:len,thesummary.umbc_mean(:,10),'bo-',1:len,thesummary.merra2_mean(:,10),'g',1:len,thesummary.era5_mean(:,10),'r','linewidth',2); xlim([1 length(comment)]); 
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',10);
  ax = gca; ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2 3]
   tafov(ii).FontSize = 8;
  end
if exist('iBestChoice'); arrayfun(@(ax)xline(ax,[iBestChoice iBestChoice], 'k','LineWidth',2), tafov); end

figure(5); clf; 
ta = tiledlayout(3,1);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tafov(1) = nexttile; plot(1:len,thesummary.umbc_mean(:,5),'bo-',1:len,thesummary.merra2_mean(:,5),'g',1:len,thesummary.era5_mean(:,5),'r','linewidth',2);    
    xlim([1 length(comment)]); hl = legend('UMBC','MERRA2','ERA5','location','best','fontsize',8);
    t = title('T/ML/P absolute dmmw/dt'); t.FontSize = 12; t.FontWeight = 'normal';
  tafov(2) = nexttile; plot(1:len,thesummary.umbc_mean(:,8),'bo-',1:len,thesummary.merra2_mean(:,8),'g',1:len,thesummary.era5_mean(:,8),'r','linewidth',2);    xlim([1 length(comment)]); 
  tafov(3) = nexttile; plot(1:len,thesummary.umbc_mean(:,11),'bo-',1:len,thesummary.merra2_mean(:,11),'g',1:len,thesummary.era5_mean(:,11),'r','linewidth',2); xlim([1 length(comment)]); 
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',10);
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',10);
  ax = gca; ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2 3]
   tafov(ii).FontSize = 8;
  end
if exist('iBestChoice'); arrayfun(@(ax)xline(ax,[iBestChoice iBestChoice], 'k','LineWidth',2), tafov); end

figure(6); clf; 
ta = tiledlayout(3,1);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tafov(1) = nexttile; plot(1:len,thesummary.umbc_mean(:,6),'bo-',1:len,thesummary.merra2_mean(:,6),'g',1:len,thesummary.era5_mean(:,6),'r','linewidth',2); 
    xlim([1 length(comment)]); hl = legend('UMBC','MERRA2','ERA5','location','best','fontsize',8);
    t = title('T/ML/P percent dmmw/dt'); t.FontSize = 12; t.FontWeight = 'normal';
  tafov(2) = nexttile; plot(1:len,thesummary.umbc_mean(:,9),'bo-',1:len,thesummary.merra2_mean(:,9),'g',1:len,thesummary.era5_mean(:,9),'r','linewidth',2);    xlim([1 length(comment)]); 
  tafov(3) = nexttile; plot(1:len,thesummary.umbc_mean(:,12),'bo-',1:len,thesummary.merra2_mean(:,12),'g',1:len,thesummary.era5_mean(:,12),'r','linewidth',2); xlim([1 length(comment)]); 
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',10);
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',10);
  ax = gca; ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2 3]
   xlim([1 length(comment)]); 
   tafov(ii).FontSize = 8;
  end
if exist('iBestChoice'); arrayfun(@(ax)xline(ax,[iBestChoice iBestChoice], 'k','LineWidth',2), tafov); end
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7); clf; 
ta = tiledlayout(2,1);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tafov(1) = nexttile; plot(1:len,thesummary.umbc_mean(:,13),'bo-',1:len,thesummary.merra2_mean(:,13),'g',1:len,thesummary.era5_mean(:,13),'r','linewidth',2);    title('O/L dSKT/dt'); hl = legend('UMBC','MERRA2','ERA5','location','best')';
    xlim([1 length(comment)]); 
  tafov(2) = nexttile; plot(1:len,thesummary.umbc_mean(:,16),'bo-',1:len,thesummary.merra2_mean(:,16),'g',1:len,thesummary.era5_mean(:,16),'r','linewidth',2);
    xlim([1 length(comment)]); 
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2]
   tafov(ii).FontSize = 8;
   xlim([1 length(comment)]); 
  end
if exist('iBestChoice'); arrayfun(@(ax)xline(ax,[iBestChoice iBestChoice], 'k','LineWidth',2), tafov); end

figure(8); clf; 
ta = tiledlayout(2,1);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tafov(1) = nexttile; plot(1:len,thesummary.umbc_mean(:,14),'bo-',1:len,thesummary.merra2_mean(:,14),'g',1:len,thesummary.era5_mean(:,14),'r','linewidth',2);    title('O/L absolute dmmw/dt'); hl = legend('UMBC','MERRA2','ERA5','location','best')';
    xlim([1 length(comment)]); 
  tafov(2) = nexttile; plot(1:len,thesummary.umbc_mean(:,17),'bo-',1:len,thesummary.merra2_mean(:,17),'g',1:len,thesummary.era5_mean(:,17),'r','linewidth',2);
    xlim([1 length(comment)]); 
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2]
   tafov(ii).FontSize = 8;
   xlim([1 length(comment)]); 
  end
if exist('iBestChoice'); arrayfun(@(ax)xline(ax,[iBestChoice iBestChoice], 'k','LineWidth',2), tafov); end
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9); clf; plot(1:len,thesummary.umbc_era5_corr(:,[1 4 7]),'--',1:len,thesummary.umbc_era5_corr(:,[2 5 8]),':',1:len,thesummary.umbc_era5_corr(:,[3 6 9 ]),'linewidth',2);
  title('Global Correlations between UMBC and ERA5'); hl = legend('200 mb RH','500 mb RH','800 mb RH','200 mb WV','500 mb WV','800 mb WV','200 mb T','500 mb T','800 mb T','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

figure(10); clf; plot(1:len,thesummary.umbc_era5_mean(:,[1 4 7]),'--',1:len,thesummary.umbc_era5_mean(:,[2 5 8]),':',1:len,thesummary.umbc_era5_mean(:,[3 6 9 ]),'linewidth',2);
  title('Global MEAN BIAS between UMBC and ERA5'); hl = legend('200 mb RH','500 mb RH','800 mb RH','200 mb WV','500 mb WV','800 mb WV','200 mb T','500 mb T','800 mb T','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

figure(11); clf; plot(1:len,thesummary.umbc_era5_std(:,[1 4 7]),'--',1:len,thesummary.umbc_era5_std(:,[2 5 8]),':',1:len,thesummary.umbc_era5_std(:,[3 6 9 ]),'linewidth',2);
  title('Global STD BIAS between UMBC and ERA5'); hl = legend('200 mb RH','500 mb RH','800 mb RH','200 mb WV','500 mb WV','800 mb WV','200 mb T','500 mb T','800 mb T','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end

figure(12); clf; plot(1:len,thesummary.umbc_era5_frac(:,[1 4 7]),'--',1:len,thesummary.umbc_era5_frac(:,[2 5 8]),':',1:len,thesummary.umbc_era5_frac(:,[3 6 9 ]),'linewidth',2);
  title('Global FRAC AGREEMENT between UMBC and ERA5'); hl = legend('200 mb RH','500 mb RH','800 mb RH','200 mb WV','500 mb WV','800 mb WV','200 mb T','500 mb T','800 mb T','location','best')';
  set(gca,'xtick',[1:len],'xticklabel',xstrstr,'fontsize',8);
    if exist('iBestChoice'); ax = axis; line([iBestChoice iBestChoice],[ax(3) ax(4)],'color','k','linewidth',2); end