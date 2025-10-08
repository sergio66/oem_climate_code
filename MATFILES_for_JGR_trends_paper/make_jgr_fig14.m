load fig14.mat

load llsmap5

sizefont = 14;
ccc = [1 1 1]*0.5;  %% gray color dots

iSkip_1_2_3_4 = +1;
if iSkip_1_2_3_4 < 0
  figure(1); clf
    pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.pavg,fig14_ztest_unc.Tunc); set(gca,'ydir','reverse'); ylim([1 1000]); colormap jet; shading interp; % title('\sigma T'); ; plotaxis2;
    set(gca,'yscale','log'); ylim([10 1000])
    ylim([10 1000]); caxis([0 0.02]);  xlabel('Latitude [deg]'); ylabel('Pressure [mb]');
    % colorbar('vertical'); text(90,08,'K/yr'); set(gca,'fontsize',12)
    c = colorbar('horizontal'); text(-105,6000,'K/yr','fontsize',12); set(gca,'fontsize',sizefont)
    c.Ruler.Exponent = 0;
    c.Ruler.TickLabelFormat = '%0.3f';
  figure(3); clf
    pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.pavg,fig14_ztest_unc.WVunc); set(gca,'ydir','reverse'); ylim([100 1000]); colormap jet; shading interp; % title('\sigma WV'); plotaxis2;
    set(gca,'yscale','linear'); ylim([100 1000])
    ylim([100 1000]); caxis([0 0.003]); xlabel('Latitude [deg]'); ylabel('Pressure [mb]');
    % colorbar('vertical'); text(90,50,'1/yr'); set(gca,'fontsize',12)
    c = colorbar('horizontal'); text(-105,1350,'1/yr'); set(gca,'fontsize',sizefont)
    c.Ruler.Exponent = 0;
    c.Ruler.TickLabelFormat = '%0.3f';
  
  figure(2); clf
    pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.Ttrend_p,fig14_ztest_unc.Ttrend); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); caxis([-1 +1]*0.101); xlabel('Latitude [deg]'); ylabel('Pressure [mb]');
    colorbar('horizontal');
    ht = reshape(fig14_ztest_unc.Ttrend_ztest,64,101);
    for jj = 1 : 100
      for ii = 1 : 64
        if ht(ii,jj) == 1
          hold on; plot(fig14_ztest_unc.rlat(ii),fig14_ztest_unc.Ttrend_p(jj),'.','color',ccc,'Markersize',4,'linewidth',4);
        end
      end
    end
    hold off
    set(gca,'yscale','log'); ylim([10 1000])
  figure(4); clf
    pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.WVtrend_p,fig14_ztest_unc.WVtrend); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); caxis([-1 +1]*0.0101); xlabel('Latitude [deg]'); ylabel('Pressure [mb]');
    colorbar('horizontal');
    ht = reshape(fig14_ztest_unc.WVtrend_ztest,64,101);
    for jj = 1 : 100
      for ii = 1 : 64
        if ht(ii,jj) == 1
          hold on; plot(fig14_ztest_unc.rlat(ii),fig14_ztest_unc.WVtrend_p(jj),'.','color',ccc,'Markersize',4,'linewidth',4);
        end
      end
    end
    hold off
    ylim([100 1000])
  
  figure(1); set(gca,'fontsize',sizefont); text(-80,20, '(A)','fontsize',12); caxis([0 0.015]);
  figure(3); set(gca,'fontsize',sizefont); text(-80,200,'(C)','fontsize',12); caxis([0 2e-3]);
  figure(2); text(-105,6000,'K/yr','fontsize',12);; set(gca,'fontsize',sizefont); text(-80,20, '(B)','fontsize',12);
  figure(4); text(-105,1350,'1/yr');set(gca,'fontsize',sizefont); text(-80,200,'(D)','fontsize',12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2, and as tall as Fig 1 + Fig 3')
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2, and as tall as Fig 1 + Fig 3')
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2, and as tall as Fig 1 + Fig 3')

ta = tiledlayout(2,2);

tafov(1) = nexttile; 
  pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.pavg,fig14_ztest_unc.Tunc);
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); ylim([10 1000])
  colormap(tafov(1),jet); shading interp; % title('\sigma T'); ; plotaxis2;
  ylabel('Pressure [mb]'); % xlabel('Latitude [deg]'); 

  c = colorbar('horizontal'); text(-105,2000,'K/yr','fontsize',12); set(gca,'fontsize',sizefont)
  c.Ruler.Exponent = 0;
  c.Ruler.TickLabelFormat = '%0.3f';
  colorbarPosition = c.Position;
  newHeight = colorbarPosition(4) * 0.5; % Make it half its original height
  newBottom = colorbarPosition(2) + (colorbarPosition(4) - newHeight) / 2; % Center it vertically
  %c.Position = [colorbarPosition(1)+0.10, newBottom, colorbarPosition(3)-0.20, newHeight];
  c.Position = [colorbarPosition(1)+0.01, colorbarPosition(2), colorbarPosition(3)-0.02, colorbarPosition(4)];
  set(gca,'fontsize',sizefont); text(-80,20, '(A)','fontsize',12); caxis([0 0.015]);  

tafov(2) = nexttile; 
  pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.Ttrend_p,fig14_ztest_unc.Ttrend); 
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); ylim([10 1000])  
  colormap(tafov(2),llsmap5); shading interp; % title('dT/dt'); plotaxis2;
  caxis([-1 +1]*0.101); %xlabel('Latitude [deg]'); %ylabel('Pressure [mb]');
  set(gca,'fontsize',sizefont); text(-80,20, '(B)','fontsize',12);
  
  c = colorbar('horizontal'); text(-105,2000,'K/yr','fontsize',12);; set(gca,'fontsize',sizefont)
  colorbarPosition = c.Position; 
  newHeight = colorbarPosition(4) * 0.5; % Make it half its original height
  newBottom = colorbarPosition(2) + (colorbarPosition(4) - newHeight) / 2; % Center it vertically
  %c.Position = [colorbarPosition(1)+0.10, newBottom, colorbarPosition(3)-0.20, newHeight];
  c.Position = [colorbarPosition(1)+0.01, colorbarPosition(2), colorbarPosition(3)-0.02, colorbarPosition(4)];  

  ht = reshape(fig14_ztest_unc.Ttrend_ztest,64,101);
  for jj = 1 : 100
    for ii = 1 : 64
      if ht(ii,jj) == 1
        hold on; plot(fig14_ztest_unc.rlat(ii),fig14_ztest_unc.Ttrend_p(jj),'.','color',ccc,'Markersize',4,'linewidth',4);
      end
    end
  end
  hold off

tafov(3) = nexttile;   
  pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.pavg,fig14_ztest_unc.WVunc); 
  set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000])
  colormap(tafov(3),jet); shading interp; % title('\sigma WV'); plotaxis2;  
  xlabel('Latitude [deg]'); ylabel('Pressure [mb]');

  c = colorbar('horizontal'); text(-105,1200,'1/yr','fontsize',12); set(gca,'fontsize',sizefont)
  c.Ruler.Exponent = 0;
  c.Ruler.TickLabelFormat = '%0.3f';
  colorbarPosition = c.Position;
  newHeight = colorbarPosition(4) * 0.5; % Make it half its original height
  newBottom = colorbarPosition(2) + (colorbarPosition(4) - newHeight) / 2; % Center it vertically
  %c.Position = [colorbarPosition(1)+0.10, newBottom, colorbarPosition(3)-0.20, newHeight];
  c.Position = [colorbarPosition(1)+0.01, colorbarPosition(2)-0.0225, colorbarPosition(3)-0.02, colorbarPosition(4)];
  set(gca,'fontsize',sizefont); text(-80,200,'(C)','fontsize',12); caxis([0 2e-3]); 
  
tafov(4) = nexttile;   
  pcolor(fig14_ztest_unc.rlat,fig14_ztest_unc.WVtrend_p,fig14_ztest_unc.WVtrend); 
  set(gca,'ydir','reverse');  set(gca,'yscale','linear'); ylim([100 1000])
  colormap(tafov(4),llsmap5); shading interp; %title('dWVfrac/dt')
  caxis([-1 +1]*0.0101); xlabel('Latitude [deg]'); %ylabel('Pressure [mb]');
  set(gca,'fontsize',sizefont); text(-80,200,'(D)','fontsize',12);
  
  c = colorbar('horizontal'); text(-105,1200,'1/yr','fontsize',12);; set(gca,'fontsize',sizefont)
  colorbarPosition = c.Position; 
  newHeight = colorbarPosition(4) * 0.5; % Make it half its original height
  newBottom = colorbarPosition(2) + (colorbarPosition(4) - newHeight) / 2; % Center it vertically
  %c.Position = [colorbarPosition(1)+0.10, newBottom, colorbarPosition(3)-0.20, newHeight];
  c.Position = [colorbarPosition(1)+0.01, colorbarPosition(2)-0.0225, colorbarPosition(3)-0.01, colorbarPosition(4)];  

  ht = reshape(fig14_ztest_unc.WVtrend_ztest,64,101);
  for jj = 1 : 100
    for ii = 1 : 64
      if ht(ii,jj) == 1
        hold on; plot(fig14_ztest_unc.rlat(ii),fig14_ztest_unc.WVtrend_p(jj),'.','color',ccc,'Markersize',4,'linewidth',4);
      end
    end
  end
  hold off

  tafov(2).YTickLabel = ' ';
  tafov(4).YTickLabel = ' ';  

  ta.Padding = 'none';
  ta.TileSpacing = 'compact';
  ta.Padding = 'compact';
  ta.TileSpacing = 'tight';

  ta.TileSpacing = 'loose';
  ta.Padding = 'loose';

% sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/newfig14');
