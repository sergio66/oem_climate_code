figure(41); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    %aslmap(41,rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,1),72,64)',1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]); 
    aslmapSergio(rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,2),72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(1),llsmap5); caxis([-1 +1]*0.15); colorbar
    text(-2.5,1.15,'(a)')
  tfov(2) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,3),72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(2),llsmap5); caxis([ 0 +1]/5*0.15); colorbar
    text(-2.5,1.15,'(b)')
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';

figure(42); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_RH(:,2),length(i100),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]*0.5); colormap(usa2); 
    text(-80,+200,'(a)');  ylabel('Pressure [mb]'); 
  tfov(2) = nexttile;
    pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_RH(:,3),length(i100),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([0 +1]/5*0.5); colormap(jet); 
    text(-80,+200,'(b)'); ylabel('Pressure [mb]'); xlabel('Latitude [deg]');
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = []; 

figure(43); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_WVfrac(:,2),length(i100),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]*0.015); colormap(usa2); 
    text(-80,+200,'(a)'); ylabel('Pressure [mb]');
  tfov(2) = nexttile;
    pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_WVfrac(:,3),length(i100),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([0 +1]/5*0.015); colormap(jet); 
    text(-80,+200,'(b)'); ylabel('Pressure [mb]'); xlabel('Latitude [deg]');
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = []; 

figure(44); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    pcolor(trend_rlat64,plays100(i10),reshape(frac_neg0pos_mean_std_T(:,2),length(i10),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([-1 +1]*0.15); colormap(usa2); 
    text(-80,+20,'(a)'); ylabel('Pressure [mb]');
  tfov(2) = nexttile;
    pcolor(trend_rlat64,plays100(i10),reshape(frac_neg0pos_mean_std_T(:,3),length(i10),64)); 
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([0 +1]/5*0.15); colormap(jet); 
    text(-80,+20,'(b)'); ylabel('Pressure [mb]'); xlabel('Latitude [deg]');
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(45); clf; aslmap(45,rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,1),72,64)',1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]);      title('dSKT/dt : sign agreement');
figure(46); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_RH(:,1),length(i100),64));     title('Fractional Agreement : dRH/dt');     
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]); colormap(usa2); 
figure(47); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_WVfrac(:,1),length(i100),64)); title('Fractional Agreement : dWVfrac/dt'); 
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]); colormap(usa2); 
figure(48); clf; pcolor(trend_rlat64,plays100(i10),reshape(frac_neg0pos_mean_std_T(:,1),length(i10),64));      title('Fractional Agreement : dT/dt');      
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([-1 +1]); colormap(usa2); 

fprintf('Sign agreement for dSKT/dt    = %8.6f out of %5i \n',nanmean(abs(frac_neg0pos_mean_std_stemp(:,1))),length(frac_neg0pos_mean_std_stemp));
fprintf('Sign agreement for dRH/dt     = %8.6f out of %5i \n',nanmean(abs(frac_neg0pos_mean_std_RH(:,1))),length(frac_neg0pos_mean_std_RH));
fprintf('Sign agreement for dWVfrac/dt = %8.6f out of %5i \n',nanmean(abs(frac_neg0pos_mean_std_WVfrac(:,1))),length(frac_neg0pos_mean_std_WVfrac));
fprintf('Sign agreement for dT/dt      = %8.6f out of %5i \n',nanmean(abs(frac_neg0pos_mean_std_T(:,1))),length(frac_neg0pos_mean_std_T));

%{
figure(45); clf; aslmap(45,rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,2),72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.15); title('dSKT/dt : mean over 6 results');
figure(46); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_RH(:,2),length(i100),64));     title('mean over 5 results :  dRH/dt');     
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]*0.5); colormap(usa2); 
figure(47); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_WVfrac(:,2),length(i100),64)); title('mean over 5 results :  dWVfrac/dt'); 
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]*0.015); colormap(usa2); 
figure(48); clf; pcolor(trend_rlat64,plays100(i10),reshape(frac_neg0pos_mean_std_T(:,2),length(i10),64));      title('mean over 5 results :  dT/dt');      
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([-1 +1]*0.15); colormap(usa2); 

figure(45); clf; aslmap(45,rlat65,rlon73,smoothn(reshape(frac_neg0pos_mean_std_stemp(:,3),72,64)',1), [-90 +90],[-180 +180]); colormap(jet); caxis([ 0 +1]/5*0.15); title('dSKT/dt : std over 6 results');
figure(46); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_RH(:,3),length(i100),64));     title('stddev over 5 results :  dRH/dt');     
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([0 +1]/5*0.5); colormap(jet); 
figure(47); clf; pcolor(trend_rlat64,plays100(i100),reshape(frac_neg0pos_mean_std_WVfrac(:,3),length(i100),64)); title('stddev over 5 results :  dWVfrac/dt'); 
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([0 +1]/5*0.015); colormap(jet); 
figure(48); clf; pcolor(trend_rlat64,plays100(i10),reshape(frac_neg0pos_mean_std_T(:,3),length(i10),64));      title('stddev over 5 results :  dT/dt');      
  colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([0 +1]/5*0.15); colormap(jet); 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%

