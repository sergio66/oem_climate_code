figure(2);
disp('reading in the ERA figures to get ERA rates')
figname = 'WORKS_Nov19_2019_fit_straight_spectral_rate/Figs/era_raw_o3_decadal_trends.fig';
[latera,pressera,cera_o3] = get_fig_image_data(figname);
figname = 'WORKS_Nov19_2019_fit_straight_spectral_rate/Figs/era_raw_wv_decadal_trends.fig';
[latera,pressera,cera_wv] = get_fig_image_data(figname);
figname = 'WORKS_Nov19_2019_fit_straight_spectral_rate/Figs/era_raw_tz_decadal_trends.fig';
[latera,pressera,cera_tz] = get_fig_image_data(figname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jett = jet; jett(1,:) = 1;
llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); 

disp('ret to see T(z),WV(z),O3(z)'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('latbin stemp smoothed over 10 years')

  figure(4); pcolor(okdates,playsN,smoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('latbin WV smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar; 
  colormap(llsmap4.llsmap4);

  figure(5); pcolor(okdates,playsN,smoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('latbin Tz smoothed over 10 years')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  figure(6); pcolor(okdates,playsN,smootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('latbin O3 smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar
  colormap(llsmap4.llsmap4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iLoadERAanom = input('Enter (+1) to load ERA raw anom, (-1) no, default : ');
if length(iLoadERAanom) == 0
  iLoadERAanom = -1;
end

if iLoadERAanom > 0

%{
a = load('era_ptempanom.mat');
b = load('OutputAnomaly_OBS/21/anomtest_timestep300.mat');
aa = squeeze(a.era_ptempanom(21,:,300));
bb = b.oem.finalrates(b.jacobian.temp_i)';
sum(aa-bb)
%}

  load era_ptempanom.mat
  whos era_ptempanom smoothtz tz
  era_ptempanom = permute(era_ptempanom,[1 3 2]);

  figure(10); clf
  boo = tz - era_ptempanom;
  pcolor(squeeze(boo(20,:,:))); shading flat; colorbar; caxis([-10 +10]);

  mean_era_ptempanom = squeeze(nanmean(era_ptempanom(iaTropics,:,:),1));
  for ii = 1 : nlays
    smooth_era_ptempanom(:,ii) = smooth(mean_era_ptempanom(:,ii),2*5);
  end

  figure(7); pcolor(okdates,playsN,smooth_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('raw ERA latbin Tz smoothed over 10 years')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  figure(8); pcolor(okdates,playsN,smoothtz' - smooth_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('(latbin - raw ERA ) Tz smoothed over 10 years')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

%{
addpath /asl/matlib/plotutils
figure(7); hgsave rawERAanoms.fig
figure(5); hgsave fittedOBSanoms.fig
figure(8); hgsave diff_rawERA_fittedanoms.fig
%}

  figure(9); clf
  for ii = 1 : nlays
    B = Math_tsfit_lin_robust((1:365)*16,mean_era_ptempanom(:,jj),4); trend_mean_era_ptempanom(jj) = B(2);
    B = Math_tsfit_lin_robust((1:365)*16,smooth_era_ptempanom(:,jj),4); trend_mean_era_ptempanom(jj) = B(2);
  end
  for ii = 1 : 40
    for jj = 1 : nlays
      B = Math_tsfit_lin_robust((1:365)*16,squeeze(era_ptempanom(ii,:,jj)),4); trend_era_ptempanom(ii,jj) = B(2);
    end 
  end

  figure(9); pcolor(latbins,playsN,10*trend_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends for fixed Tanom')
  caxis([-0.05 +0.05]*10); colorbar
  colormap(jett)
  colormap(usa2)
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp

  semilogy(trend_mean_era_ptempanom,playsN,'b',nanmean(trendtz2(iaTropics,:)'),playsN,'k',trendtz2(iaTropics,:),playsN,'r'); grid
    hl = legend('<latbin raw ERA anom>','<latbin retrieval>','location','best');

  disp('ret'); pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('ret to see T(z),WV(z),O3(z) decadal trends'); pause
  %figure(7); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  %ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('latbin stemp smoothed over 10 years')

  figure(8); plot(10*trendwv2',playsN,'b',1.0*cera_wv(:,iReadLatbin),pressera,'r'); shading interp; colorbar
  ax = axis; ax(1) = -0.2; ax(2) = +0.2; ax(3) = 10; ax(4) = 1000;; axis(ax); grid; set(gca,'ydir','reverse'); title('decadal WV trends')
  caxis([-0.01 +0.01]*10); colorbar; 
  hl = legend('retr','ERA','location','best');

  figure(9); plot(10*trendtz2',playsN,'b',1.0*cera_tz(:,iReadLatbin),pressera,'r'); shading interp; colorbar
  ax = axis; ax(1) = -1; ax(2) = +1; ax(3) = 10; ax(4) = 1000;; axis(ax); grid; set(gca,'ydir','reverse'); title('decadal Tz trends')
  caxis([-0.05 +0.05]*10); colorbar
  hl = legend('retr','ERA','location','best');

  figure(10); plot(10*trendo32',playsN,'b',1.0*cera_o3(:,iReadLatbin),pressera,'r'); shading interp; colorbar
  ax = axis; ax(1) = -0.2; ax(2) = +0.2; ax(3) = 10; ax(4) = 1000;; axis(ax); grid; set(gca,'ydir','reverse'); title('decadal O3 trends')
  caxis([-0.02 +0.02]*10); colorbar
  hl = legend('retr','ERA','location','best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12); colormap(llsmap4.llsmap4);
figure(13); colormap(llsmap4.llsmap4);
i400 = abs(playsN-400); i400 = find(i400 == min(i400));
tzanom400    = tz(:,i400);
tzanomsig400 = tz_sig(:,i400);
figure(12); plot(okdates,tzanom400);          title('Tz anom at 400 mb');     ylim([-5 +5]);
figure(13); plot(okdates,real(tzanomsig400)); title('Tz unc anom at 400 mb'); ylim([0 0.1]);

