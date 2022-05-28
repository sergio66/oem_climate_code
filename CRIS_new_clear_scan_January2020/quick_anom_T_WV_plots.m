latbins = latbins(1:38);

jett = jet; jett(1,:) = 1;
llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); 

disp('ret to see T(z),WV(z),O3(z)'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 1 year')

  figure(4); pcolor(okdates,playsN,smoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical WV smoothed over 1 year')
  caxis([-0.2 +0.2]); colorbar; 
  colormap(llsmap4.llsmap4);
  hold on
    plot(okdates,500-300*smooth(meanstemp,2*5),'k',okdates,500*ones(size(okdates)),'k--',okdates,500-300*smooth(mean(save_dat_1231(:,iaTropics)'),10),'rx-','linewidth',2);
  hold off
  axis(ax);

  figure(5); pcolor(okdates,playsN,smoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical Tz smoothed over 1 year')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);
  hold on
    plot(okdates,500-300*smooth(meanstemp,2*5),'k',okdates,500*ones(size(okdates)),'k--',okdates,500-300*smooth(mean(save_dat_1231(:,iaTropics)'),10),'rx-','linewidth',2);
  hold off
  axis(ax);

  figure(6); pcolor(okdates,playsN,smootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 1; ax(4) = 100;; grid; set(gca,'ydir','reverse'); title('tropical O3 smoothed over 1 year')
  caxis([-0.2 +0.2]); colorbar
  colormap(llsmap4.llsmap4);
  hold on
    plot(okdates,50-30*smooth(meanstemp,2*5),'k',okdates,50*ones(size(okdates)),'k--',okdates,50-30*smooth(mean(save_dat_1231(:,iaTropics)'),10),'rx-','linewidth',2);
  hold off
  axis(ax);

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
  era_ptempanom = era_ptempanom(1:38,:,:);

  figure(10); clf
  boo = tz - era_ptempanom;
  boo = squeeze(boo(iaTropics,:,:));
  boo = squeeze(mean(boo,1));
  pcolor(okdates,1:20,boo'); shading flat; colorbar; caxis([-10 +10]);

  mean_era_ptempanom = squeeze(nanmean(era_ptempanom(iaTropics,:,:),1));
  for ii = 1 : nlays
    smooth_era_ptempanom(:,ii) = smooth(mean_era_ptempanom(:,ii),2*5);
  end

  figure(7); pcolor(okdates,playsN,smooth_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('raw ERA tropical Tz smoothed over 1 year')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  figure(8); pcolor(okdates,playsN,smoothtz' - smooth_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('(tropical - raw ERA ) Tz smoothed over 1 year')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  %{
  addpath /asl/matlib/plotutils
  figure(7); hgsave rawERAanoms.fig
  figure(5); hgsave fittedOBSanoms.fig
  figure(8); hgsave diff_rawERA_fittedanoms.fig
  %}

  figure(11); clf
  for ii = 1 : nlays
    B = Math_tsfit_lin_robust((1:157)*16,mean_era_ptempanom(:,jj),4); trend_mean_era_ptempanom(jj) = B(2);
    B = Math_tsfit_lin_robust((1:157)*16,smooth_era_ptempanom(:,jj),4); trend_mean_era_ptempanom(jj) = B(2);
  end
  for ii = 1 : 38
    for jj = 1 : nlays
      B = Math_tsfit_lin_robust((1:157)*16,squeeze(era_ptempanom(ii,:,jj)),4); trend_era_ptempanom(ii,jj) = B(2);
    end 
  end

  figure(11); pcolor(latbins,playsN,10*trend_era_ptempanom'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends for ERA')
  caxis([-0.05 +0.05]*10); colorbar
  caxis([-0.1  +0.1 ]*10); colorbar
  colormap(jett)
  colormap(usa2)
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp

  %semilogy(trend_mean_era_ptempanom,playsN,'b',nanmean(trendtz2(iaTropics,:)'),playsN,'k',trendtz2(iaTropics,:),playsN,'r');grid
  %  hl = legend('<tropical raw ERA anom>','<tropical retrieval>','location','best');

  disp('ret to go back to retrievals'); pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latbins = latbins(1:38);

disp('ret to see T(z),WV(z),O3(z) decadal trends'); pause
  figure(7); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 1 year')

  figure(8); pcolor(latbins,playsN,10*trendwv2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal WV trends')
  caxis([-0.01 +0.01]*10); colorbar; 

  figure(9); pcolor(latbins,playsN,10*trendtz2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends')
  caxis([-0.05 +0.05]*10); colorbar

  figure(10); pcolor(latbins,playsN,10*trendo32'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal O3 trends')
  caxis([-0.02 +0.02]*10); colorbar

for ii = 7 : 10
  figure(ii)
  colormap(jett)
  colormap(usa2)
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPrint = +1;
iPrint = input('print these (default = NO -1/+1) : ');
if length(iPrint) == 0
  iPrint = -1;
end
if iPrint > 0
  addpath /asl/matlib/plotutils/

  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; 
  %title('tropical stemp smoothed over 1 year')

  figure(4); pcolor(okdates,playsN,smoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); 
  %title('tropical WV smoothed over 1 year')
  caxis([-0.2 +0.2]); colorbar; 
  set(gca,'yscale','log'); axis([2012.35 2018.75 100 1000]); xlabel('time'); ylabel('P(mb)')
  colormap(llsmap4.llsmap4);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  

  figure(5); pcolor(okdates,playsN,smoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); 
  %title('tropical Tz smoothed over 1 year')
  caxis([-1 +1]); colorbar
  set(gca,'yscale','log'); axis([2012.35 2018.75 100 1000]); xlabel('time'); ylabel('P(mb)')
  colormap(llsmap4.llsmap4);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  

  figure(6); pcolor(okdates,playsN,smootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 1; ax(4) = 100;; grid; set(gca,'ydir','reverse'); 
  %title('tropical O3 smoothed over 1 year')
  caxis([-0.2 +0.2]); colorbar
  set(gca,'yscale','log'); axis([2012.35 2018.75 100 1000]); xlabel('time'); ylabel('P(mb)')
  colormap(llsmap4.llsmap4);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'1','10','100'})  

  figure(8); pcolor(latbins,playsN,10*trendwv2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal WV trends')
  caxis([-0.01 +0.01]*10); colorbar; 
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp

  figure(9); pcolor(latbins,playsN,10*trendtz2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends')
  caxis([-0.05 +0.05]*10); colorbar
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp

  figure(10); pcolor(latbins,playsN,10*trendo32'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal O3 trends')
  caxis([-0.02 +0.02]*10); colorbar
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp

  if iOBSorCAL == 0
    figure(5); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_ptemp_anom_200209_201808.pdf');
    figure(4); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_wv_anom_200209_201808.pdf');
    figure(6); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_o3_anom_200209_201808.pdf');
    figure(2); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_stemp_anom_200209_201808.pdf');

    figure(8); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_ptemp_rate_200209_201808.pdf');
    figure(9); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_wv_rate_200209_201808.pdf');
    figure(10); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_obs_o3_rate_200209_201808.pdf');

  else
    figure(5); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_ptemp_anom_200209_201808.pdf');
    figure(4); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_wv_anom_200209_201808.pdf');
    figure(6); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_o3_anom_200209_201808.pdf');
    figure(2); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_stemp_anom_200209_201808.pdf');

    figure(8); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_ptemp_rate_200209_201808.pdf');
    figure(9); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_wv_rate_200209_201808.pdf');
    figure(10); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sep19/Figs/ClearAnom/umbc_clr_retr_cal_o3_rate_200209_201808.pdf');

  end

end

figure(12); colormap(llsmap4.llsmap4);
figure(13); colormap(llsmap4.llsmap4);
i400 = abs(playsN-400); i400 = find(i400 == min(i400));
tzanom400    = squeeze(tz(:,:,i400));
tzanomsig400 = squeeze(tz_sig(:,:,i400));
figure(12); pcolor(okdates,latbins,tzanom400);          shading interp; title('Tz anom at 400 mb'); caxis([-5 +5]); colorbar
figure(13); pcolor(okdates,latbins,real(tzanomsig400)); shading interp; title('Tz unc anom at 400 mb'); caxis([0 0.1]); colorbar

figure(14); colormap(llsmap4.llsmap4);
figure(15); colormap(llsmap4.llsmap4);
i400 = abs(playsN-400); i400 = find(i400 == min(i400));
wvanom400    = squeeze(wv(:,:,i400));
wvanomsig400 = squeeze(wv_sig(:,:,i400));
figure(14); pcolor(okdates,latbins,wvanom400);          shading interp; title('Wv anom at 400 mb'); caxis([-1 +1]); colorbar
figure(15); pcolor(okdates,latbins,real(wvanomsig400)); shading interp; title('Wv unc anom at 400 mb'); caxis([0 0.005]); colorbar

