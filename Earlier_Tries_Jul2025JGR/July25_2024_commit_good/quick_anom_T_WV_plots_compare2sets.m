jett = jet; jett(1,:) = 1;
llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); 

disp('ret to see T(z),WV(z),O3(z)'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(4); pcolor(okdates,playsN,smoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical WV smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar; 
  colormap(llsmap4.llsmap4);

  figure(5); pcolor(okdates,playsN,smoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical Tz smoothed over 10 years')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  figure(6); pcolor(okdates,playsN,smootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical O3 smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar
  colormap(llsmap4.llsmap4);

disp('ret to see xT(z),xWV(z),xO3(z)'); pause
  figure(7); plot(okdates,smooth(xmeanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(8); pcolor(okdates,playsN,xsmoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical WV smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar; 
  colormap(llsmap4.llsmap4);

  figure(9); pcolor(okdates,playsN,xsmoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical Tz smoothed over 10 years')
  caxis([-1 +1]); colorbar
  colormap(llsmap4.llsmap4);

  figure(10); pcolor(okdates,playsN,xsmootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical O3 smoothed over 10 years')
  caxis([-0.2 +0.2]); colorbar
  colormap(llsmap4.llsmap4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('ret to see T(z),WV(z),O3(z) decadal trends'); pause
  figure(3); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(4); pcolor(latbins,playsN,10*trendwv2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal WV trends')
  caxis([-0.01 +0.01]*10); colorbar; 

  figure(5); pcolor(latbins,playsN,10*trendtz2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends')
  caxis([-0.05 +0.05]*10); colorbar

  figure(6); pcolor(latbins,playsN,10*trendo32'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal O3 trends')
  caxis([-0.02 +0.02]*10); colorbar

for ii = 4 : 6
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

%%%%%%%%%%%%
  figure(7); plot(okdates,smooth(xmeanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical xstemp smoothed over 10 years')

  figure(8); pcolor(latbins,playsN,10*xtrendwv2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal xWV trends')
  caxis([-0.01 +0.01]*10); colorbar; 

  figure(9); pcolor(latbins,playsN,10*xtrendtz2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal xTz trends')
  caxis([-0.05 +0.05]*10); colorbar

  figure(10); pcolor(latbins,playsN,10*xtrendo32'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal xO3 trends')
  caxis([-0.02 +0.02]*10); colorbar

for ii = 8 : 10
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

