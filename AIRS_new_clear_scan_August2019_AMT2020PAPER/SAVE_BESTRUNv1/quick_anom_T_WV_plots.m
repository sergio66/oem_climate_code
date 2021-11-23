disp('ret to see 400 mb T(lat),WV(lat),O3(lat) vs t'); pause

  figure(4); pcolor(okdates,latbins,smoothwv_400mb); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('400 mb WV smoothed over 0.5 years')
  caxis([-0.2 +0.2]); colorbar; 

  figure(5); pcolor(okdates,latbins,smoothtz_400mb); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('400 mb Tz smoothed over 0.5 years')
  caxis([-1 +1]); colorbar

  figure(6); pcolor(okdates,latbins,smootho3_400mb); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('400 mb O3 smoothed over 0.5 years')
  caxis([-0.2 +0.2]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to see T(z),WV(z),O3(z) vs t'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 0.5 years')

  figure(4); pcolor(okdates,playsN,smoothwv'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical WV smoothed over 0.5 years')
  caxis([-0.2 +0.2]); colorbar; 

  figure(5); pcolor(okdates,playsN,smoothtz'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical Tz smoothed over 0.5 years')
  caxis([-1 +1]); colorbar

  figure(6); pcolor(okdates,playsN,smootho3'); shading interp; colorbar
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('tropical O3 smoothed over 0.5 years')
  caxis([-0.2 +0.2]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to see T(z),WV(z),O3(z) decadal trends'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 0.5 years')

  figure(4); pcolor(latbins,playsN,10*trendwv2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal WV trends')
  caxis([-0.01 +0.01]*10); colorbar; 

  figure(5); pcolor(latbins,playsN,10*trendtz2'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends')
  caxis([-0.05 +0.05]*10); colorbar

  figure(6); pcolor(latbins,playsN,10*trendo32'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal O3 trends')
  caxis([-0.02 +0.02]*10); colorbar

jett = jet; jett(1,:) = 1;
for ii = 4 : 6
  figure(ii)
  colormap(jett)
  if ii <= 6
    colormap usa2  
    llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  end
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to see T(z),WV(z),O3(z) unc decadal trends'); pause
  figure(2); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 0.5 years')

  figure(4); pcolor(latbins,playsN,10*trendwv2_unc'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal WV trends unc')
  caxis([0 +0.01]*10); colorbar; 

  figure(5); pcolor(latbins,playsN,10*trendtz2_unc'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal Tz trends unc')
  caxis([0 +0.05]*10); colorbar

  figure(6); pcolor(latbins,playsN,10*trendo32_unc'); shading interp; colorbar
  ax = axis; ax(1) = min(latbins); ax(2) = max(latbins); ax(3) = 10; ax(4) = 1000;; grid; set(gca,'ydir','reverse'); title('decadal O3 trends unc')
  caxis([0 +0.02]*10); colorbar

jett = jet; jett(1,:) = 1;
for ii = 4 : 6
  figure(ii)
  colormap(jett)
  if ii <= 6
    colormap usa2  
    llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  end
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp
end
