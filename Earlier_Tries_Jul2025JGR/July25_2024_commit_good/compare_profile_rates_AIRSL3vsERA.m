%addpath ~/Matlab/Sergio/COLORMAP/LLS
addpath ~/Matlab/Sergio/PLOTTER

% for ii = 1 : 10; figure(ii); colormap jet; end
load ~/Matlab/Cmaps/old_llsmap5.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ../Data/all40latbins.mat

figure(1); clf; pcolor(latbins,airsp,all.waterrate'); colorbar; title('water rate ERA frac/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([-0.01 +0.01]); colorbar;

figure(2); clf; pcolor(latbins,airsp,all.waterratestd_lag'); colorbar; title('water rate ERA err LAG frac/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([0 4e-3]); colorbar
  
figure(3); clf; pcolor(latbins,airsp,all.ptemprate'); colorbar; title('ptemp rate ERA K/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);  
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([-0.15 +0.15]); colorbar;
  
figure(4); clf; pcolor(latbins,airsp,all.ptempratestd_lag'); colorbar; title('ptemp rate ERA err LAG K/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([0 0.15]); colorbar
  
figure(5); clf; errorbar(latbins,all.stemprate,all.stempratestd_lag,'color','blue','linewidth',2)
grid; title('Stemp rate ERA K/yr'); ax = axis; axis([-90 +90 ax(3) ax(4)])
axis([-90 +90 -0.05 +0.05])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ../Data/airsL3_v6_rates_stats_Mar2016_13yr.mat

figure(6); clf; pcolor(thestats.lats,Qlevs,thestats.waterrate'); colorbar; title('water rate AIRS L3 frac/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([-0.01 +0.01]); colorbar;

figure(7); clf; pcolor(thestats.lats,Qlevs,thestats.waterratestd_lag'); colorbar; title('water rate AIRS L3 err LAG frac/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([0 4e-3]); colorbar
  
figure(8); clf; pcolor(thestats.lats,Tlevs,thestats.ptemprate'); colorbar; title('ptemp rate AIRS L3 K/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);  
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([-0.15 +0.15]); colorbar;

figure(9); clf; pcolor(thestats.lats,Tlevs,thestats.ptempratestd_lag'); colorbar; title('ptemp rate AIRS L3 err LAG K/yr')
  colormap(llsmap5); shading flat; set(gca,'ydir','reverse'); axis([-90 +90 1 3]);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  caxis([0 0.15]); colorbar
  
figure(10); clf; errorbar(thestats.lats,thestats.stemprate,thestats.stempratestd_lag,'color','blue','linewidth',2)
grid; title('Stemp rate AIRS L3 K/yr'); ax = axis; axis([-90 +90 ax(3) ax(4)])
axis([-90 +90 -0.05 +0.05])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for ii = 1 : 5
%   figure(ii); figure(ii+5); disp('ret'); pause
% end
