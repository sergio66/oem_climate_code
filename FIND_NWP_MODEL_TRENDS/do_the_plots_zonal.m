addpath /home/sergio/MATLABCODE/COLORMAP/LLS
lala = load('llsmap5'); 

Airs_PQ = Qlevs;
Airs_PT = Tlevs;

[x,waterlen] = size(thestats.waterrate);
figure(1); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PQ),double(thestats.waterrate(iXStart:iXEnd,:)')); 
  title('water rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.01 +0.01]); colorbar
colormap(lala.llsmap5)

[x,ptemplen] = size(thestats.ptemprate);
figure(2); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestats.ptemprate(iXStart:iXEnd,:)')); 
  title('tempr rate (K/yr)')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp; colorbar
  caxis([-0.15 +0.15]); colorbar
colormap(lala.llsmap5)

[x,ozonelen] = size(thestats.ozonerate);
figure(3); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestats.ozonerate(iXStart:iXEnd,:)')); 
  title('ozone rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.02 +0.02]); colorbar
colormap(lala.llsmap5)

[x,ch4len] = size(thestats.ch4rate);
figure(4); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestats.ch4rate(iXStart:iXEnd,:)')); 
  title('ch4 rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.01 +0.01]); colorbar
colormap(lala.llsmap5)

[x,colen] = size(thestats.corate);
figure(5); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestats.corate(iXStart:iXEnd,:)')); 
  title('co rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.01 +0.01]); colorbar
colormap(lala.llsmap5)

figure(6); clf;
plot(save_lat(iXStart:iXEnd),thestats.stemprate,'k',...
     save_lat(iXStart:iXEnd),thestats_other.olrrate,'r',...
     save_lat(iXStart:iXEnd),thestats_other.clrolrrate,'b','linewidth',2); grid on
xlabel('latitude'); hl = legend('Stemp','OLR','Clr OLR','location','north'); set(hl,'fontsize',10);
title('rates in K/yr and W/m2/yr')

figure(7); clf;
plot(save_lat(iXStart:iXEnd),thestats.stemprate,'k',...
     save_lat(iXStart:iXEnd),thestats_other.olrrate,'r',...
     save_lat(iXStart:iXEnd),thestats_other.clrolrrate,'b'); hold on
errorbar(save_lat(iXStart:iXEnd),thestats.stemprate,thestats.stempratestd,'color','k','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),thestats_other.olrrate,thestats_other.olrratestd,'color','r','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),thestats_other.clrolrrate,thestats_other.clrolrratestd,'color','b','linewidth',2); hold off
xlabel('latitude'); hl = legend('Stemp','OLR','Clr OLR','location','north'); set(hl,'fontsize',10);
title('rates in K/yr and W/m2/yr')

[x,waterlen] = size(thestats.waterrate);
figure(8); clf
  booW = thestats.waterrate(iXStart:iXEnd,:);
  booT = thestats.ptemprate(iXStart:iXEnd,:);  
  booWT = booW .* booT(:,1:12);
  pcolor(save_lat(iXStart:iXEnd),(Airs_PQ),double(booWT'))
  title('water * temp rate (K/yr * frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.001 +0.001]); colorbar
colormap(lala.llsmap5)

figure(9);
  plot(save_lat(iXStart:iXEnd),thestats_other.ice_od_rate,'b',save_lat(iXStart:iXEnd),thestats_other.liq_water_rate,'r',...
       save_lat(iXStart:iXEnd),thestats_other.icesze_rate,'c','linewidth',2)
  hl = legend('ice OD frac','liq water frac','ice size frac');
