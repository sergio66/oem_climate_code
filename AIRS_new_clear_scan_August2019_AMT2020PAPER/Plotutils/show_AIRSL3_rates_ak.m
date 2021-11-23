disp('for AIRS L3 showing /asl/s1/sergio/AIRS_L3/airsL3_v6_rates_stats_Mar2016_13yr.mat')
airsL2 = load('/asl/s1/sergio/AIRS_L3/airsL3_v6_rates_stats_Mar2016_13yr.mat');

clear thestats

save_lat = airsL2.save_lat;
iXStart = 1;
iXEnd   = 40;

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
lala = load('llsmap5'); 

Airs_PQ = airsL2.Qlevs;
Airs_PT = airsL2.Tlevs;

playsX = plays; plays(end) = 1001;
%% true profile = AK * Retr_profile +   (I - AK) apriori_profile
[xplays2,yplays2]   = meshgrid(playsX);
[xAirs_PQ2,yAirs_PQ2] = meshgrid(Airs_PQ);
[xAirs_PT2,yAirs_PT2] = meshgrid(Airs_PT);

for ii = 1 : 40
  fname = ['/home/sergio//MATLABCODE/oem_pkg_run/AIRS_AllSky_97_O3jacs/Output/test' num2str(ii) '.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader)
  clear ak_*

  ak_wv    = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PQ2),log10(yAirs_PQ2));
  ak_t     = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PT2),log10(yAirs_PT2));
  ak_ozone = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PT2),log10(yAirs_PT2));

  ak_wv(isnan(ak_wv)) = 0;  ak_wv(isinf(ak_wv)) = 0;
  ak_t(isnan(ak_t)) = 0;  ak_t(isinf(ak_t)) = 0;
  ak_ozone(isnan(ak_ozone)) = 0;  ak_ozone(isinf(ak_ozone)) = 0;
  
  thestatsL2.akwaterrate(ii,:) = ak_wv    * airsL2.thestats.waterrate(ii,:)';
  thestatsL2.akptemprate(ii,:) = ak_t     * airsL2.thestats.ptemprate(ii,:)';
  thestatsL2.akozonerate(ii,:) = ak_ozone * airsL2.thestats_other.ozonerate(ii,:)';    

  thestatsL2.akwaterrate_stdlag(ii,:)  = airsL2.thestats.waterratestd_lag(ii,:)';
  thestatsL2.akptemprate_stdlag(ii,:)  = airsL2.thestats.ptempratestd_lag(ii,:)';
  %thestatsL2.akwaterrate_stdlag(ii,:)  = ak_wv    * airsL2.thestats.waterratestd_lag(ii,:)';
  %thestatsL2.akptemprate_stdlag(ii,:)  = ak_t     * airsL2.thestats.ptempratestd_lag(ii,:)';
  %%%thestatsL2.akozonerate_stdlag(ii,:) = ak_ozone * airsL2.thestats_ther.ozoneratestd_lag(ii,:)';    

%  thestatsL2.ozonerate(ii,:)    = xstats.ozonerate(ii,:);
%  thestatsL2.ozoneratestd(ii,:) = xstats.ozoneratestd(ii,:);
%  thestatsL2.mmwrate(ii)        = xstats.mmwrate(ii);
%  thestatsL2.mmwratestd(ii)     = xstats.mmwratestd(ii);    
end

[x,waterlen] = size(thestatsL2.akwaterrate);
figure(7); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PQ),double(thestatsL2.akwaterrate(iXStart:iXEnd,:)')); 
  title('AIRS L3 water rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.01 +0.01]); colorbar
colormap(lala.llsmap5)

[x,ptemplen] = size(thestatsL2.akptemprate);
figure(8); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestatsL2.akptemprate(iXStart:iXEnd,:)')); 
  title('AIRS L3 tempr rate (K/yr)')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp; colorbar
  caxis([-0.15 +0.15]); colorbar
colormap(lala.llsmap5)

[x,ozonelen] = size(thestatsL2.akozonerate);
figure(9); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),double(thestatsL2.akozonerate(iXStart:iXEnd,:)')); 
  title('AIRS L3 ozone rate (frac/yr)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.02 +0.02]); colorbar
colormap(lala.llsmap5)

figure(10); clf;
plot(save_lat(iXStart:iXEnd),airsL2.thestats.stemprate,'k',...
     save_lat(iXStart:iXEnd),airsL2.thestats_other.olrrate,'r',...
     save_lat(iXStart:iXEnd),airsL2.thestats_other.clrolrrate,'b'); hold on
errorbar(save_lat(iXStart:iXEnd),airsL2.thestats.stemprate,airsL2.thestats.stempratestd,'color','k','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),airsL2.thestats_other.olrrate,airsL2.thestats_other.olrratestd,'color','r','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),airsL2.thestats_other.clrolrrate,airsL2.thestats_other.clrolrratestd,'color','b','linewidth',2); hold off
xlabel('latitude'); hl = legend('Stemp','OLR','Clr OLR','location','north'); set(hl,'fontsize',10);
title('AIRS L3 rates in K/yr and W/m2/yr')

figure(11);
  plot(save_lat(iXStart:iXEnd),airsL2.thestats_other.ice_od_rate,'b',save_lat(iXStart:iXEnd),airsL2.thestats_other.liq_water_rate,'r',...
       save_lat(iXStart:iXEnd),airsL2.thestats_other.icesze_rate,'c','linewidth',2)
  title('AIRS L3 fractional cloud rates'); grid       
  hl = legend('ice OD frac','liq water frac','ice size frac');
