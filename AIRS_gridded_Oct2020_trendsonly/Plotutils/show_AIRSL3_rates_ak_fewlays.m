disp('for AIRS L3 showing /asl/s1/sergio/AIRS_L3/airsL3_v6_rates_stats_Mar2016_13yr.mat')
airsL3 = load('/asl/s1/sergio/AIRS_L3/airsL3_v6_rates_stats_Mar2016_13yr.mat');

clear thestats

save_lat = airsL3.save_lat;
iXStart = 1;
iXEnd   = 40;

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
lala = load('llsmap5'); 

Airs_PQ = airsL3.Qlevs;
Airs_PT = airsL3.Tlevs;

playsX = plays;    playsX(end) = 1001;
playsX = playsRET; playsX(end) = 1001;
%% true profile = AK * Retr_profile +   (I - AK) apriori_profile
[xplays2,yplays2]   = meshgrid(playsX);
[xAirs_PQ2,yAirs_PQ2] = meshgrid(Airs_PQ);
[xAirs_PT2,yAirs_PT2] = meshgrid(Airs_PT);

iMultAK = -1;
iMultAK = +1;

for ii = 1 : 40
  fname = ['/home/sergio//MATLABCODE/oem_pkg_run/AIRS_AllSky_97_O3jacs/Output/test' num2str(ii) '.mat'];
  fname = [outputdir '/test' num2str(ii) '.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader)
  clear ak_*

  ak_wv    = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PQ2),log10(yAirs_PQ2));
  ak_t     = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PT2),log10(yAirs_PT2));
  ak_ozone = interp2(log10(xplays2),log10(yplays2),a.oem.ak_water,log10(xAirs_PT2),log10(yAirs_PT2));

  ak_wv(isnan(ak_wv)) = 0;  ak_wv(isinf(ak_wv)) = 0;
  ak_t(isnan(ak_t)) = 0;  ak_t(isinf(ak_t)) = 0;
  ak_ozone(isnan(ak_ozone)) = 0;  ak_ozone(isinf(ak_ozone)) = 0;

  if iMultAK > 0  
    thestatsL2.akwaterrate(ii,:) = ak_wv    * airsL3.thestats.waterrate(ii,:)';
    thestatsL2.akptemprate(ii,:) = ak_t     * airsL3.thestats.ptemprate(ii,:)';
    thestatsL2.akozonerate(ii,:) = ak_ozone * airsL3.thestats_other.ozonerate(ii,:)';    
  else
    thestatsL2.akwaterrate(ii,:) = airsL3.thestats.waterrate(ii,:)';
    thestatsL2.akptemprate(ii,:) = airsL3.thestats.ptemprate(ii,:)';
    thestatsL2.akozonerate(ii,:) = airsL3.thestats_other.ozonerate(ii,:)';    
  end

  thestatsL2.akwaterrate_stdlag(ii,:)  = airsL3.thestats.waterratestd_lag(ii,:)';
  thestatsL2.akptemprate_stdlag(ii,:)  = airsL3.thestats.ptempratestd_lag(ii,:)';
  %thestatsL2.akwaterrate_stdlag(ii,:)  = ak_wv    * airsL3.thestats.waterratestd_lag(ii,:)';
  %thestatsL2.akptemprate_stdlag(ii,:)  = ak_t     * airsL3.thestats.ptempratestd_lag(ii,:)';
  %%%thestatsL2.akozonerate_stdlag(ii,:) = ak_ozone * airsL3.thestats_ther.ozoneratestd_lag(ii,:)';    

end

[x,waterlen] = size(thestatsL2.akwaterrate);
figure(7); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PQ),10*double(thestatsL2.akwaterrate(iXStart:iXEnd,:)')); 
  title('AIRS L3 water rate (frac/decade)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.1 +0.1]); colorbar
colormap(lala.llsmap5)

[x,ptemplen] = size(thestatsL2.akptemprate);
figure(8); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),10*double(thestatsL2.akptemprate(iXStart:iXEnd,:)')); 
  title('AIRS L3 tempr rate (K/decade)')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp; colorbar
  caxis([-0.5 +0.5]); colorbar
colormap(lala.llsmap5)

[x,ozonelen] = size(thestatsL2.akozonerate);
figure(9); clf
  pcolor(save_lat(iXStart:iXEnd),(Airs_PT),10*double(thestatsL2.akozonerate(iXStart:iXEnd,:)')); 
  title('AIRS L3 ozone rate (frac/decade)')
  set(gca,'ydir','reverse');
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})    
  shading interp; colorbar
  caxis([-0.2 +0.2]); colorbar
colormap(lala.llsmap5)

figure(10); clf;
plot(save_lat(iXStart:iXEnd),airsL3.thestats.stemprate,'k',...
     save_lat(iXStart:iXEnd),airsL3.thestats_other.olrrate,'r',...
     save_lat(iXStart:iXEnd),airsL3.thestats_other.clrolrrate,'b'); hold on
errorbar(save_lat(iXStart:iXEnd),airsL3.thestats.stemprate,airsL3.thestats.stempratestd,'color','k','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),airsL3.thestats_other.olrrate,airsL3.thestats_other.olrratestd,'color','r','linewidth',2); hold on
errorbar(save_lat(iXStart:iXEnd),airsL3.thestats_other.clrolrrate,airsL3.thestats_other.clrolrratestd,'color','b','linewidth',2); hold off
xlabel('latitude'); hl = legend('Stemp','OLR','Clr OLR','location','north'); set(hl,'fontsize',10);
title('AIRS L3 rates in K/yr and W/m2/yr')

figure(11);
  plot(save_lat(iXStart:iXEnd),airsL3.thestats_other.ice_od_rate,'b',save_lat(iXStart:iXEnd),airsL3.thestats_other.liq_water_rate,'r',...
       save_lat(iXStart:iXEnd),airsL3.thestats_other.icesze_rate,'c','linewidth',2)
  title('AIRS L3 fractional cloud rates'); grid       
  hl = legend('ice OD frac','liq water frac','ice size frac');
