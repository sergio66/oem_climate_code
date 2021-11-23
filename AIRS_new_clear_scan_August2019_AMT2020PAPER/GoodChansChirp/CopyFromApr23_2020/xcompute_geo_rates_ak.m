% load('/home/sergio/MATLABCODE/ROSES_2013_GRANT/NewClearPDFs/era_geo_rates.mat');
% load('/home/sergio/MATLABCODE/ROSES_2013_GRANT/NewClearPDFs/era_geo_rates_allsky_AIRIBRAD.mat');
%load('/home/sergio/MATLABCODE/ROSES_2013_GRANT/NewClearPDFs/daily_era_geo_rates_allsky_AIRIBRAD.mat');

disp('for ERA using /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS_ANOM_STM_Sept2016/LatbinProfRates/all40latbins.mat')
clear thestats xthestats
xthestats = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS_ANOM_STM_Sept2016/LatbinProfRates/all40latbins.mat');
  thestats.waterrate = xthestats.all.waterrate;
  thestats.ptemprate = xthestats.all.ptemprate;
  thestats.waterratestd_lag = xthestats.all.waterratestd_lag;
  thestats.ptempratestd_lag = xthestats.all.ptempratestd_lag;
%  xtats.ozonerate    = xthestats.all.ozonerate;
%  xtats.ozoneratestd = xthestats.all.ozoneratestd;    
  xstats.ozonerate    = xthestats.all.waterrate * 0;
  xstats.ozoneratestd = xthestats.all.waterratestd * 0;    
  thestats.stemprate    = xthestats.all.stemprate;
  thestats.stempratestd = xthestats.all.stempratestd;    
  
%% true profile = AK * Retr_profile +   (I - AK) apriori_profile
for ii = 1 : 40
  fname = ['/home/sergio//MATLABCODE/oem_pkg_run/AIRS_AllSky_97_O3jacs/Output/test' num2str(ii) '.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader)
  clear ak_*
  %{
  ak_wv = zeros(100,100);    ak_wv(1:97,1:97)    = a.oem.ak_water;
  ak_t = zeros(100,100);     ak_t(1:97,1:97)     = a.oem.ak_temp;
  ak_ozone = zeros(100,100); ak_ozone(1:97,1:97) = a.oem.ak_ozone;
  
  thestats.akwaterrate(ii,:) = ak_wv    * thestats.waterrate(ii,:)';
  thestats.akptemprate(ii,:) = ak_t     * thestats.ptemprate(ii,:)';
  thestats.akozonerate(ii,:) = ak_ozone * xstats.ozonerate(ii,:)';    
  thestats.ozonerate(ii,:)    = xstats.ozonerate(ii,:);
  thestats.ozoneratestd(ii,:) = xstats.ozoneratestd(ii,:);
  thestats.mmwrate(ii)        = xstats.mmwrate(ii);
  thestats.mmwratestd(ii)     = xstats.mmwratestd(ii);    
  %}
  
  ak_wv    = a.oem.ak_water;
  ak_t     = a.oem.ak_temp;
  ak_ozone = a.oem.ak_ozone;
  thestats.akwaterrate(ii,:) = ak_wv    * thestats.waterrate(ii,1:97)';
  thestats.akptemprate(ii,:) = ak_t     * thestats.ptemprate(ii,1:97)';
  thestats.akozonerate(ii,:) = ak_ozone * xstats.ozonerate(ii,1:97)';    

  thestats.akwaterrate_stdlag(ii,:)  = thestats.waterratestd_lag(ii,1:97)';
  thestats.akptemprate_stdlag(ii,:)  = thestats.ptempratestd_lag(ii,1:97)';
  %thestats.akwaterrate_stdlag(ii,:)  = ak_wv    * thestats.waterratestd_lag(ii,1:97)';
  %thestats.akptemprate_stdlag(ii,:)  = ak_t     * thestats.ptempratestd_lag(ii,1:97)';
  %%%thestats.akozonerate_stdlag(ii,:) = ak_ozone * xstats.ozoneratestd_lag(ii,1:97)';    

  thestats.ozonerate(ii,:)    = xstats.ozonerate(ii,:);
  thestats.ozoneratestd(ii,:) = xstats.ozoneratestd(ii,:);
  thestats.mmwrate(ii)        = xstats.mmwrate(ii);
  thestats.mmwratestd(ii)     = xstats.mmwratestd(ii);    
end

%{
figure(1); clf
  pcolor(latx,airslevels,thestats.waterrate');
  caxis([-0.01 +0.01])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('raw WV(lat,z)')
  
figure(2); clf
  pcolor(latx,airslevels,thestats.ptemprate');
  caxis([-0.15 +0.15])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('raw T(lat,z)')
  
figure(3); clf
  pcolor(latx,airslevels,xstats.ozonerate');
  caxis([-0.02 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('raw O3(lat,z)')

figure(4);
errorbar(latx,thestats.stemprate,thestats.stempratestd,'color','b','linewidth',2)
hold on
errorbar(latx,xstats.mmwrate,xstats.mmwratestd,'color','r','linewidth',2)
hold off
hl = legend('stemp','mmw','location','south'); set(hl,'fontsize',10); grid
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airslevels97 = airslevels(1:97);
figure(7); clf
  pcolor(latx,airslevels97,thestats.akwaterrate');
  caxis([-0.01 +0.01])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * WV(lat,z)')
  
figure(8); clf
  pcolor(latx,airslevels97,thestats.akptemprate');
  caxis([-0.15 +0.15])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * T(lat,z)')
  
figure(9); clf
  pcolor(latx,airslevels97,thestats.akozonerate');
  caxis([-0.02 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * O3(lat,z)')

figure(10);
errorbar(latx,thestats.stemprate,thestats.stempratestd,'color','b','linewidth',2)
hold on
errorbar(latx,xstats.mmwrate,xstats.mmwratestd,'color','r','linewidth',2)
hold off
title('ERA')
hl = legend('stemp','mmw','location','south'); set(hl,'fontsize',10); grid

jett = jet; jett(1,:) = 1;

figure(11); clf
  pcolor(latx,airslevels97,thestats.akwaterrate_stdlag');
  caxis([0 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * WV(lat,z)stdlag')
  colormap(jett)
  
figure(12); clf
  pcolor(latx,airslevels97,thestats.akptemprate_stdlag');
  caxis([0 +0.1])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * T(lat,z)stdlag')
  colormap(jett)
  
figure(9); clf
  pcolor(latx,airslevels97,thestats.akozonerate');
  caxis([-0.02 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading flat
  colorbar
  title('ERA AK * O3(lat,z)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1); figure(7); disp('ret'); pause
figure(2); figure(8); disp('ret'); pause
figure(3); figure(9); disp('ret'); pause
figure(4); figure(10); disp('ret'); pause
figure(4); axis([-90 +90 -0.5 +0.5]);
figure(10); axis([-90 +90 -0.5 +0.5]);
%}
