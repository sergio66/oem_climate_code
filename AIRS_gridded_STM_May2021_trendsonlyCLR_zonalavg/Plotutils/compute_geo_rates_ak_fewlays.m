fnameERArates = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS_ANOM_STM_Sept2016/LatbinProfRates/all40latbins.mat';
fnameERArates = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Mar29_2019/Desc/all_latbins_rates.mat';
fnameERArates = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Mar29_2019_Clr/Desc/all_latbins_rates.mat';

fnameERArtp = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Mar29_2019_Clr/Desc/xlatbin_oneyear_1_40.op.rtp';

disp(['for ERA rates using ' fnameERArates])
disp(['mean ERA profile    ' fnameERArtp])
clear thestats xthestats
xthestats = load(fnameERArates);

addpath /asl/matlib/h4tools
[hjunk,ha,pjunk,pa] = rtpread(fnameERArtp);

i1231 = find(fairs >= 1231.15,1);

  thestats.waterrate0 = xthestats.thestats.waterrate;
  thestats.ptemprate0 = xthestats.thestats.ptemprate;
  thestats.ozonerate0 = xthestats.thestats.ozonerate;
  thestats.waterratestd_lag0 = xthestats.thestats.waterratestd;
  thestats.ptempratestd_lag0 = xthestats.thestats.ptempratestd;
  thestats.ozoneratestd_lag0 = xthestats.thestats.ozoneratestd;    
  thestats.stemprate0    = xthestats.thestats.stemprate;
  thestats.stempratestd0 = xthestats.thestats.stempratestd;    
  
%% true profile = AK * Retr_profile +   (I - AK) apriori_profile
for ii = 1 : 40
  fname = [outputdir '/test' num2str(ii) '.mat'];
  loader = ['a = load(''' fname ''');'];
  eval(loader)
  clear ak_*

  ak_wv    = a.oem.ak_water;
  ak_t     = a.oem.ak_temp;
  ak_ozone = a.oem.ak_ozone;

  if ii == 1
    for jj = 1 : length(driver.jacobian.wvjaclays_used)
      junk = driver.jacobian.wvjaclays_used{jj}-6;
      thestats.waterrate1(:,jj)        = mean(thestats.waterrate0(:,junk)');
      thestats.waterratestd_lag1(:,jj) = mean(thestats.waterratestd_lag0(:,junk)');
      thestats.ptemprate1(:,jj)        = mean(thestats.ptemprate0(:,junk)');
      thestats.ptempratestd_lag1(:,jj) = mean(thestats.ptempratestd_lag0(:,junk)');
      thestats.ozonerate1(:,jj)        = mean(thestats.ozonerate0(:,junk)');
      thestats.ozoneratestd_lag1(:,jj) = mean(thestats.ozoneratestd_lag0(:,junk)');

      thestats.mean_era_ptemp_0(:,jj) = mean(pjunk.ptemp(junk,:));
      thestats.mean_era_gas_1_0(:,jj) = mean(pjunk.gas_1(junk,:));
      thestats.mean_era_gas_3_0(:,jj) = mean(pjunk.gas_3(junk,:));
    end
  end

  thestats.akwaterrate(ii,:) = (ak_wv    * thestats.waterrate1(ii,:)')';
  thestats.akptemprate(ii,:) = (ak_t     * thestats.ptemprate1(ii,:)')';
  thestats.akozonerate(ii,:) = (ak_ozone * thestats.ozonerate1(ii,:)')';    

  thestats.akwaterrate_stdlag(ii,:)  = (ak_wv    * thestats.waterratestd_lag1(ii,:)')';
  thestats.akptemprate_stdlag(ii,:)  = (ak_t     * thestats.ptempratestd_lag1(ii,:)')';
  thestats.akozonerate_stdlag(ii,:)  = (ak_ozone * thestats.ozoneratestd_lag1(ii,:)')';

  thestats.rate1231(ii)    = a.rateset.rates(i1231);
  thestats.rate1231unc(ii) = a.rateset.unc_rates(i1231);

  %% see how profiles change using ERA fractional rates
  thestats.mean_era_ptemp_F(ii,:) = thestats.mean_era_ptemp_0(ii,:) + thestats.akptemprate(ii,:);
  thestats.mean_era_gas_1_F(ii,:) = thestats.mean_era_gas_1_0(ii,:) .* ( 1 + thestats.akwaterrate(ii,:));
  thestats.mean_era_gas_3_F(ii,:) = thestats.mean_era_gas_3_0(ii,:) .* ( 1 + thestats.akozonerate(ii,:));

  %% see how profiles change using Sergio retrieved rates
  thestats.mean_era_ptemp_S(ii,:) = thestats.mean_era_ptemp_0(ii,:) + a.oem.finalrates(a.jacobian.temp_i)';
  thestats.mean_era_gas_1_S(ii,:) = thestats.mean_era_gas_1_0(ii,:) .* ( 1 + a.oem.finalrates(a.jacobian.water_i)');
  thestats.mean_era_gas_3_S(ii,:) = thestats.mean_era_gas_3_0(ii,:) .* ( 1 + a.oem.finalrates(a.jacobian.ozone_i)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9); clf
plot(1 - thestats.mean_era_gas_1_0 ./ thestats.mean_era_gas_1_F,playsRET,'b',1 - thestats.mean_era_gas_1_0 ./ thestats.mean_era_gas_1_S,playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.025 +0.025 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('frac \delta WV(z) (b) ERA (r) UMBC')

figure(10); clf
plot(thestats.mean_era_ptemp_0-thestats.mean_era_ptemp_F,playsRET,'b',thestats.mean_era_ptemp_0-thestats.mean_era_ptemp_S,playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.2 +0.2 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('\delta T(z) (b) ERA (r) UMBC')

figure(11); clf
plot(1 - thestats.mean_era_gas_3_0 ./ thestats.mean_era_gas_3_F,playsRET,'b',1 - thestats.mean_era_gas_3_0 ./ thestats.mean_era_gas_3_S,playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.1 +0.1 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('frac \delta O3(z) (b) ERA (r) UMBC')

figure(9); clf
plot(1 - thestats.mean_era_gas_1_0(22,:) ./ thestats.mean_era_gas_1_F(22,:),playsRET,'b',1 - thestats.mean_era_gas_1_0(22,:) ./ thestats.mean_era_gas_1_S(22,:),playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.025 +0.025 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('frac \delta WV(z) (b) ERA (r) UMBC')

figure(10); clf
plot(thestats.mean_era_ptemp_0(22,:)-thestats.mean_era_ptemp_F(22,:),playsRET,'b',thestats.mean_era_ptemp_0(22,:)-thestats.mean_era_ptemp_S(22,:),playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.2 +0.2 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('\delta T(z) (b) ERA (r) UMBC')

figure(11); clf
plot(1 - thestats.mean_era_gas_3_0(22,:) ./ thestats.mean_era_gas_3_F(22,:),playsRET,'b',1 - thestats.mean_era_gas_3_0(22,:) ./ thestats.mean_era_gas_3_S(22,:),playsRET,'r')
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-0.1 +0.1 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  title('frac \delta O3(z) (b) ERA (r) UMBC')

%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to go on ERA(lat,z)')
pause

figure(9); clf
  pcolor(latx,playsRET,thestats.waterrate1'*10);
  caxis([-0.1 +0.1])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA  WV(lat,z)')
  
figure(10); clf
  pcolor(latx,playsRET,thestats.ptemprate1'*10);
  caxis([-0.5 +0.5])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA  T(lat,z)')
  
figure(11); clf
  pcolor(latx,playsRET,thestats.ozonerate1'*10);
  caxis([-0.2 +0.2])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA  O3(lat,z)')

%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to go on AK * ERA(lat,z)')
pause

figure(9); clf
  pcolor(latx,playsRET,thestats.akwaterrate'*10);
  caxis([-0.1 +0.1])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * WV(lat,z)')
  
figure(10); clf
  pcolor(latx,playsRET,thestats.akptemprate'*10);
  caxis([-0.5 +0.5])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * T(lat,z)')
  
figure(11); clf
  pcolor(latx,playsRET,thestats.akozonerate'*10);
  caxis([-0.2 +0.2])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * O3(lat,z)')

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12); clf
errorbar(latx,100*traceNstemp(6,:),100*traceNstemp_sigs(6,:),'b','linewidth',2); grid
hold on
errorbar(latx,100*thestats.stemprate0,100*thestats.stempratestd0,'color','r','linewidth',2)
errorbar(latx,100*thestats.rate1231,100*thestats.rate1231unc,'g','linewidth',2);
title('(b) UMBC (r)ERA stemp/century ')
hl = legend('UMBC stemp','ERA stemp','1231 cm-1 rate','location','south'); set(hl,'fontsize',10); grid
hold off
grid on

jett = jet; jett(1,:) = 1;

disp('skipping error bars for ERA std devs')
%{
figure(13); clf
  pcolor(latx,playsRET,thestats.akwaterrate_stdlag');
  caxis([0 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * WV(lat,z)stdlag')
  colormap(jett)
  
figure(14); clf
  pcolor(latx,playsRET,thestats.akptemprate_stdlag');
  caxis([0 +0.1])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * T(lat,z)stdlag')
  colormap(jett)
  
figure(11); clf
  pcolor(latx,playsRET,thestats.akozonerate'*10);
  caxis([-0.02 +0.02])
  llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})
  shading interp
  colorbar
  title('ERA AK * O3(lat,z)')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); figure(9); disp('ret'); pause
figure(2); figure(10); disp('ret'); pause
figure(3); figure(11); disp('ret'); pause
figure(4); figure(12); disp('ret'); pause
figure(4); axis([-90 +90 -0.5 +0.5]);
%figure(10); axis([-90 +90 -0.5 +0.5]);

