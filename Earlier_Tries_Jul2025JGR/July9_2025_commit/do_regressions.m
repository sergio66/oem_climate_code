addpath /home/sergio/MATLABCODE/NANROUTINES/

figure(9); clf; 
% [r,chisqr,P,sigP,numpts] = nanlinearcorrelation(y1in,y2in)

[theregress.umbc_airsL3.stemp.r,theregress.umbc_airsL3.stemp.chisqr,theregress.umbc_airsL3.stemp.P,theregress.umbc_airsL3.stemp.sigP,theregress.umbc_airsL3.stemp.numpts] = ...
    nanlinearcorrelation(results(:,6)',reshape(airsL3.thestats64x72.stemprate,1,72*64));
[theregress.umbc_cmip6.stemp.r,theregress.umbc_cmip6.stemp.chisqr,theregress.umbc_cmip6.stemp.P,theregress.umbc_cmip6.stemp.sigP,theregress.umbc_cmip6.stemp.numpts] = nanlinearcorrelation(results(:,6)',cmip6.trend_stemp);
[theregress.umbc_era5.stemp.r,theregress.umbc_era5.stemp.chisqr,theregress.umbc_era5.stemp.P,theregress.umbc_era5.stemp.sigP,theregress.umbc_era5.stemp.numpts] = nanlinearcorrelation(results(:,6)',era5.trend_stemp);

for ii = 1 : 97
  [theregress.umbc_airsL3.ptemp.r(ii),theregress.umbc_airsL3.ptemp.chisqr(ii),theregress.umbc_airsL3.ptemp.P(ii,:),theregress.umbc_airsL3.ptemp.sigP(ii,:),theregress.umbc_airsL3.ptemp.numpts(ii,:)] = ...
    nanlinearcorrelation(deltaT(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(ii,:));
  [theregress.umbc_cmip6.ptemp.r(ii),theregress.umbc_cmip6.ptemp.chisqr(ii),theregress.umbc_cmip6.ptemp.P(ii,:),theregress.umbc_cmip6.ptemp.sigP(ii,:),theregress.umbc_cmip6.ptemp.numpts(ii,:)] = ...
    nanlinearcorrelation(deltaT(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(ii,:));
  [theregress.umbc_era5.ptemp.r(ii),theregress.umbc_era5.ptemp.chisqr(ii),theregress.umbc_era5.ptemp.P(ii,:),theregress.umbc_era5.ptemp.sigP(ii,:),theregress.umbc_era5.ptemp.numpts(ii,:)] = ...
    nanlinearcorrelation(deltaT(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(ii,:));
end

warning off
for ii = 1 : 97
  [theregress.umbc_airsL3.fracWV.r(ii),theregress.umbc_airsL3.fracWV.chisqr(ii),theregress.umbc_airsL3.fracWV.P(ii,:),theregress.umbc_airsL3.fracWV.sigP(ii,:),theregress.umbc_airsL3.fracWV.numpts(ii,:)] = ...
    nanlinearcorrelation(fracWV(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(ii,:));
  [theregress.umbc_cmip6.fracWV.r(ii),theregress.umbc_cmip6.fracWV.chisqr(ii),theregress.umbc_cmip6.fracWV.P(ii,:),theregress.umbc_cmip6.fracWV.sigP(ii,:),theregress.umbc_cmip6.fracWV.numpts(ii,:)] = ...
    nanlinearcorrelation(fracWV(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1(ii,:));
  [theregress.umbc_era5.fracWV.r(ii),theregress.umbc_era5.fracWV.chisqr(ii),theregress.umbc_era5.fracWV.P(ii,:),theregress.umbc_era5.fracWV.sigP(ii,:),theregress.umbc_era5.fracWV.numpts(ii,:)] = ...
    nanlinearcorrelation(fracWV(ii,:),nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.gas_1(ii,:));
end

boo1 = (compute_deltaRH.umbc.final(1:97,:)'-compute_deltaRH.umbc.orig(1:97,:)');
boo2 = (compute_deltaRH.airsL3.final(1:97,:)'-compute_deltaRH.airsL3.orig(1:97,:)');
boo3 = (compute_deltaRH.cmip6.final(1:97,:)'-compute_deltaRH.cmip6.orig(1:97,:)');
boo4 = (compute_deltaRH.era5.final(1:97,:)'-compute_deltaRH.era5.orig(1:97,:)');

for ii = 1 : 97
  [theregress.umbc_airsL3.deltaRH.r(ii),theregress.umbc_airsL3.deltaRH.chisqr(ii),theregress.umbc_airsL3.deltaRH.P(ii,:),theregress.umbc_airsL3.deltaRH.sigP(ii,:),theregress.umbc_airsL3.deltaRH.numpts(ii,:)] = ...
    nanlinearcorrelation(boo1(:,ii),boo2(:,ii));
  [theregress.umbc_cmip6.deltaRH.r(ii),theregress.umbc_cmip6.deltaRH.chisqr(ii),theregress.umbc_cmip6.deltaRH.P(ii,:),theregress.umbc_cmip6.deltaRH.sigP(ii,:),theregress.umbc_cmip6.deltaRH.numpts(ii,:)] = ...
    nanlinearcorrelation(boo1(:,ii),boo3(:,ii));
  [theregress.umbc_era5.deltaRH.r(ii),theregress.umbc_era5.deltaRH.chisqr(ii),theregress.umbc_era5.deltaRH.P(ii,:),theregress.umbc_era5.deltaRH.sigP(ii,:),theregress.umbc_era5.deltaRH.numpts(ii,:)] = ...
    nanlinearcorrelation(boo1(:,ii),boo4(:,ii));
end
warning on

figure(10);
  plot(theregress.umbc_airsL3.ptemp.r,playsjunk,theregress.umbc_cmip6.ptemp.r,playsjunk,theregress.umbc_era5.ptemp.r,playsjunk,'linewidth',2)
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dT(z)/dt')
figure(11);
  plot(theregress.umbc_airsL3.fracWV.r,playsjunk,theregress.umbc_cmip6.fracWV.r,playsjunk,theregress.umbc_era5.fracWV.r,playsjunk,'linewidth',2)
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dfracWV(z)/dt')
figure(12);
  plot(theregress.umbc_airsL3.deltaRH.r,playsjunk,theregress.umbc_cmip6.deltaRH.r,playsjunk,theregress.umbc_era5.deltaRH.r,playsjunk,'linewidth',2)
  set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dRH/dt')
for ii = 10:12;  figure(ii); ylim([10 1000]); end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to look zonally'); pause
bonk = reshape(airsL3.thestats64x72.stemprate,1,72*64);
for ll = 1 : 64
  ind = (1:72) + (ll-1)*72;
  [theregressY.umbc_airsL3.stemp.r(ll),theregressY.umbc_airsL3.stemp.chisqr(ll),theregressY.umbc_airsL3.stemp.P(ll,:),theregressY.umbc_airsL3.stemp.sigP(ll,:),theregressY.umbc_airsL3.stemp.numpts(ll,:)] = ...
    nanlinearcorrelation(results(ind,6)',bonk(ind));
  [theregressY.umbc_cmip6.stemp.r(ll),theregressY.umbc_cmip6.stemp.chisqr(ll),theregressY.umbc_cmip6.stemp.P(ll,:),theregressY.umbc_cmip6.stemp.sigP(ll,:),theregressY.umbc_cmip6.stemp.numpts(ll,:)] = ...
    nanlinearcorrelation(results(ind,6)',cmip6.trend_stemp(ind));
  [theregressY.umbc_era5.stemp.r(ll),theregressY.umbc_era5.stemp.chisqr(ll),theregressY.umbc_era5.stemp.P(ll,:),theregressY.umbc_era5.stemp.sigP(ll,:),theregressY.umbc_era5.stemp.numpts(ll,:)] = ...
    nanlinearcorrelation(results(ind,6)',era5.trend_stemp(ind));
end
figure(9); plot(theregressY.umbc_airsL3.stemp.r,rlat,theregressY.umbc_cmip6.stemp.r,rlat,theregressY.umbc_era5.stemp.r,rlat,'linewidth',2)
  hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('lat'); xlabel('R^2'); title('Correlation dStemp/dt')

for ll = 1 : 64
  ind = (1:72) + (ll-1)*72;
  for ii = 1 : 97
    [theregressY.umbc_airsL3.ptemp.r(ll,ii),theregressY.umbc_airsL3.ptemp.chisqr(ll,ii),theregressY.umbc_airsL3.ptemp.P(ll,ii,:),theregressY.umbc_airsL3.ptemp.sigP(ll,ii,:),theregressY.umbc_airsL3.ptemp.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(deltaT(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(ii,ind));
    [theregressY.umbc_cmip6.ptemp.r(ll,ii),theregressY.umbc_cmip6.ptemp.chisqr(ll,ii),theregressY.umbc_cmip6.ptemp.P(ll,ii,:),theregressY.umbc_cmip6.ptemp.sigP(ll,ii,:),theregressY.umbc_cmip6.ptemp.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(deltaT(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(ii,ind));
    [theregressY.umbc_era5.ptemp.r(ll,ii),theregressY.umbc_era5.ptemp.chisqr(ll,ii),theregressY.umbc_era5.ptemp.P(ll,ii,:),theregressY.umbc_era5.ptemp.sigP(ll,ii,:),theregressY.umbc_era5.ptemp.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(deltaT(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(ii,ind));

    [theregressY.umbc_airsL3.fracWV.r(ll,ii),theregressY.umbc_airsL3.fracWV.chisqr(ll,ii),theregressY.umbc_airsL3.fracWV.P(ll,ii,:),theregressY.umbc_airsL3.fracWV.sigP(ll,ii,:),theregressY.umbc_airsL3.fracWV.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(fracWV(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(ii,ind));
    [theregressY.umbc_cmip6.fracWV.r(ll,ii),theregressY.umbc_cmip6.fracWV.chisqr(ll,ii),theregressY.umbc_cmip6.fracWV.P(ll,ii,:),theregressY.umbc_cmip6.fracWV.sigP(ll,ii,:),theregressY.umbc_cmip6.fracWV.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(fracWV(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1(ii,ind));
    [theregressY.umbc_era5.fracWV.r(ll,ii),theregressY.umbc_era5.fracWV.chisqr(ll,ii),theregressY.umbc_era5.fracWV.P(ll,ii,:),theregressY.umbc_era5.fracWV.sigP(ll,ii,:),theregressY.umbc_era5.fracWV.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(fracWV(ii,ind),nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.gas_1(ii,ind));

    [theregressY.umbc_airsL3.deltaRH.r(ll,ii),theregressY.umbc_airsL3.deltaRH.chisqr(ll,ii),theregressY.umbc_airsL3.deltaRH.P(ll,ii,:),theregressY.umbc_airsL3.deltaRH.sigP(ll,ii,:),theregressY.umbc_airsL3.deltaRH.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(boo1(ind,ii),boo2(ind,ii));
    [theregressY.umbc_cmip6.deltaRH.r(ll,ii),theregressY.umbc_cmip6.deltaRH.chisqr(ll,ii),theregressY.umbc_cmip6.deltaRH.P(ll,ii,:),theregressY.umbc_cmip6.deltaRH.sigP(ll,ii,:),theregressY.umbc_cmip6.deltaRH.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(boo1(ind,ii),boo3(ind,ii));
    [theregressY.umbc_era5.deltaRH.r(ll,ii),theregressY.umbc_era5.deltaRH.chisqr(ll,ii),theregressY.umbc_era5.deltaRH.P(ll,ii,:),theregressY.umbc_era5.deltaRH.sigP(ll,ii,:),theregressY.umbc_era5.deltaRH.numpts(ll,ii,:)] = ...
      nanlinearcorrelation(boo1(ind,ii),boo4(ind,ii));
  end
end
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPlotAllTogether_or_ByModel = input('Enetr (-1) to plot all models together or (+1,default) by model : ');
if length(iPlotAllTogether_or_ByModel) == 0
  iPlotAllTogether_or_ByModel = +1;
end

if iPlotAllTogether_or_ByModel == -1
  figure(9); plot(theregressY.umbc_airsL3.stemp.r,rlat,theregressY.umbc_cmip6.stemp.r,rlat,theregressY.umbc_era5.stemp.r,rlat,'linewidth',2)
    hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('lat'); xlabel('R^2'); title('Correlation dStemp/dt')

%{
  figure(10);
    plot(theregressY.umbc_airsL3.ptemp.r,playsjunk,theregressY.umbc_cmip6.ptemp.r,playsjunk,theregressY.umbc_era5.ptemp.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dT(z)/dt')
  figure(11);
    plot(theregressY.umbc_airsL3.fracWV.r,playsjunk,theregressY.umbc_cmip6.fracWV.r,playsjunk,theregressY.umbc_era5.fracWV.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dfracWV(z)/dt')
  figure(12);
    plot(theregressY.umbc_airsL3.deltaRH.r,playsjunk,theregressY.umbc_cmip6.deltaRH.r,playsjunk,theregressY.umbc_era5.deltaRH.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R^2'); title('Correlation dRH/dt')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
%}
    
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(9); plot(theregressY.umbc_airsL3.stemp.P(:,1),rlat,theregressY.umbc_cmip6.stemp.P(:,1),rlat,theregressY.umbc_era5.stemp.P(:,1),rlat,'linewidth',2)
    hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('lat'); xlabel('S'); title('Slope dStemp/dt')
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  figure(9); plot(theregressY.umbc_airsL3.stemp.r,rlat,theregressY.umbc_cmip6.stemp.r,rlat,theregressY.umbc_era5.stemp.r,rlat,'linewidth',2)
    hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('lat'); xlabel('R'); title('Correlation dStemp/dt')

%{
  figure(10);
    plot(theregressY.umbc_airsL3.ptemp.r,playsjunk,theregressY.umbc_cmip6.ptemp.r,playsjunk,theregressY.umbc_era5.ptemp.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R'); title('Correlation dT(z)/dt')
  figure(11);
    plot(theregressY.umbc_airsL3.fracWV.r,playsjunk,theregressY.umbc_cmip6.fracWV.r,playsjunk,theregressY.umbc_era5.fracWV.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R'); title('Correlation dfracWV(z)/dt')
  figure(12);
    plot(theregressY.umbc_airsL3.deltaRH.r,playsjunk,theregressY.umbc_cmip6.deltaRH.r,playsjunk,theregressY.umbc_era5.deltaRH.r,playsjunk,'linewidth',2)
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); plotaxis2; hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('P(mb)'); xlabel('R'); title('Correlation dRH/dt')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
%}
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt AIRS L3')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt AIRS L3')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause

  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt CMIP6')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt CMIP6')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause

  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.ptemp.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dT(z)/dt ERA5')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.fracWV.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dfracWV(z)/dt ERA5')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.deltaRH.r'); shading interp; caxis([-1.01 +1.01]); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Correlation dRH(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(9); plot(theregressY.umbc_airsL3.stemp.P(:,1),rlat,theregressY.umbc_cmip6.stemp.P(:,1),rlat,theregressY.umbc_era5.stemp.P(:,1),rlat,'linewidth',2)
    hl = legend('AIRS L3','CMIP6','ERA5','location','best','fontsize',8); ylabel('lat'); xlabel('R'); title('Slope dStemp/dt')
  
  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt AIRS L3')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt AIRS L3')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_airsL3.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt AIRS L3')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause

  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt CMIP6')
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt CMIP6')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_cmip6.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt CMIP6')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause

  figure(10);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.ptemp.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dT(z)/dt ERA5')  
  figure(11);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.fracWV.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dfracWV(z)/dt ERA5')
  figure(12);
    pcolor(rlat,playsjunk(1:97),theregressY.umbc_era5.deltaRH.P(:,:,1)'); shading interp; caxis([-1 +1]*2); colormap(llsmap5); colorbar
    set(gca,'yscale','log'); set(gca,'ydir','reverse'); ylabel('P(mb)'); xlabel('lat'); title('Slope dRH(z)/dt ERA5')
  for ii = 10:12;  figure(ii); ylim([10 1000]); end 
  disp('ret to continue'); pause
end
