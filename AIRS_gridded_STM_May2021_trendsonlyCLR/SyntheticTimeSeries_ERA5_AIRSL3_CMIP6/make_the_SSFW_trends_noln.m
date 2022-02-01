figure(1); clf
figure(2); clf
figure(3); clf

figure(1); 
plot(nwp_trends.era5_100_layertrends.gas_1(ixN,ind),ixN,'o-',quickgas_1rate(ixN),ixN); set(gca,'ydir','reverse')
semilogy(nwp_trends.era5_100_layertrends.gas_1(ixN,ind),p72avg.plevs(ixN),'o-',quickgas_1rate(ixN),p72avg.plevs(ixN)); set(gca,'ydir','reverse'); ylim([0.1 1000]); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ixN = 1 : p72avg.nlevs-1;
raReconstructWithUnc = zeros(1,2645);

if iRampCO2_CH4_N2O > 0
  raReconstructWithUnc = raReconstructWithUnc + jCO2';
  raReconstructWithUnc = raReconstructWithUnc + jN2O';
  raReconstructWithUnc = raReconstructWithUnc + jCH4';
end

raReconstructWithUnc = raReconstructWithUnc + jST'*(quickstemprate+errstemprate);
bah = ones(2645,1) * ((quickgas_1rate(ixN)'+errgas_1rate(ixN)')-0); raReconstructWithUnc = raReconstructWithUnc + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * ((quickgas_3rate(ixN)'+errgas_3rate(ixN)')-0); raReconstructWithUnc = raReconstructWithUnc + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * (quickptemprate(ixN)'+errptemprate(ixN)'); raReconstructWithUnc = raReconstructWithUnc + sum(bah.*jT(:,ixN),2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ixN = 1 : p72avg.nlevs-1;
raReconstruct = zeros(1,2645);

if iRampCO2_CH4_N2O > 0
  raReconstruct = raReconstruct + jCO2';
  raReconstruct = raReconstruct + jN2O';
  raReconstruct = raReconstruct + jCH4';
end

%{
%% orig
raReconstruct = raReconstruct + jST'*nwp_trends.era5_100_layertrends.stemp(ind);
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_1(ixN,ind)'); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_3(ixN,ind)'); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ixN,ind)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';
%}

%{
%% very similar if rates are tiny
raReconstruct = raReconstruct + jST'*nwp_trends.era5_100_layertrends.stemp(ind);
bah = ones(2645,1) * ((nwp_trends.era5_100_layertrends.gas_1(ixN,ind)')-0); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * ((nwp_trends.era5_100_layertrends.gas_3(ixN,ind)')-0); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ixN,ind)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';
%}

figure(2); plot(errgas_1rate(ixN),ixN,'k--',quickgas_1rate(ixN)',ixN,'bo-')
  set(gca,'ydir','reverse')

raReconstruct = raReconstruct + jST'*quickstemprate;
bah = ones(2645,1) * ((quickgas_1rate(ixN)')-0); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * ((quickgas_3rate(ixN)')-0); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * quickptemprate(ixN)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c'); 
   hl = legend('SARTA NWP trend','quicktrend','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-'); 
   hl = legend('Reconstructed jac x dX/dt','NEW RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c',fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-',fKc,obsx,'k.-'); 
   hl = legend('SARTA NWP trend','quicktrend','Reconstructed jac x dX/dt','NEW RECONSTRUCT','ACTUAL AIRS OBS','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c',fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-'); 
   hl = legend('SARTA NWP trend','quicktrend','Reconstructed jac x dX/dt','NEW RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,raReconstruct,'r-'); 
   plotaxis2;
   hl = legend('SARTA NWP trend','RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])
title('Using all 19yrs x 12 months')
%}

figure(3)
plot(quickptemprate(ixN)',ixN); set(gca,'ydir','reverse'); plotaxis2;
plot(((quickgas_1rate(ixN)')-0),ixN); set(gca,'ydir','reverse'); plotaxis2;

plot(fKc,quickBTtrend,'b.-',fKc,raReconstruct,'r-',fKc,raReconstructWithUnc,'g-');
hold on; plot(fKc,quickBTtrend+errBTtrend,'color',[1 1 1]*0.7)
hold on; plot(fKc,quickBTtrend-errBTtrend,'color',[1 1 1]*0.7)
hold off
   plotaxis2;
   hl = legend('SARTA QUICKTREND','RECONSTRUCT','RECONTRUCT WITH UNC','location','best','fontsize',8);
   xlim([640 1640])
title(['Using all 19yrs Season ' num2str(iSeason)])

sartatrendX = sartatrend(:,iLonBin);
saver = ['save reconstruct_SSFW_season' num2str(iSeason) '.mat fKc quickBTtrend raReconstruct raReconstructWithUnc errBTtrend sartatrendX'];
eval(saver);
