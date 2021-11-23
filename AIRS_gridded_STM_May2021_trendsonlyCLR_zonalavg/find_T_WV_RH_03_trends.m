%% see make_avgprof64x72_to_64.m
[havg64,~,pavg64,~] = rtpread('pavg64.op.rtp');

load('llsmap5');

h = havg64; p = pavg64;

pjunkN = p.plevs(1:100,:)-p.plevs(2:101,:);
pjunkD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
pavgLAY = pjunkN./pjunkD;

[xRH0,xRH1km0,xcolwater0] = layeramt2RH(h,p);
[Tw0,Tw1km0,Tdew0,WBGT0,RH0,RH1km0,colwater0,TwSurf0,RHSurf0,TdewSurf0] = layeramt2RH_wet_bulb_dew_point_temperature(h,p);

if ~isfield(p,'plays')
  playsN = p.plevs(1:100,:)-p.plevs(2:101,:);
  playsD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
  p.plays = zeros(size(p.plevs));
  p.plays(1:100,:) = playsN./playsD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pert     = p;
pert_unc = p; 

for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  roo = interp1(log(pavg),resultsT(ii,:),log(playsjunk),[],'extrap');
  pert.ptemp(1:nlays,ii) =  pert.ptemp(1:nlays,ii) + roo;
  roounc = interp1(log(pavg),resultsTunc(ii,:),log(playsjunk),[],'extrap');
  pert_unc.ptemp(1:nlays,ii) =  pert_unc.ptemp(1:nlays,ii) + roo + roounc;

  roo = interp1(log(pavg),resultsWV(ii,:),log(playsjunk),[],'extrap');
  pert.gas_1(1:nlays,ii) =  pert.gas_1(1:nlays,ii) .* (1+roo);
  roounc = interp1(log(pavg),resultsWVunc(ii,:),log(playsjunk),[],'extrap');
  pert_unc.gas_1(1:nlays,ii) =  pert_unc.gas_1(1:nlays,ii) .* (1+roo+roounc);

  roo = interp1(log(pavg),resultsO3(ii,:),log(playsjunk),[],'extrap');
  pert.gas_3(1:nlays,ii) =  pert.gas_3(1:nlays,ii) .* (1+roo);
  roounc = interp1(log(pavg),resultsO3unc(ii,:),log(playsjunk),[],'extrap');
  pert_unc.gas_3(1:nlays,ii) =  pert_unc.gas_3(1:nlays,ii) .* (1+roo+roounc);

  pert.stemp(ii) = pert.stemp(ii) + results(ii,6);
  pert_unc.stemp(ii) = pert_unc.stemp(ii) + resultsunc(ii,6);

  pert.gas_2(1:nlays,ii) =  pert.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert.gas_4(1:nlays,ii) =  pert.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert.gas_6(1:nlays,ii) =  pert.gas_6(1:nlays,ii) .* (1+4.5/1700);
  pert_unc.gas_2(1:nlays,ii) =  pert_unc.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert_unc.gas_4(1:nlays,ii) =  pert_unc.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert_unc.gas_6(1:nlays,ii) =  pert_unc.gas_6(1:nlays,ii) .* (1+4.5/1700);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rtpwrite('summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp',h,ha,pert,pa);

[xRHpert,xRH1kmpert,xcolwaterpert] = layeramt2RH(h,pert);
[Twpert,Tw1kmpert,Tdewpert,WBGTpert,RHpert,RH1kmpert,colwaterpert,TwSurfpert,RHSurfpert,TdewSurfpert] = layeramt2RH_wet_bulb_dew_point_temperature(h,pert);

[xRHpert_unc,xRH1kmpert_unc,xcolwaterpert_unc] = layeramt2RH(h,pert_unc);
[Twpert_unc,Tw1kmpert_unc,Tdewpert_unc,WBGTpert_unc,RHpert_unc,RH1kmpert_unc,colwaterpert_unc,TwSurfpert_unc,RHSurfpert_unc,TdewSurfpert_unc] = layeramt2RH_wet_bulb_dew_point_temperature(h,pert_unc);

scatter_coast(pert.rlon,pert.rlat,50,colwaterpert-colwater0); title('\delta colwater mm/yr');

i005 = find(p.plevs(1:90,32) >= 005,1);
i010 = find(p.plevs(1:90,32) >= 010,1);
i025 = find(p.plevs(1:90,32) >= 025,1);
i050 = find(p.plevs(1:90,32) >= 050,1);
i100 = find(p.plevs(1:90,32) >= 100,1);
i200 = find(p.plevs(1:90,32) >= 200,1);
i500 = find(p.plevs(1:90,32) >= 500,1);
i600 = find(p.plevs(1:90,32) >= 600,1);
i800 = find(p.plevs(1:90,32) >= 800,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  ptempjunk = p.ptemp(1:nlays,ii);
  roo0 = interp1(log(playsjunk),ptempjunk,log(p.spres(ii)),[],'extrap');
  roo  = interp1(log(pavg),resultsT(ii,:),log(p.spres(ii)),[],'extrap');
  pert.a2mtemp_orig(ii) =  roo0;
  pert.a2mtemp(ii) =  roo0 + roo;

  gas1junk = p.gas_1(1:nlays,ii);
  roo0 = interp1(log(playsjunk),log(gas1junk),log(p.spres(ii)),[],'extrap');
  roo0 = exp(roo0);
  pert.a2mgas_1_orig(ii) =  roo0;
  roo  = interp1(log(pavg),resultsT(ii,:),log(p.spres(ii)),[],'extrap');
  pert.a2mgas_1(ii) =  roo0 * (1+roo);

  rhjunk = RH0(1:nlays,ii);
  roo0 = interp1(log(playsjunk),rhjunk,log(p.spres(ii)),[],'extrap');
  rhjunk = RHpert(1:nlays,ii);
  roo = interp1(log(playsjunk),rhjunk,log(p.spres(ii)),[],'extrap');
  pert.a2mRH_orig(ii) =  roo0;
  pert.a2mRH(ii)      =  roo;
end

maskLF = ones(1,64);

figure(60); plot(p.rlat,maskLF.*(pert.a2mgas_1./pert.a2mgas_1_orig-1)); caxis([-0.1 +0.1]);   title('2m NewWV/OrigWV-1'); colormap(usa2);
figure(61); plot(p.rlat,maskLF.*(pert.a2mtemp-pert.a2mtemp_orig));      caxis([-0.15 +0.15]); title('2m NewT/OrigT'); colormap(usa2);
figure(62); plot(p.rlat,maskLF.*(pert.gas_1(90,:)./p.gas_1(90,:)-1));   caxis([-0.01 +0.01]); title('815 mb layer 90 NewWV/OrigWV-1'); colormap(usa2);
figure(63); plot(p.rlat,maskLF.*(pert.ptemp(90,:)-p.ptemp(90,:)));      caxis([-0.15 +0.15]); title('815 mb layer 90 NewT-OrigT'); colormap(usa2);
figure(64); plot(p.rlat,maskLF.*(pert.gas_1(76,:)./p.gas_1(76,:)-1));   caxis([-0.01 +0.01]); title('500 mb layer 76 NewWV/OrigWV-1'); colormap(usa2);
figure(65); plot(p.rlat,maskLF.*(pert.ptemp(76,:)-p.ptemp(76,:)));      caxis([-0.15 +0.15]); title('500 mb layer 76 NewT-OrigT'); colormap(usa2);
figure(66); plot(p.rlat,maskLF.*(pert.gas_1(56,:)./p.gas_1(56,:)-1));   caxis([-0.01 +0.01]); title('200 mb layer 56 NewWV/OrigWV-1'); colormap(usa2);
figure(67); plot(p.rlat,maskLF.*(pert.ptemp(56,:)-p.ptemp(56,:)));      caxis([-0.15 +0.15]); title('200 mb layer 56 NewT-OrigT'); colormap(usa2);

pert.RH      = RHpert;
pert.RH_orig = RH0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fracWV = pert.gas_1 ./ p.gas_1 - 1;
fracWV = fracWV .* (ones(101,1)*maskLF);

fracWVunc = pert_unc.gas_1 ./ p.gas_1 - 1;
fracWVunc = fracWVunc .* (ones(101,1)*maskLF);

fracWVlat = fracWV';
fracWVlatunc = fracWVunc';

cmap = llsmap5;

figure(2); clf 
pcolor(rlat,pavgLAY(1:97,32),fracWVlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracWVlat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracWV = nansum(fracWVlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgt_fracWV(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-2 +2]*1e-3); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
colormap(cmap)

figure(41); clf; plot(mean(fracWV'),1:101,mean(fracWVunc'),1:101,'b--'); ; plotaxis2; set(gca,'ydir','reverse'); ylim([0 100])

%%%%%%%%%%%%%%%%%%%%%%%%%

deltaRH    = (RHpert-RH0) .* (ones(100,1)*maskLF);;     %% DO NOT USE deltaRH(97,:) as surface change in RH, use  RHSurfpert-RHSurf0
deltaRHunc = (RHpert_unc-RH0) .* (ones(100,1)*maskLF);; %% DO NOT USE deltaRH(97,:) as surface change in RH, use  RHSurfpert-RHSurf0

deltaRHlat = deltaRH'; 
deltaRHlatunc = deltaRHunc'; 

figure(3); clf 
pcolor(rlat,pavgLAY(1:97,32),deltaRHlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracWVlat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracWV = nansum(fracWVlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgt_fracWV(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-25 +25]*1e-2); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt RH UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
colormap(cmap)


figure(41); clf; plot(mean(deltaRH'),1:100,mean(deltaRHunc'),1:100,'b--'); plotaxis2; set(gca,'ydir','reverse'); ylim([0 100])

%%%%%%%%%%%%%%%%%%%%%%%%%

deltaT = pert.ptemp-p.ptemp;
deltaT = deltaT .* (ones(101,1) * maskLF);

deltaTunc = pert_unc.ptemp-p.ptemp;
deltaTunc = deltaTunc .* (ones(101,1) * maskLF);

figure(41); clf; plot(nanmean(deltaT'),1:101,nanmean(deltaTunc'),1:101,'b--'); plotaxis2; set(gca,'ydir','reverse'); ylim([0 100])

deltaTlat = deltaT';
deltaTlatunc = deltaTunc';

figure(1); clf
pcolor(rlat,pavgLAY(1:97,32),deltaTlat(:,1:97)'); 
pcolor(rlat,pavgLAY(1:97,32),smoothn(deltaTlat(:,1:97)',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaTlat)); junk = cos(rlat) * ones(1,101);
%area_wgtT = nansum(deltaTlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtT(1:97)*1000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-0.15 +0.15]); colorbar('horizontal'); colormap(cmap); title(['Zonal d/dt T UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_T_zonal_trends.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%

fracO3 = pert.gas_3 ./ p.gas_3 - 1;
fracO3 = fracO3 .* (ones(101,1)*maskLF);

fracO3unc = pert_unc.gas_3 ./ p.gas_3 - 1;
fracO3unc = fracO3unc .* (ones(101,1)*maskLF);

fracO3lat = fracO3';
fracO3latunc = fracO3unc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); ylim([10 1000]);  caxis([-0.15 +0.15]);
figure(2); ylim([100 1000]); caxis([-0.15 +0.15]/10);
figure(3); ylim([100 1000]); caxis([-0.25 +0.25]);

