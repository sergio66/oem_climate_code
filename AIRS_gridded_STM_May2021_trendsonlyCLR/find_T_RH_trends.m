load('llsmap5');

h = hMean17years; p = pMean17years;
h = hTimeStep1; p = pTimeStep1;

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

%{
px = p;
px.ptemp = px.ptemp + 0.1;
px.gas_1 = px.gas_1 * (1.01);
[RHx,RH1kmx,colwaterx] = layeramt2RH(h,px);
scatter_coast(px.rlon,px.rlat,50,colwaterx-colwater0); title('\delta colwater mm/yr');

i200 = find(p.plevs(1:90,2300) >= 200,1);
i500 = find(p.plevs(1:90,2300) >= 500,1);
i800 = find(p.plevs(1:90,2300) >= 800,1);
scatter_coast(px.rlon,px.rlat,50,RHx(i800,:)-RH0(i800,:)); title('\delta RH 800 /yr');
scatter_coast(px.rlon,px.rlat,50,RHx(i600,:)-RH0(i600,:)); title('d/dt UMBC RH(600 mb)');
figure(25); aslmap(25,rlat65,rlon73,smoothn(reshape(RHx(i600,:)-RH0(i600,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]); colormap(cmap);  title('d/dt UMBC RH(600 mb)'); 
clear px RHx RH1kmx colwaterx
clear px RHx RH1kmx colwaterx
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pert     = p;
pert_unc = p; 

for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  roo = interp1(log(pavg),resultsT(ii,:),log(playsjunk),[],'extrap');
    pert.ptemp(1:nlays,ii)     = pert.ptemp(1:nlays,ii) + roo;
    pert.ptemp_unc(:,ii)       = 0 * pert.ptemp(:,ii);
    roounc = interp1(log(pavg),resultsTunc(ii,:),log(playsjunk),[],'extrap');  
    pert.ptemp_unc(1:nlays,ii) = roounc;
    pert_unc.ptemp(1:nlays,ii) = pert_unc.ptemp(1:nlays,ii) + roo + roounc;

  roo = interp1(log(pavg),resultsWV(ii,:),log(playsjunk),[],'extrap');
    pert.gas_1(1:nlays,ii)     = pert.gas_1(1:nlays,ii) .* (1+roo);
    pert.gas_1_unc(:,ii)       = 0 * pert.gas_1(:,ii);
    roounc = interp1(log(pavg),resultsWVunc(ii,:),log(playsjunk),[],'extrap');
    pert.gas_1_unc(1:nlays,ii) = pert.gas_1(1:nlays,ii) .* (0+roounc);
    pert_unc.gas_1(1:nlays,ii) = pert_unc.gas_1(1:nlays,ii) .* (1+roo+roounc);

  roo = interp1(log(pavg),resultsO3(ii,:),log(playsjunk),[],'extrap');
    pert.gas_3(1:nlays,ii)     = pert.gas_3(1:nlays,ii) .* (1+roo);
    pert.gas_3_unc(:,ii)       = 0 * pert.gas_3(:,ii);
    roounc = interp1(log(pavg),resultsO3unc(ii,:),log(playsjunk),[],'extrap');
    pert.gas_3_unc(1:nlays,ii) = pert.gas_3(1:nlays,ii) .* (0+roounc);
    pert_unc.gas_3(1:nlays,ii) =  pert_unc.gas_3(1:nlays,ii) .* (1+roo+roounc);

  pert.stemp(ii) = pert.stemp(ii) + results(ii,6);
  pert.stemp_unc(ii) = results(ii,6);
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

i005 = find(p.plevs(1:90,2300) >= 005,1);
i010 = find(p.plevs(1:90,2300) >= 010,1);
i025 = find(p.plevs(1:90,2300) >= 025,1);
i050 = find(p.plevs(1:90,2300) >= 050,1);
i100 = find(p.plevs(1:90,2300) >= 100,1);
i200 = find(p.plevs(1:90,2300) >= 200,1);
i500 = find(p.plevs(1:90,2300) >= 500,1);
i600 = find(p.plevs(1:90,2300) >= 600,1);
i800 = find(p.plevs(1:90,2300) >= 800,1);

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

figure(60); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.a2mgas_1./pert.a2mgas_1_orig-1)); caxis([-0.1 +0.1]);   title('2m NewWV/OrigWV-1'); colormap(usa2);
figure(61); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.a2mtemp-pert.a2mtemp_orig));      caxis([-0.15 +0.15]); title('2m NewT/OrigT'); colormap(usa2);
figure(62); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.gas_1(90,:)./p.gas_1(90,:)-1));   caxis([-0.01 +0.01]); title('815 mb layer 90 NewWV/OrigWV-1'); colormap(usa2);
figure(63); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.ptemp(90,:)-p.ptemp(90,:)));      caxis([-0.15 +0.15]); title('815 mb layer 90 NewT-OrigT'); colormap(usa2);
figure(64); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.gas_1(76,:)./p.gas_1(76,:)-1));   caxis([-0.01 +0.01]); title('500 mb layer 76 NewWV/OrigWV-1'); colormap(usa2);
figure(65); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.ptemp(76,:)-p.ptemp(76,:)));      caxis([-0.15 +0.15]); title('500 mb layer 76 NewT-OrigT'); colormap(usa2);
figure(66); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.gas_1(56,:)./p.gas_1(56,:)-1));   caxis([-0.01 +0.01]); title('200 mb layer 56 NewWV/OrigWV-1'); colormap(usa2);
figure(67); scatter_coast(p.rlon,p.rlat,50,maskLF.*(pert.ptemp(56,:)-p.ptemp(56,:)));      caxis([-0.15 +0.15]); title('200 mb layer 56 NewT-OrigT'); colormap(usa2);

figure(60); junk = pert.a2mgas_1./pert.a2mgas_1_orig-1; aslmap(60,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.10 +0.10]);   title('2m NewWV/OrigWV-1');     colormap(usa2);
figure(61); junk = pert.a2mtemp-pert.a2mtemp_orig;      aslmap(61,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);   title('2m NewT-OrigT');         colormap(usa2);
figure(62); junk = pert.gas_1(90,:)./p.gas_1(90,:)-1;   aslmap(62,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.005 +0.005]); title('815 mb NewWV/OrigWV-1'); colormap(usa2);
figure(63); junk = pert.ptemp(90,:)-p.ptemp(90,:);      aslmap(63,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);   title('815 mb NewT-OrigT');     colormap(usa2);
figure(64); junk = pert.gas_1(76,:)./p.gas_1(76,:)-1;   aslmap(64,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.005 +0.005]); title('500 mb NewWV/OrigWV-1'); colormap(usa2);
figure(65); junk = pert.ptemp(76,:)-p.ptemp(76,:);      aslmap(65,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);   title('500 mb NewT-OrigT');     colormap(usa2);
figure(66); junk = pert.gas_1(56,:)./p.gas_1(56,:)-1;   aslmap(66,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.005 +0.005]); title('200 mb NewWV/OrigWV-1'); colormap(usa2);
figure(67); junk = pert.ptemp(56,:)-p.ptemp(56,:);      aslmap(67,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);   title('200 mb NewT-OrigT');     colormap(usa2);

pert.RH      = RHpert;
pert.RH_orig = RH0;
figure(60); junk = pert.a2mRH-pert.a2mRH_orig;          aslmap(60,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.05 +0.05]*10);   title('2m NewRH-OrigRH');     colormap(llsmap5);
figure(61); junk = pert.a2mtemp-pert.a2mtemp_orig;      aslmap(61,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);      title('2m NewT-OrigT');       colormap(llsmap5);
figure(62); junk = pert.RH(90,:)-pert.RH_orig(90,:);    aslmap(62,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.05 +0.05]*10);   title('815 mb NewRH-OrigRH'); colormap(llsmap5);
figure(63); junk = pert.ptemp(90,:)-p.ptemp(90,:);      aslmap(63,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);      title('815 mb NewT-OrigT');   colormap(llsmap5);
figure(64); junk = pert.RH(76,:)-pert.RH_orig(76,:);    aslmap(64,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.05 +0.05]*10);   title('500 mb NewRH-OrigRH'); colormap(llsmap5);
figure(65); junk = pert.ptemp(76,:)-p.ptemp(76,:);      aslmap(65,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);      title('500 mb NewT-OrigT');   colormap(llsmap5);
figure(66); junk = pert.RH(56,:)-pert.RH_orig(56,:);    aslmap(66,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.05 +0.05]*10);   title('200 mb NewRH-OrigRH'); colormap(llsmap5);
figure(67); junk = pert.ptemp(56,:)-p.ptemp(56,:);      aslmap(67,rlat65,rlon73,maskLFmatr.*smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); caxis([-0.15 +0.15]);      title('200 mb NewT-OrigT');   colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(nanmean(abs(resultsWV)),1:iNumLay,nanmean(resultsWVunc),1:iNumLay); set(gca,'ydir','reverse');
fracWV = pert.gas_1 ./ p.gas_1 - 1;
fracWV = fracWV .* (ones(101,1)*maskLF);

fracWVunc = pert_unc.gas_1 ./ p.gas_1 - 1;
fracWVunc = pert.gas_1_unc ./ p.gas_1;
plot(pert.gas_1 ./ p.gas_1 - 1,1:101,'b',pert.gas_1_unc./p.gas_1,1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101])
plot(nanmean(abs(pert.gas_1' ./ p.gas_1' - 1)),1:101,nanmean(pert.gas_1_unc'./p.gas_1'),1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101]); grid; plotaxis2;
plot(nanmean(abs(pert.gas_1' ./ p.gas_1' - 1)),pert.plevs(:,2345),nanmean(pert.gas_1_unc'./p.gas_1'),pert.plevs(:,2345),'r'); set(gca,'ydir','reverse'); ylim([100 1000]); grid; plotaxis2;

fracWVunc = fracWVunc .* (ones(101,1)*maskLF);

figure(25); simplemap(Y(:),X(:),100*fracWV(i050,:)'.*maskLF',5); colorbar; title(['percent d(fracWV)/dt yr-1 050 mb'])
figure(26); simplemap(Y(:),X(:),100*fracWV(i100,:)'.*maskLF',5); colorbar; title(['percent d(fracWV)/dt yr-1 100 mb'])
figure(27); simplemap(Y(:),X(:),100*fracWV(i200,:)'.*maskLF',5); colorbar; title(['percent d(fracWV)/dt yr-1 200 mb'])
for ii = 25 : 27; figure(ii); caxis([-0.25 +0.25]); colormap(llsmap5); plotaxis2; end
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  fracWVlat(ii,:) = nanmean(fracWV(:,boo),2);
end

cmap = llsmap5;

figure(30); 
pcolor(rlat,pavgLAY(1:97,1000),fracWVlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracWVlat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracWV = nansum(fracWVlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgt_fracWV(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-2 +2]*1e-3); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
colormap(cmap)

%%%%%%%%%%%%%%%%%%%%%%%%%

deltaRH    = (RHpert-RH0) .* (ones(100,1)*maskLF);;     %% DO NOT USE deltaRH(97,:) as surface change in RH, use  RHSurfpert-RHSurf0
deltaRHunc = (RHpert_unc-RH0) .* (ones(100,1)*maskLF);; %% DO NOT USE deltaRH(97,:) as surface change in RH, use  RHSurfpert-RHSurf0

figure(41); clf; plot(mean(deltaRH'),1:100,mean(deltaRHunc'),1:100,'b--'); plotaxis2; set(gca,'ydir','reverse'); ylim([0 100])

figure(25); clf; simplemap(Y(:),X(:),deltaRH(i200,:)'.*maskLF',5); colorbar; title(['d(RH)/dt yr-1 200 mb'])
figure(26); clf; simplemap(Y(:),X(:),deltaRH(i500,:)'.*maskLF',5); colorbar; title(['d(RH)/dt yr-1 500 mb'])
figure(27); clf; simplemap(Y(:),X(:),deltaRH(i800,:)'.*maskLF',5); colorbar; title(['d(RH)/dt yr-1 800 mb'])
figure(27); clf; simplemap(Y(:),X(:),deltaRH(i600,:)'.*maskLF',5); colorbar; title(['d(RH)/dt yr-1 600 mb'])
for ii = 25 : 27; figure(ii); caxis([-1 +1]); colormap(llsmap5); plotaxis2; end
for ii = 25 : 27; figure(ii); caxis([-0.25 +0.25]); colormap(llsmap5); plotaxis2; end

iUT = find(pavgLAY(1:97,3000) >= 200 & pavgLAY(1:97,3000) <= 500);
% Trends in Upper-Tropospheric Humidity: Expansion of the Subtropical Dry Zones? DOI: 10.1175/JCLI-D-19-0046.1
% MIRIAM TIVIG AND VERENA GRUTZUN. J. Clim 2020
boo = nanmean(deltaRH(iUT,:));
figure(27); aslmap(27,rlat65,rlon73,10*maskLFmatr.*smoothn(reshape(boo,72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]/2); caxis([-10 +10]/4); colormap(cmap);  title('d/dt UMBC RH(UT 200-500 mb) %/decade'); 
addpath /home/sergio/MATLABCODE/COLORMAP/COLORBREWER/cbrewer/cbrewer
ct = cbrewer('div','BrBG',8); colormap(ct);
ct = cbrewer('div','BrBG',63,'spline');ct(ct < 0) = 0; ct(ct > 1) = 1;  colormap(ct);

%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_rh_ut_500_200mb_global_trends.pdf');

figure(27); aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(deltaRH(i600,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]); colormap(cmap);  title('d/dt UMBC RH(600 mb)'); 
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_rh_600mb_global_trends.pdf');

figure(27); aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(deltaRH(i500,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]); colormap(cmap);  title('d/dt UMBC RH(500 mb)'); 
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_rh_500mb_global_trends.pdf');

for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  deltaRHlat(ii,:) = nanmean(deltaRH(:,boo),2);
end
figure(28);
pcolor(rlat,pavgLAY(1:97,3000),deltaRHlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaRHlat)); junk = cos(rlat) * ones(1,100);
area_wgtRH = nansum(deltaRHlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtRH(1:97)*100,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-0.15 +0.15]); colorbar('horizontal'); colormap(cmap); title(['Zonal d/dt RH UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_rh_zonal_trends.pdf');

%{
clear data dataMap
data = maskLFmatr.*smoothn(deltaRHlat(:,1:97)',1); p97 = pavgLAY(1:97,1000); 
dataMap = maskLFmatr.*smoothn(reshape(deltaRH(i600,:),72,64)',1);
save umbc_RH_zonal_trends.mat rlat p97 data dataMap rlat65 rlon73
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(nanmean(abs(resultsT)),1:iNumLay,nanmean(resultsTunc),1:iNumLay); set(gca,'ydir','reverse');
deltaT = pert.ptemp-p.ptemp;
deltaT = deltaT .* (ones(101,1) * maskLF);

deltaTunc = pert_unc.ptemp-p.ptemp;
deltaTunc = pert.ptemp_unc;
plot(pert.ptemp-p.ptemp,1:101,'b',pert.ptemp_unc,1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101])
plot(nanmean(abs(pert.ptemp'-p.ptemp')),1:101,nanmean(pert.ptemp_unc'),1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101]); grid; plotaxis2;
semilogy(nanmean(abs(pert.ptemp'-p.ptemp')),p.plevs(:,2345),nanmean(pert.ptemp_unc'),p.plevs(:,2345),'r'); set(gca,'ydir','reverse'); ylim([10 1000]); grid; plotaxis2;

deltaTunc = deltaTunc .* (ones(101,1) * maskLF);
figure(41); clf; plot(nanmean(deltaT'),1:101,nanmean(deltaTunc'),1:101,'b--'); plotaxis2; set(gca,'ydir','reverse'); ylim([0 100])

figure(25); simplemap(Y(:),X(:),deltaT(i200,:)'.*maskLF',5); colorbar; title(['d(T)/dt yr-1 200 mb'])
figure(26); simplemap(Y(:),X(:),deltaT(i500,:)'.*maskLF',5); colorbar; title(['d(T)/dt yr-1 500 mb'])
figure(27); simplemap(Y(:),X(:),deltaT(i800,:)'.*maskLF',5); colorbar; title(['d(T)/dt yr-1 800 mb'])
%% figure(27); simplemap(Y(:),X(:),deltaT(i600,:)'.*maskLF',5); colorbar; title(['d(T)/dt yr-1 600 mb'])
for ii = 25 : 27; figure(ii); caxis([-1 +1]); colormap(llsmap5); plotaxis2; end
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  deltaTlat(ii,:) = nanmean(deltaT(:,boo),2);
end

figure(27); clf; aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(deltaT(i500,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]/10); colormap(cmap);  title('d/dt UMBC T(600 mb)'); 
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_T_600mb_global_trends.pdf');

figure(29); clf
pcolor(rlat,pavgLAY(1:97,3000),deltaTlat(:,1:97)'); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(deltaTlat(:,1:97)',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaTlat)); junk = cos(rlat) * ones(1,101);
%area_wgtT = nansum(deltaTlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtT(1:97)*1000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([10 1000]); caxis([-0.15 +0.15]); colorbar('horizontal'); colormap(cmap); title(['Zonal d/dt T UMBC Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_T_zonal_trends.pdf');

%{
clear data dataMap
data = maskLFmatr.*smoothn(deltaTlat(:,1:97)',1); p97 = pavgLAY(1:97,1000); 
dataMap = maskLFmatr.*smoothn(reshape(deltaT(i500,:),72,64)',1);
save umbc_T_zonal_trends.mat rlat p97 data dataMap rlat65 rlon73
%}

%%%%%%%%%%%%%%%%%%%%%%%%%
