% get_deltaO3   %%% already done in find_T_RH_trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(25); clf; simplemap(Y(:),X(:),duO3_col0'.*maskLF',5); colorbar; title(['O3 column du']); caxis([250 350]); colormap jet; plotaxis2;
if nlayO3 > 97
  ppmvLAY3         = ppmvLAY3(1:97,:);
  ppmvLAYpert3     = ppmvLAYpert3(1:97,:);
  ppmvLAYpert_unc3 = ppmvLAYpert_unc3(1:97,:);
  nlayO3 = 97;
end
boo = zeros(nlayO3,72,64); for ijunk = 1 : nlayO3; boo(ijunk,:,:) = maskLFmatr'; end
junk0    = boo.*reshape(ppmvLAY3,nlayO3,72,64);
junkpert = boo.*reshape(ppmvLAYpert3,nlayO3,72,64);
pcolor(1:64,1:nlayO3,squeeze(nanmean(junk0,2)));
pcolor(unique(Y(:)),playsjunk(1:nlayO3),squeeze(nanmean(junk0,2))); shading flat; colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('O3 ppmv'); xlabel('Latitude'); ylabel('p(mb)')
  ylim([1 100])
pcolor(unique(Y(:)),p2h(playsjunk(1:nlayO3))/1000,squeeze(nanmean(junk0,2))); shading flat; colorbar; title('O3 ppmv'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); 

figure(26); clf; pcolor(unique(Y(:)),p2h(playsjunk(1:nlayO3))/1000,squeeze(nanmean(junkpert-junk0,2))); shading flat; colorbar; title('\delta O3 ppmv/yr'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); colormap(llsmap5); caxis([-5 +5]*1e-3)

figure(27); clf; pcolor(unique(Y(:)),p2h(playsjunk(1:nlayO3))/1000,100*squeeze(nanmean(junkpert-junk0,2))./squeeze(nanmean(junk0,2))); 
  shading flat; colorbar; title('\delta O3 percent/yr'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); colormap(llsmap5); caxis([-2 +2]*1e-1)

junk = (duO3_colpert - duO3_col0).*maskLF;
figure(25); clf; simplemap(Y(:),X(:),junk'.*maskLF',5); colorbar; title(['d(O3)/dt column du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;
junk = (duO3_300mbpert - duO3_300mb0).*maskLF;
figure(26); clf; simplemap(Y(:),X(:),junk'.*maskLF',5); colorbar; title(['d(O3)/dt TOA->300mb  du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;

%% https://acp.copernicus.org/articles/19/3257/2019/   tropospheric column O3

plot(nanmean(abs(resultsO3)),1:iNumLay,nanmean(resultsO3unc),1:iNumLay); set(gca,'ydir','reverse');

plot(pert.gas_3 ./ p.gas_3 - 1,1:101,'b',pert.gas_3_unc./p.gas_3,1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101])
plot(nanmean(abs(pert.gas_3' ./ p.gas_3' - 1)),1:101,nanmean(pert.gas_3_unc'./p.gas_3'),1:101,'r'); set(gca,'ydir','reverse'); ylim([1 101]); grid; plotaxis2;

fracO3unc = fracO3unc .* (ones(101,1) * maskLF);

figure(25); clf; simplemap(Y(:),X(:),100*fracO3(i050,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 050 mb'])
figure(26); clf; simplemap(Y(:),X(:),100*fracO3(i100,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 100 mb'])
figure(27); clf; simplemap(Y(:),X(:),100*fracO3(i200,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 200 mb'])
for ii = 25 : 27; figure(ii); caxis([-0.25 +0.25]); colormap(llsmap5); plotaxis2; end
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  fracO3lat(ii,:) = nanmean(fracO3(:,boo),2);
end
pcolor(rlat,pavgLAY(1:nlayO3,1000),fracO3lat(:,1:nlayO3)'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(fracO3lat)); junk = cos(rlat) * ones(1,101);
area_wgt_fracO3 = nansum(fracO3lat.*junk,1)./nansum(junk,1);
hold on; plot(area_wgt_fracO3(1:nlayO3)*10000,pavgLAY(1:nlayO3,1000),'color','r','linewidth',2); hold off
ylim([1 100]); caxis([-2 +2]*1e-3); colorbar; plotaxis2;

figure(25); clf; simplemap(Y(:),X(:),deltaO3(i050,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 050 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(26); clf; simplemap(Y(:),X(:),deltaO3(i100,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 100 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(27); clf; simplemap(Y(:),X(:),deltaO3(i200,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 200 mb']); caxis([-0.50 +0.50]/100); plotaxis2;

figure(27); clf; aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(deltaO3(i025,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]/100); colormap(cmap);  title('d/dt UMBC O3(025 mb)'); 
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_o3_025mb_global_trends.pdf');

for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  deltaO3lat(ii,:) = nanmean(deltaO3(:,boo),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(31); clf; aslmap(31,rlat65,rlon73,maskLFmatr.*smoothn(reshape(duO3_F-duO3_0,72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]*1); colormap(cmap);  title('d/dt colO3 (du/year)'); 

%% Global, regional and seasonal analysis of total ozone trends derived from the 1995–2020 GTO-ECV climate data record
%% Melanie Coldewey-Egbers1 , Diego G. Loyola1, Christophe Lerot2 , and Michel van Roozendael2
%% Deutsches Zentrum für Luft- und Raumfahrt (DLR), Institut für Methodik der Fernerkundung, Oberpfaffenhofen, Germ
%% https://doi.org/10.5194/acp-2021-1047
disp('Fig 32 : see Figure 3 in https://acp.copernicus.org/preprints/acp-2021-1047/acp-2021-1047.pdf')
figure(32); clf; aslmap(32,rlat65,rlon73,10*maskLFmatr.*smoothn(reshape(100*(duO3_F-duO3_0)./duO3_0,72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]*2.5); colormap(cmap);  title('d/dt colO3 (percent/decade)'); 
usa22 = usa2; usa22 = usa22(60-40:60+40,:); colormap(usa22);

deltaO3lat = deltaO3lat(:,1:nlayO3);
figure(33); clf; pcolor(rlat,pavgLAY(1:nlayO3,1000),smoothn(deltaO3lat(:,1:nlayO3)',1)); colorbar('horizontal'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaO3lat)); junk = cos(rlat) * ones(1,nlayO3);
area_wgtO3 = nansum(deltaO3lat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtO3(1:nlayO3)*10000,pavgLAY(1:nlayO3,1000),'color','r','linewidth',2); hold off
ylim([0.1 100]); caxis([-0.5 +0.5]*4e-2); title('Zonal d/dt O3 UMBC'); 
colormap(cmap)
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_o3_zonal_trends.pdf');

figure(34); clf; 
aslmap_polar(33,rlat65,rlon73,smoothn((reshape(duO3_F-duO3_0,72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*1); colormap(cmap);  title('d/dt colO3 (du/year)'); 
clear junk
junk.color = 'k'; junk.iNorS = +1; aslmap_polar(33,rlat65,rlon73,smoothn((reshape(duO3_F-duO3_0,72,64)') ,1), [-90 +90],[-180 +180],junk); title('UMBC dcolO3/dt du/yr');  caxis([-1 +1]*1); colormap(llsmap5)
junk.color = 'k'; junk.iNorS = -1; aslmap_polar(33,rlat65,rlon73,smoothn((reshape(duO3_F-duO3_0,72,64)') ,1), [-90 +90],[-180 +180],junk); title('UMBC dcolO3/dt du/yr');  caxis([-1 +1]*1); colormap(llsmap5)

if exist('era5')
  jjunk = era5.trend_gas_3(1:nlayO3,:);
  jjunk = reshape(jjunk,nlayO3,72,64);
  jjunk = squeeze(nanmean(jjunk,2))';
  figure(33); clf; pcolor(rlat,pavgLAY(1:nlayO3,1000),smoothn(jjunk',1)); colorbar('horizontal'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
  junk = zeros(size(deltaO3lat)); junk = cos(rlat) * ones(1,nlayO3);
  area_wgtO3 = nansum(deltaO3lat.*junk,1)./nansum(junk,1);
  %hold on; plot(area_wgtO3(1:nlayO3)*10000,pavgLAY(1:nlayO3,1000),'color','r','linewidth',2); hold off
  ylim([0.1 100]); caxis([-0.5 +0.5]*4e-2); title('Zonal d/dt O3 ERA5'); 
  colormap(cmap)
  %% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era5_o3_zonal_trends.pdf');
end

%{
clear data dataMap
data = smoothn(deltaO3lat(:,1:nlayO3)',1); pnlayO3 = pavgLAY(1:nlayO3,1000); 
dataMap = smoothn(reshape(deltaO3(i025,:),72,64)',1);
save umbc_o3_zonal_trends.mat rlat pnlayO3 data dataMap rlat65 rlon73
%}
