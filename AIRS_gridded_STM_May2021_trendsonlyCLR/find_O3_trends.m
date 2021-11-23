duO3_col0   = dobson_gas_rtp(h, p, 3);
duO3_300mb0 = dobson_gas_rtp(h, p, 3, 300);
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,p,1:length(p.stemp),3);

duO3_colpert   = dobson_gas_rtp(h, pert, 3);
duO3_300mbpert = dobson_gas_rtp(h, pert, 3, 300);
[ppmvLAYpert,ppmvAVGpert,ppmvMAXpert,pavgLAYpert,tavgLAYpert,ppmv500pert,ppmv75pert,ppmvSURFpert] = layers2ppmv(h,pert,1:length(p.stemp),3);
[ppmvLAYpert_unc,ppmvAVGpert_unc,ppmvMAXpert_unc,pavgLAYpert_unc,tavgLAYpert_unc,ppmv500pert_unc,ppmv75pert_unc,ppmvSURFpert_unc] = layers2ppmv(h,pert_unc,1:length(p.stemp),3);

figure(25); simplemap(Y(:),X(:),duO3_col0'.*maskLF',5); colorbar; title(['O3 column du']); caxis([250 350]); colormap jet; plotaxis2;
[nlayO3,~] = size(ppmvLAY);
if nlayO3 > 97
  ppmvLAY     = ppmvLAY(1:97,:);
  ppmvLAYpert = ppmvLAYpert(1:97,:);
  ppmvLAYpert_unc = ppmvLAYpert_unc(1:97,:);
  nlayO3 = 97;
end
boo = zeros(nlayO3,72,64); for ijunk = 1 : nlayO3; boo(ijunk,:,:) = maskLFmatr'; end
junk0    = boo.*reshape(ppmvLAY,nlayO3,72,64);
junkpert = boo.*reshape(ppmvLAYpert,nlayO3,72,64);
pcolor(1:64,1:nlayO3,squeeze(nanmean(junk0,2)));
pcolor(unique(Y(:)),playsjunk,squeeze(nanmean(junk0,2))); shading flat; colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('O3 ppmv'); xlabel('Latitude'); ylabel('p(mb)')
  ylim([1 100])
pcolor(unique(Y(:)),p2h(playsjunk)/1000,squeeze(nanmean(junk0,2))); shading flat; colorbar; title('O3 ppmv'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); 

figure(26); pcolor(unique(Y(:)),p2h(playsjunk)/1000,squeeze(nanmean(junkpert-junk0,2))); shading flat; colorbar; title('\delta O3 ppmv/yr'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); colormap(llsmap5); caxis([-5 +5]*1e-3)

figure(27); pcolor(unique(Y(:)),p2h(playsjunk)/1000,100*squeeze(nanmean(junkpert-junk0,2))./squeeze(nanmean(junk0,2))); 
  shading flat; colorbar; title('\delta O3 percent/yr'); xlabel('Latitude'); ylabel('h(km)')
  ylim([16 48]); ylim([10 60]); colormap(llsmap5); caxis([-2 +2]*1e-1)

deltaO3 = (duO3_colpert - duO3_col0).*maskLF;
figure(25); simplemap(Y(:),X(:),deltaO3'.*maskLF',5); colorbar; title(['d(O3)/dt column du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;
deltaO3 = (duO3_300mbpert - duO3_300mb0).*maskLF;
figure(26); simplemap(Y(:),X(:),deltaO3'.*maskLF',5); colorbar; title(['d(O3)/dt TOA->300mb  du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;

%% https://acp.copernicus.org/articles/19/3257/2019/   tropospheric column O3

fracO3 = pert.gas_3 ./ p.gas_3 - 1;
fracO3 = fracO3 .* (ones(101,1) * maskLF);

fracO3unc = pert_unc.gas_3 ./ p.gas_3 - 1;
fracO3unc = fracO3unc .* (ones(101,1) * maskLF);

figure(25); simplemap(Y(:),X(:),100*fracO3(i050,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 050 mb'])
figure(26); simplemap(Y(:),X(:),100*fracO3(i100,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 100 mb'])
figure(27); simplemap(Y(:),X(:),100*fracO3(i200,:)'.*maskLF',5); colorbar; title(['percent d(fracO3)/dt yr-1 200 mb'])
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

deltaO3 = ppmvLAYpert - ppmvLAY;
deltaO3 = deltaO3 .* (ones(nlayO3,1) * maskLF);

deltaO3unc = ppmvLAYpert_unc - ppmvLAY;
deltaO3unc = deltaO3unc .* (ones(nlayO3,1) * maskLF);

figure(25); simplemap(Y(:),X(:),deltaO3(i050,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 050 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(26); simplemap(Y(:),X(:),deltaO3(i100,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 100 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(27); simplemap(Y(:),X(:),deltaO3(i200,:)'.*maskLF',5); colorbar; title(['d(O3)/dt ppm/yr 200 mb']); caxis([-0.50 +0.50]/100); plotaxis2;

figure(27); aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(deltaO3(i025,:),72,64)',1), [-90 +90],[-180 +180]);
caxis([-1 +1]/100); colormap(cmap);  title('d/dt UMBC O3(025 mb)'); 
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_o3_025mb_global_trends.pdf');

for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  deltaO3lat(ii,:) = nanmean(deltaO3(:,boo),2);
end

figure(31); pcolor(rlat,pavgLAY(1:nlayO3,1000),smoothn(deltaO3lat(:,1:nlayO3)',1)); colorbar('horizontal'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(deltaO3lat)); junk = cos(rlat) * ones(1,nlayO3);
area_wgtO3 = nansum(deltaO3lat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtO3(1:nlayO3)*10000,pavgLAY(1:nlayO3,1000),'color','r','linewidth',2); hold off
ylim([0.1 100]); caxis([-0.5 +0.5]*4e-2); title('Zonal d/dt O3 UMBC'); %colorbar; plotaxis2; title('d/dt UMBC O3(025 mb)'); 
colormap(cmap)
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_o3_zonal_trends.pdf');

%{
clear data dataMap
data = smoothn(deltaO3lat(:,1:nlayO3)',1); pnlayO3 = pavgLAY(1:nlayO3,1000); 
dataMap = smoothn(reshape(deltaO3(i025,:),72,64)',1);
save umbc_o3_zonal_trends.mat rlat pnlayO3 data dataMap rlat65 rlon73
%}
