[RH0,RH1km0,colwater0] = layeramt2RH(h,p);
duO3_col0   = dobson_gas_rtp(h, p, 3);
duO3_300mb0 = dobson_gas_rtp(h, p, 3, 300);
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,p,1:length(p.stemp),3);

px = p;
px.ptemp = px.ptemp + 0.1;
px.gas_1 = px.gas_1 * (1.01);
[RHx,RH1kmx,colwaterx] = layeramt2RH(h,px);
scatter_coast(px.rlon,px.rlat,50,colwaterx-colwater0); title('\delta colwater mm/yr');

i200 = find(p.plevs(1:90,2300) >= 200,1);
i500 = find(p.plevs(1:90,2300) >= 500,1);
i800 = find(p.plevs(1:90,2300) >= 800,1);
scatter_coast(px.rlon,px.rlat,50,RHx(i800,:)-RH0(i800,:)); title('\delta RH 800 /yr');

if ~isfield(p,'plays')
  playsN = p.plevs(1:100,:)-p.plevs(2:101,:);
  playsD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
  p.plays = zeros(size(p.plevs));
  p.plays(1:100,:) = playsN./playsD;
end

pert = p;
for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  roo = interp1(log(pavg),resultsT(ii,:),log(playsjunk),[],'extrap');
  pert.ptemp(1:nlays,ii) =  pert.ptemp(1:nlays,ii) + roo;

  roo = interp1(log(pavg),resultsWV(ii,:),log(playsjunk),[],'extrap');
  pert.gas_1(1:nlays,ii) =  pert.gas_1(1:nlays,ii) .* (1+roo);

  roo = interp1(log(pavg),resultsO3(ii,:),log(playsjunk),[],'extrap');
  pert.gas_3(1:nlays,ii) =  pert.gas_3(1:nlays,ii) .* (1+roo);

  pert.stemp(ii) = pert.stemp(ii) + results(ii,6);

  pert.gas_2(1:nlays,ii) =  pert.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert.gas_4(1:nlays,ii) =  pert.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert.gas_6(1:nlays,ii) =  pert.gas_6(1:nlays,ii) .* (1+4.5/1700);
end
% rtpwrite('summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp',h,ha,pert,pa);

[RHpert,RH1kmpert,colwaterpert] = layeramt2RH(h,pert);
duO3_colpert   = dobson_gas_rtp(h, pert, 3);
duO3_300mbpert = dobson_gas_rtp(h, pert, 3, 300);
[ppmvLAYpert,ppmvAVGpert,ppmvMAXpert,pavgLAYpert,tavgLAYpert,ppmv500pert,ppmv75pert,ppmvSURFpert] = layers2ppmv(h,pert,1:length(p.stemp),3);
scatter_coast(pert.rlon,pert.rlat,50,colwaterpert-colwater0); title('\delta colwater mm/yr');

i050 = find(p.plevs(1:90,2300) >= 050,1);
i100 = find(p.plevs(1:90,2300) >= 100,1);
i200 = find(p.plevs(1:90,2300) >= 200,1);
i500 = find(p.plevs(1:90,2300) >= 500,1);
i800 = find(p.plevs(1:90,2300) >= 800,1);

deltaRH = RHpert-RH0;
figure(25); simplemap(Y(:),X(:),deltaRH(i200,:)',5); colorbar; title(['d(RH)/dt yr-1 200 mb'])
figure(26); simplemap(Y(:),X(:),deltaRH(i500,:)',5); colorbar; title(['d(RH)/dt yr-1 500 mb'])
figure(27); simplemap(Y(:),X(:),deltaRH(i800,:)',5); colorbar; title(['d(RH)/dt yr-1 800 mb'])
for ii = 25 : 27; figure(ii); caxis([-1 +1]); colormap(usa2); plotaxis2; end

deltaO3 = ppmvLAYpert - ppmvLAY;
figure(25); simplemap(Y(:),X(:),deltaO3(i050,:)',5); colorbar; title(['d(O3)/dt ppm/yr 050 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(26); simplemap(Y(:),X(:),deltaO3(i100,:)',5); colorbar; title(['d(O3)/dt ppm/yr 100 mb']); caxis([-0.50 +0.50]/100); plotaxis2;
figure(27); simplemap(Y(:),X(:),deltaO3(i200,:)',5); colorbar; title(['d(O3)/dt ppm/yr 200 mb']); caxis([-0.50 +0.50]/100); plotaxis2;

deltaO3 = duO3_colpert - duO3_col0;
figure(25); simplemap(Y(:),X(:),deltaO3',5); colorbar; title(['d(O3)/dt column du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;
deltaO3 = duO3_col0;
figure(25); simplemap(Y(:),X(:),deltaO3',5); colorbar; title(['O3 column du']); caxis([250 350]); colormap jet; plotaxis2;
deltaO3 = duO3_300mbpert - duO3_300mb0;
figure(26); simplemap(Y(:),X(:),deltaO3',5); colorbar; title(['d(O3)/dt TOA->300mb  du/yr']); caxis([-0.05 +0.05]*10); plotaxis2;

%% https://acp.copernicus.org/articles/19/3257/2019/   tropospheric column O3

fluxX = load('olr_PERTv1.mat');
flux0 = load('olr_CLEAR_UNPERT.mat');
figure(27); simplemap(Y(:),X(:),fluxX.p2x.olrclr'-flux0.p2x.olrclr',5); colorbar; title(['d(clrOLR)/dt W/m2/yr']); caxis([-0.50 +0.50]); plotaxis2;
