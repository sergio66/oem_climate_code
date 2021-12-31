addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');
iAorOorL = input('Enter region : ');
if length(iAorOorL) == 0 
  iAorOorL = 0;
end

clear maskLF
maskLF = zeros(1,4608);
if iAorOorL == -7
  maskLF(abs(Ylat) > 60) = 1;
elseif iAorOorL == -6
  maskLF(abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -5
  maskLF(abs(Ylat) < 30) = 1;
elseif iAorOorL == 0
  maskLF = ones(1,4608);
elseif iAorOorL == -1
  maskLF(landfrac == 1) = 1;
elseif iAorOorL == +1
  maskLF(landfrac == 0) = 1;
elseif iAorOorL == -2
  maskLF(landfrac == 1 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == +2
  maskLF(landfrac == 0 & abs(Ylat) <= 30) = 1;
elseif iAorOorL == -3
  maskLF(landfrac == 1 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == +3
  maskLF(landfrac == 0 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
elseif iAorOorL == -4
  maskLF(landfrac == 1 & abs(Ylat) > 60) = 1;
elseif iAorOorL == +4
  maskLF(landfrac == 0 & abs(Ylat) > 60) = 1;
end
maskLFmatr = reshape(maskLF,72,64)';
mask = find(maskLF == 1);

figure(1); scatter_coast(X(:),Y(:),50,maskLF'.*landfrac); colorbar; title('landfrac');  caxis([-1 1]); colormap jet; hold on; plot(Xlon(maskLF == 1),Ylat(maskLF == 1),'.'); hold off

%{
figure(1); pcolor(X,Y,(reshape(results(:,1),72,64)).*maskLFmatr'); shading interp; colorbar; title('CO2');
figure(2); pcolor(X,Y,(reshape(results(:,2),72,64)).*maskLFmatr'); shading interp; colorbar; title('N2O');
figure(3); pcolor(X,Y,(reshape(results(:,3),72,64)).*maskLFmatr'); shading interp; colorbar; title('CH4');
figure(4); pcolor(X,Y,(reshape(results(:,4),72,64)).*maskLFmatr'); shading interp; colorbar; title('CFC11');
figure(5); pcolor(X,Y,(reshape(results(:,5),72,64)).*maskLFmatr'); shading interp; colorbar; title('CFC12');
figure(6); pcolor(X,Y,(reshape(results(:,6),72,64)).*maskLFmatr'); shading interp; colorbar; title('ST');
if iNorD > 0
  figure(7); pcolor(X,Y,data_trends.b_desc(:,:,1520).*maskLFmatr');  shading interp; colorbar; title('d/dt BT1231');
else
  figure(7); pcolor(X,Y,data_trends.b_asc(:,:,1520).*maskLFmatr');  shading interp; colorbar; title('d/dt BT1231');
end

figure(1); scatter_coast(X(:),Y(:),50,results(:,1).*maskLF'); shading interp; colorbar; title('CO2');  caxis([1.5 2.5])
figure(2); scatter_coast(X(:),Y(:),50,results(:,2).*maskLF'); shading interp; colorbar; title('N2O');  caxis([0 1])
figure(3); scatter_coast(X(:),Y(:),50,results(:,3).*maskLF'); shading interp; colorbar; title('CH4');  caxis([3 5])
figure(4); scatter_coast(X(:),Y(:),50,results(:,4).*maskLF'); shading interp; colorbar; title('CFC11'); caxis([-2 +2]*1e-10)
figure(5); scatter_coast(X(:),Y(:),50,results(:,5).*maskLF'); shading interp; colorbar; title('CFC12'); caxis([-5 +5]*1e-10)
figure(6); scatter_coast(X(:),Y(:),50,results(:,6).*maskLF'); shading interp; colorbar; title('ST');    caxis([-0.1 +0.1]/10)
if iNorD > 0
  figure(7); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,1520),72*64,1).*maskLF');  shading interp; colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])
else
  figure(7); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_asc(:,:,1520),72*64,1).*maskLF');  shading interp; colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])
end
figure(1); simplemap(Y(:),X(:),results(:,1).*maskLF',5); colorbar; title('CO2');  caxis([1.5 2.5])
figure(2); simplemap(Y(:),X(:),results(:,2).*maskLF',5); colorbar; title('N2O');  caxis([0 1])
figure(3); simplemap(Y(:),X(:),results(:,3).*maskLF',5); colorbar; title('CH4');  caxis([3 5])
figure(4); simplemap(Y(:),X(:),results(:,4).*maskLF',5); colorbar; title('CFC11'); caxis([-2 +2]*1e-3)
figure(5); simplemap(Y(:),X(:),results(:,5).*maskLF',5); colorbar; title('CFC12'); caxis([-5 +5]*1e-1)
figure(6); simplemap(Y(:),X(:),results(:,6).*maskLF',5); colorbar; title('ST');    caxis([-0.1 +0.1])
if iNorD > 0
  figure(7); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,1520),72*64,1).*maskLF',5);  colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])
else
  figure(7); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,1520),72*64,1).*maskLF',5);  colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])
end
%}

%opts = 'mercator';
figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,1),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt CO2');  caxis([1.5 2.5])
figure(2); aslmap(2,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,2),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt N2O');  caxis([0 1])
figure(3); aslmap(3,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,3),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt CH4');  caxis([3 5])
figure(4); aslmap(4,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,4),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt CFC11'); caxis([-2 +2]*1e-3)
figure(5); aslmap(5,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,5),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt CFC12'); caxis([-5 +5]*1e-1)
figure(6); aslmap(6,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,6),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt ST');    caxis([-0.15 +0.15])
if iNorD > 0
  figure(7); aslmap(7,rlat65,rlon73,maskLFmatr.*smoothn((data_trends.b_desc(:,:,1520)'),1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt BT1231'); caxis([-0.1 +0.1])
else
  figure(7); aslmap(7,rlat65,rlon73,maskLFmatr.*smoothn((data_trends.b_asc(:,:,1520)'),1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt BT1231'); caxis([-0.1 +0.1])
end
for ii = 1 : 7; figure(ii); plotaxis2; end

%{
umbc_st_trend4608 = reshape(results(:,6),72,64);
airs_quantile16_bt1231_trend4608 = data_trends.b_desc(:,:,1520);
comment = 'see driver_gather_gridded_retrieval_results.m';
save umbc_trends.mat airs_quantile16_bt1231_trend4608 umbc_st_trend4608 X Y comment rlat65 rlon73

figure(6); aslmap(6,rlat65,rlon73,maskLFmatr.*smoothn(umbc_st_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt UMBC'); caxis([-0.1 +0.1])
figure(7); aslmap(7,rlat65,rlon73,maskLFmatr.*smoothn(airs_quantile16_bt1231_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt BT1231'); caxis([-0.1 +0.1])

giss = load('giss_trends.mat');
figure(8); aslmap(8,rlat65,rlon73,maskLFmatr.*smoothn(giss.giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt GISS'); caxis([-0.1 +0.1])

airsL3 = load('airsL3_trends.mat');
figure(9); aslmap(9,rlat65,rlon73,maskLFmatr.*smoothn(airsL3.airsL3_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt AIRSL3 v7'); caxis([-0.1 +0.1])
%}

figure(1); plot(rlat,nanmean(reshape(results(:,1),72,64),1)); title('zonal CO2')
figure(2); plot(rlat,nanmean(reshape(results(:,2),72,64),1)); title('zonal N2O')
figure(3); plot(rlat,nanmean(reshape(results(:,3),72,64),1)); title('zonal CH4')
figure(4); plot(rlat,nanmean(reshape(results(:,4),72,64),1)); title('zonal CFC11')
figure(5); plot(rlat,nanmean(reshape(results(:,5),72,64),1)); title('zonal CFC12')
figure(6); plot(rlat,nanmean(reshape(results(:,6),72,64),1)); title('zonal Stemp')
if iNorD > 0
  figure(7); plot(rlat,nanmean(data_trends.b_desc(:,:,1520),1)); title('zonal BT1231')
else
  figure(7); plot(rlat,nanmean(data_trends.b_asc(:,:,1520),1)); title('zonal BT1231')
end
for ii = 1 : 7; figure(ii); colormap(llsmap5); plotaxis2; end

%{
figure(25); simplemap(Y(:),X(:),resultsWV(:,12),5); colorbar; title('WV 200 mb'); 
figure(26); simplemap(Y(:),X(:),resultsWV(:,16),5); colorbar; title('WV 500 mb');
figure(27); simplemap(Y(:),X(:),resultsWV(:,18),5); colorbar; title('WV 800 mb');
%}
figure(25); aslmap(25,rlat65,rlon73,maskLFmatr.*smoothn(reshape(resultsWV(:,12),72,64)',1),[-90 +90],[-180 +180]); colorbar; title('WV 200 mb'); 
figure(26); aslmap(26,rlat65,rlon73,maskLFmatr.*smoothn(reshape(resultsWV(:,16),72,64)',1),[-90 +90],[-180 +180]); colorbar; title('WV 500 mb');
figure(27); aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn(reshape(resultsWV(:,18),72,64)',1),[-90 +90],[-180 +180]); colorbar; title('WV 800 mb');
for ii = 25 : 27; figure(ii); caxis([-1e-2 +1e-2]); colormap(llsmap5); plotaxis2; end
%{
playsX = flipud(plays);
playsX = playsX(4:100);
wvjac = [];
junk = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');   wvjac = junk.jout;
junk = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g101_jac_new.mat'); wvjac = wvjac + junk.jout;
junk = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g102_jac_new.mat'); wvjac = wvjac + junk.jout;
junk = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g103_jac_new.mat'); wvjac = wvjac + junk.jout;
plot(junk.fout,sum(wvjac'))
for ii = 1 : 2378;
  boo = wvjac(ii,:);
  peak(ii) = playsX(find(abs(boo) == max(abs(boo)),1));
end
plot(junk.fout,peak,'x-'); set(gca,'ydir','reverse'); xlim([min(junk.fout) max(junk.fout)]); grid on
%}
ix200 = find(data_trends.h.vchan >= 1419,1);
ix500 = find(data_trends.h.vchan >= 1365,1);
ix800 = find(data_trends.h.vchan >= 0900,1);
ix200 = find(data_trends.h.vchan >= 1507,1); %% IASI
ix500 = find(data_trends.h.vchan >= 1441,1); %% IASI 
ix800 = find(data_trends.h.vchan >= 0900,1);
if iNorD > 0 
  figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix200),72*64,1),5); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix500),72*64,1),5); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix800),72*64,1),5); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(llsmap5); plotaxis2; end
  figure(25); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix200),72*64,1)); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix500),72*64,1)); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix800),72*64,1)); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(llsmap5); plotaxis2; end
else
  figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix200),72*64,1),5); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix500),72*64,1),5); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix800),72*64,1),5); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(llsmap5); plotaxis2; end
  figure(25); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_asc(:,:,ix200),72*64,1)); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_asc(:,:,ix500),72*64,1)); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_asc(:,:,ix800),72*64,1)); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(llsmap5); plotaxis2; end
end

figure(25); simplemap(Y(:),X(:),resultsT(:,12),5); colorbar; title('T 200 mb'); 
figure(26); simplemap(Y(:),X(:),resultsT(:,16),5); colorbar; title('T 500 mb');
figure(27); simplemap(Y(:),X(:),resultsT(:,18),5); colorbar; title('T 800 mb');
for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]); colormap(llsmap5); plotaxis2; end
%{
playsX = flipud(plays);
playsX = playsX(4:100);
tjac = [];
junk = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac_new.mat');   tjac = junk.jtemp;
plot(junk.fout,sum(tjac'))
for ii = 1 : 2378;
  boo = tjac(ii,:);
  peak(ii) = playsX(find(abs(boo) == max(abs(boo)),1));
end
plot(junk.fout,peak,'x-'); set(gca,'ydir','reverse'); xlim([min(junk.fout) max(junk.fout)]); grid on
%}
ix200 = find(data_trends.h.vchan >= 694,1);
ix500 = find(data_trends.h.vchan >= 730.7,1);
ix800 = find(data_trends.h.vchan >= 0814,1);
if iNorD > 0
  figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix200),72*64,1),5); colorbar; title(['T 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix500),72*64,1),5); colorbar; title(['T 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix800),72*64,1),5); colorbar; title(['T 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]); colormap(llsmap5); plotaxis2; end
else
  figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix200),72*64,1),5); colorbar; title(['T 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
  figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix500),72*64,1),5); colorbar; title(['T 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
  figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_asc(:,:,ix800),72*64,1),5); colorbar; title(['T 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
  for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]); colormap(llsmap5); plotaxis2; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');
[hMean17years,ha,pMean17years,pa]     = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
% see FIND_NWP_MODEL_TRENDS/driver_computeERA_16day_trends_desc_or_asc.m
iLoad = 1; 
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  end
  era_prof = load(fin);
  hTimeStep1 = era_prof.hnew_op;
  pTimeStep1 = era_prof.pnew_op;

if dataset == 3
  disp(' ')
  disp('Extremes, so augment stemp!!!')
  h = hTimeStep1;
  p = pTimeStep1;
  gev = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/GEVresults.mat');
  scatter_coast(reshape(p.rlon,72,64)',reshape(p.rlat,72,64)',30,squeeze(gev.paramE16(:,:,3)));
  scatter_coast(reshape(p.rlon,72,64)',reshape(p.rlat,72,64)',30,squeeze(gev.paramE16(:,:,3)-gev.paramQ16(:,:,3))); caxis([-5 +5])
  junk = squeeze(gev.paramE16(:,:,3)-gev.paramQ16(:,:,3)); junk = junk'; junk = junk(:); hist(junk)
    scatter_coast(p.rlon,p.rlat,30,junk); caxis([-5 +5]); colormap(jet); caxis([0 5])
  fprintf(1,'on avg the GEV fit says Extreme-Q16 = %8.4f +/- %8.4f K on land and %8.4f +/- %8.4f on ocean \n',[nanmean(junk(p.landfrac > 0)) nanstd(junk(p.landfrac > 0)) nanmean(junk(p.landfrac == 0)) nanstd(junk(p.landfrac == 0))])
  booL = find(isnan(junk) & p.landfrac' > 0);  junk(booL) = nanmean(junk(p.landfrac > 0)) + randn(size(booL))*nanstd(junk(p.landfrac > 0));
  booO = find(isnan(junk) & p.landfrac' == 0); junk(booO) = nanmean(junk(p.landfrac == 0)) + randn(size(booO))*nanstd(junk(p.landfrac == 0));
    scatter_coast(p.rlon,p.rlat,30,junk); caxis([-5 +5]); colormap(jet); caxis([0 5])
  pTimeStep1.stemp = pTimeStep1.stemp + junk';
  p.stemp = p.stemp + junk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gather_gridded_retrieval_results_plots

find_T_RH_trends
find_O3_trends

%% OLR trends
fluxX = load('../FIND_NWP_MODEL_TRENDS/olr_PERTv1.mat');
flux0 = load('../FIND_NWP_MODEL_TRENDS/olr_CLEAR_UNPERT.mat');
figure(27); simplemap(Y(:),X(:),(fluxX.p2x.olrclr'-flux0.p2x.olrclr').*maskLF',5); colorbar; title(['d(clrOLR)/dt W/m2/yr']); caxis([-0.50 +0.50]); plotaxis2;
aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn((reshape(fluxX.p2x.olrclr,72,64)'),1), [-90 +90],[-180 +180]);
aslmap(27,rlat65,rlon73,maskLFmatr.*smoothn((reshape(fluxX.p2x.olrclr-flux0.p2x.olrclr,72,64)'),1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt clrOLR W/m2/yr'); 
caxis([-0.5 +0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now look at other model data
look_at_other_model_data
