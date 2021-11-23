addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/matlib/science/
addpath /asl/matlib/aslutil

for ii = 1 : 22
  figure(ii); colormap jet
end

iNumLay = 33;
iNumLay = 20;

data_trends = load('convert_sergio_clearskygrid_obsonly.mat');

if ~exist('iaFound')
  clear results*
  iaFound = zeros(1,4608);
  resultsWV = ones(1,iNumLay);
  resultsT  = ones(1,iNumLay);
  resultsO3 = ones(1,iNumLay);
end

iWarning = 0;
for ii = 1 : 64*72
  fname = ['Output/test' num2str(ii) '.mat'];
  if exist(fname) & iaFound(ii) == 0
    iExist = +1;
    loader = ['load ' fname];
    eval(loader);
    iaFound(ii) = +1;
    results(ii,1:6) = oem.finalrates(1:6);
    [mmn,nn] = size(oem.ak_water);

    nlays(ii) = nn;
    nn0 = min(nn,iNumLay);
    if nn0 == nn
      resultsWV(ii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
      resultsT(ii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
      resultsO3(ii,1:nn) = oem.finalrates((1:nn)+6+nn*2);
    else
      iWarning = iWarning + 1
      iaWarning(iWarning) = ii;
      wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(ii,1:nn0) = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(ii,1:nn0) = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(ii,1:nn0) = wah(1:nn0);
    end

    rates(:,ii) = rateset.rates;
    fits(:,ii)  = oem.fit';
  else
    iExist = -1;
    iaFound(ii) = 0;
    results(ii,1:6) = NaN;
    resultsWV(ii,:) = NaN;
    resultsT(ii,:)  = NaN;
    resultsO3(ii,:) = NaN;
    rates(:,ii) = NaN;
    fits(:,ii) = NaN;
  end
  if mod(ii,72) == 0
    fprintf(1,'+ %2i\n',ii/72);
  elseif iExist == 1
    fprintf(1,'.');
  elseif iExist == -1
    fprintf(1,' ');
  end
end
fprintf(1,'\n');
fprintf(1,'found %4i of %4i \n',sum(iaFound),64*72)

%%%%%%%%%%%%%%%%%%%%%%%%%
plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);

clear pavg
for ii = 1 : iNumLay
  %iavg = jacobian.wvjaclays_used{iNumLay}-6;
  iavg = jacobian.wvjaclays_used{ii}-jacobian.wvjaclays_offset;
  pavg(ii) = mean(plays(iavg));
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
figure(1); scatter_coast(X(:),Y(:),50,landfrac); colorbar; title('landfrac');  caxis([0 1])

figure(1); pcolor(X,Y,(reshape(results(:,1),72,64))); shading interp; colorbar; title('CO2');
figure(2); pcolor(X,Y,(reshape(results(:,2),72,64))); shading interp; colorbar; title('N2O');
figure(3); pcolor(X,Y,(reshape(results(:,3),72,64))); shading interp; colorbar; title('CH4');
figure(4); pcolor(X,Y,(reshape(results(:,4),72,64))); shading interp; colorbar; title('CFC11');
figure(5); pcolor(X,Y,(reshape(results(:,5),72,64))); shading interp; colorbar; title('CFC12');
figure(6); pcolor(X,Y,(reshape(results(:,6),72,64))); shading interp; colorbar; title('ST');
figure(7); pcolor(X,Y,data_trends.b_desc(:,:,1520));  shading interp; colorbar; title('d/dt BT1231');

figure(1); scatter_coast(X(:),Y(:),50,results(:,1)); shading interp; colorbar; title('CO2');  caxis([1.5 2.5])
figure(2); scatter_coast(X(:),Y(:),50,results(:,2)); shading interp; colorbar; title('N2O');  caxis([0 1])
figure(3); scatter_coast(X(:),Y(:),50,results(:,3)); shading interp; colorbar; title('CH4');  caxis([3 5])
figure(4); scatter_coast(X(:),Y(:),50,results(:,4)); shading interp; colorbar; title('CFC11'); caxis([-2 +2]*1e-10)
figure(5); scatter_coast(X(:),Y(:),50,results(:,5)); shading interp; colorbar; title('CFC12'); caxis([-5 +5]*1e-10)
figure(6); scatter_coast(X(:),Y(:),50,results(:,6)); shading interp; colorbar; title('ST');    caxis([-0.1 +0.1]/10)
figure(7); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,1520),72*64,1));  shading interp; colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])

figure(1); simplemap(Y(:),X(:),results(:,1),5); colorbar; title('CO2');  caxis([1.5 2.5])
figure(2); simplemap(Y(:),X(:),results(:,2),5); colorbar; title('N2O');  caxis([0 1])
figure(3); simplemap(Y(:),X(:),results(:,3),5); colorbar; title('CH4');  caxis([3 5])
figure(4); simplemap(Y(:),X(:),results(:,4),5); colorbar; title('CFC11'); caxis([-2 +2]*1e-3)
figure(5); simplemap(Y(:),X(:),results(:,5),5); colorbar; title('CFC12'); caxis([-5 +5]*1e-1)
figure(6); simplemap(Y(:),X(:),results(:,6),5); colorbar; title('ST');    caxis([-0.1 +0.1])
figure(7); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,1520),72*64,1),5);  colorbar; title('d/dt BT1231'); caxis([-0.1 +0.1])

figure(1); plot(rlat,nanmean(reshape(results(:,1),72,64),1)); title('zonal CO2')
figure(2); plot(rlat,nanmean(reshape(results(:,2),72,64),1)); title('zonal N2O')
figure(3); plot(rlat,nanmean(reshape(results(:,3),72,64),1)); title('zonal CH4')
figure(4); plot(rlat,nanmean(reshape(results(:,4),72,64),1)); title('zonal CFC11')
figure(5); plot(rlat,nanmean(reshape(results(:,5),72,64),1)); title('zonal CFC12')
figure(6); plot(rlat,nanmean(reshape(results(:,6),72,64),1)); title('zonal Stemp')
figure(7); plot(rlat,nanmean(data_trends.b_desc(:,:,1520),1)); title('zonal BT1231')
for ii = 1 : 7; figure(ii); colormap(usa2); plotaxis2; end

figure(25); simplemap(Y(:),X(:),resultsWV(:,12),5); colorbar; title('WV 200 mb'); 
figure(26); simplemap(Y(:),X(:),resultsWV(:,16),5); colorbar; title('WV 500 mb');
figure(27); simplemap(Y(:),X(:),resultsWV(:,18),5); colorbar; title('WV 800 mb');
for ii = 25 : 27; figure(ii); caxis([-1e-2 +1e-2]); colormap(usa2); plotaxis2; end
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
figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix200),72*64,1),5); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix500),72*64,1),5); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix800),72*64,1),5); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(usa2); plotaxis2; end
figure(25); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix200),72*64,1)); colorbar; title(['WV 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
figure(26); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix500),72*64,1)); colorbar; title(['WV 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
figure(27); scatter_coast(X(:),Y(:),50,reshape(data_trends.b_desc(:,:,ix800),72*64,1)); colorbar; title(['WV 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]/2); colormap(usa2); plotaxis2; end

figure(25); simplemap(Y(:),X(:),resultsT(:,12),5); colorbar; title('T 200 mb'); 
figure(26); simplemap(Y(:),X(:),resultsT(:,16),5); colorbar; title('T 500 mb');
figure(27); simplemap(Y(:),X(:),resultsT(:,18),5); colorbar; title('T 800 mb');
for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]); colormap(usa2); plotaxis2; end
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
figure(25); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix200),72*64,1),5); colorbar; title(['T 200 mb = ' num2str(data_trends.h.vchan(ix200)) ' cm-1']); 
figure(26); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix500),72*64,1),5); colorbar; title(['T 500 mb = ' num2str(data_trends.h.vchan(ix500)) ' cm-1']); 
figure(27); simplemap(Y(:),X(:),reshape(data_trends.b_desc(:,:,ix800),72*64,1),5); colorbar; title(['T 800 mb = ' num2str(data_trends.h.vchan(ix800)) ' cm-1']); 
for ii = 25 : 27; figure(ii); caxis([-1e-1 +1e-1]); colormap(usa2); plotaxis2; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
gather_gridded_retrieval_results_plots

find_RH_trends
