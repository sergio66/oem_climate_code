addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/science/

for ii = 1 : 15
  figure(ii); colormap jet
end

iNumLay = 33;
iNumLay = 20;

clear results
iaFound = [];
resultsWV = ones(1,iNumLay);
resultsT  = ones(1,iNumLay);
resultsO3 = ones(1,iNumLay);
iWarning = 0;
iLatbin = input('Enter latbin (1:64) : ');
if length(iLatbin) == 0
  iLatbin = 32;
end

for iii = 1 : 72
  ii = (iLatbin-1)*72 + iii;
  indexsave(iii) = ii;
  fname = ['Output/test' num2str(ii) '.mat'];
  if exist(fname)
    iExist = +1;
    loader = ['load ' fname];
    eval(loader);
    iaFound(iii) = +1;
    results(iii,1:6) = oem.finalrates(1:6);
    [mmn,nn] = size(oem.ak_water);

    nlays(ii) = nn;
    nn0 = min(nn,iNumLay);
    if nn0 == nn
      resultsWV(iii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
      resultsT(iii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
      resultsO3(iii,1:nn) = oem.finalrates((1:nn)+6+nn*2);
    else
      iWarning = iWarning + 1
      iaWarning(iWarning) = iii;
      wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(iii,1:nn0) = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(iii,1:nn0) = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(iii,1:nn0) = wah(1:nn0);
    end

    rates(:,iii) = rateset.rates;
    fits(:,iii)  = oem.fit';
  else
    iExist = -1;
    iaFound(iii) = 0;
    results(iii,1:6) = NaN;
    resultsWV(iii,:) = NaN;
    resultsT(iii,:)  = NaN;
    resultsO3(iii,:) = NaN;
    rates(:,iii) = NaN;
    fits(:,iii) = NaN;
  end
end
fprintf(1,'\n');

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
do_XX_YY_from_X_Y
[mmjunk,nnjunk] = size(XX);
if mmjunk == 1
  XX = XX';
  YY = YY';
end

addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
figure(1); scatter_coast(X(:),Y(:),50,landfrac); colorbar; title('landfrac');  caxis([0 1])

figure(1); plot(XX(indexsave),results(:,1)); title('zonal CO2')
figure(2); plot(XX(indexsave),results(:,2)); title('zonal N2O')
figure(3); plot(XX(indexsave),results(:,3)); title('zonal CH4')
figure(4); plot(XX(indexsave),results(:,4)); title('zonal Cld1')
figure(5); plot(XX(indexsave),results(:,5)); title('zonal Cld2')
figure(6); plot(XX(indexsave),results(:,6)); title('zonal ST')

junk = load('h2645structure.mat');
f           = junk.h.vchan;
figure(7);  plot(f,nanmean(rates'),'b',f,nanmean(rates'-fits'),'r'); grid
i730 = find(f >= 730,1);   figure(8);  plot(XX(indexsave),rates(i730,:)-fits(i730,:)); title('BT730 bias')
i900 = find(f >= 900,1);   figure(9);  plot(XX(indexsave),rates(i900,:)-fits(i900,:)); title('BT900 bias')
i1419 = find(f >= 1419,1); figure(10); plot(XX(indexsave),rates(i1419,:)-fits(i1419,:)); title('BT1419 bias')
i667 = find(f >= 667,1);   figure(11); plot(XX(indexsave),rates(i667,:)-fits(i667,:)); title('BT667 bias');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% have added 4 new figures, so these subsequent figures displaced by 4 eg fig 11 --> fig 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chisqr = rates(jacobian.chanset,:)-fits(jacobian.chanset,:); chisqr = sum(chisqr.*chisqr,1)/length(jacobian.chanset); 
figure(20); plot(XX(indexsave),log10(chisqr)); title('log10(\chi^2) oem chans')

figure(17); pcolor(XX(indexsave),1:iNumLay,resultsWV'); xlabel('longitude'); title('Equator WV'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(18); pcolor(XX(indexsave),1:iNumLay,resultsT');  xlabel('longitude'); title('Equator T');  colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(19); pcolor(XX(indexsave),1:iNumLay,resultsO3'); xlabel('longitude'); title('Equator O3'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');

