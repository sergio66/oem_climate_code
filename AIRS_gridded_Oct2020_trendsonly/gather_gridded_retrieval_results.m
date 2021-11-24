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
for ii = 1 : 64*72
  fname = ['Output/test' num2str(ii) '.mat'];
  if exist(fname)
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
figure(4); pcolor(X,Y,(reshape(results(:,4),72,64))); shading interp; colorbar; title('Cld1');
figure(5); pcolor(X,Y,(reshape(results(:,5),72,64))); shading interp; colorbar; title('Cld2');
figure(6); pcolor(X,Y,(reshape(results(:,6),72,64))); shading interp; colorbar; title('ST');

figure(1); scatter_coast(X(:),Y(:),50,results(:,1)); shading interp; colorbar; title('CO2');  caxis([0 2])
figure(2); scatter_coast(X(:),Y(:),50,results(:,2)); shading interp; colorbar; title('N2O');  caxis([0 1])
figure(3); scatter_coast(X(:),Y(:),50,results(:,3)); shading interp; colorbar; title('CH4');  caxis([0 4])
figure(4); scatter_coast(X(:),Y(:),50,results(:,4)); shading interp; colorbar; title('Cld1'); caxis([-2 +2]*1e-2)
figure(5); scatter_coast(X(:),Y(:),50,results(:,5)); shading interp; colorbar; title('Cld2'); caxis([-5 +5]*1e-3)
figure(6); scatter_coast(X(:),Y(:),50,results(:,6)); shading interp; colorbar; title('ST');    caxis([-0.2 +0.2])

figure(1); plot(rlat,nanmean(reshape(results(:,1),72,64),1)); title('zonal CO2')
figure(2); plot(rlat,nanmean(reshape(results(:,2),72,64),1)); title('zonal N2O')
figure(3); plot(rlat,nanmean(reshape(results(:,3),72,64),1)); title('zonal CH4')
figure(4); plot(rlat,nanmean(reshape(results(:,4),72,64),1)); title('zonal Cld1')
figure(5); plot(rlat,nanmean(reshape(results(:,5),72,64),1)); title('zonal Cld2')
figure(6); plot(rlat,nanmean(reshape(results(:,6),72,64),1)); title('zonal Stemp')

for ii = 1 : 6; figure(ii); colormap(usa2); plotaxis2; end

junk = load('h2645structure.mat');
f           = junk.h.vchan;
figure(7);  plot(f,nanmean(rates'),'b',f,nanmean(rates'-fits'),'r'); grid
i730 = find(f >= 730,1);   figure(8);  scatter_coast(X(:),Y(:),50,rates(i730,:)-fits(i730,:)); title('BT730 bias')
i900 = find(f >= 900,1);   figure(9);  scatter_coast(X(:),Y(:),50,rates(i900,:)-fits(i900,:)); title('BT900 bias')
i1419 = find(f >= 1419,1); figure(10); scatter_coast(X(:),Y(:),50,rates(i1419,:)-fits(i1419,:)); title('BT1419 bias')
i667 = find(f >= 667,1);   figure(11); scatter_coast(X(:),Y(:),50,rates(i667,:)-fits(i667,:)); title('BT667 bias')
for ii = 8:11; figure(ii); caxis([-0.02 +0.02]); colormap(usa2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% have added 4 new figures, so these subsequent figures displaced by 4 eg fig 11 --> fig 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chisqr = rates(jacobian.chanset,:)-fits(jacobian.chanset,:); chisqr = sum(chisqr.*chisqr,1)/length(jacobian.chanset); figure(20); scatter_coast(X(:),Y(:),50,log10(chisqr)); title('log10(\chi^2) oem chans')

junkBT  = rates(i900,:);
junkT20 = resultsT(:,iNumLay)';
junkW20 = resultsWV(:,iNumLay)';
figure(15); plot(rlat,nanmean(reshape(junkBT,72,64),1),'g',rlat,10*nanmean(reshape(junkW20,72,64),1),'bx-',...
                rlat,nanmean(reshape(junkT20,72,64),1),'k',rlat,nanmean(reshape(results(:,6),72,64),1),'r','linewidth',2); title('zonal rates LAND+OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

land = (landfrac > 0); 
fprintf(1,'Out of %5i grid points, the fraction over land is %8.6f \n',length(landfrac),sum(land)/length(landfrac))

resultsNAN = results;          resultsNAN(land,:) = NaN;
junkBTNAN  = rates(i900,:);    junkBTNAN(land) = NaN;
junkT20NAN = resultsT(:,iNumLay)';  junkT20NAN(land) = NaN;
junkW20NAN = resultsWV(:,iNumLay)'; junkT20NAN(land) = NaN;
figure(16); plot(rlat,nanmean(reshape(junkBTNAN,72,64),1),'g',rlat,10*nanmean(reshape(junkW20NAN,72,64),1),'bx-',...
                rlat,nanmean(reshape(junkT20NAN,72,64),1),'k',rlat,nanmean(reshape(resultsNAN(:,6),72,64),1),'r','linewidth',2); title('zonal rates OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

figure(17); pcolor(rlon,1:iNumLay,resultsWV((1:72)+30*72,:)'); xlabel('longitude'); title('Equator WV'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(18); pcolor(rlon,1:iNumLay,resultsT((1:72)+30*72,:)');  xlabel('longitude'); title('Equator T');  colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(19); pcolor(rlon,1:iNumLay,resultsO3((1:72)+30*72,:)'); xlabel('longitude'); title('Equator O3'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');

ind = 1:72:4608;
ind = ind + 31;
figure(17); pcolor(rlat,1:iNumLay,resultsWV(ind,:)'); xlabel('latitude'); title('GMT line WV'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(18); pcolor(rlat,1:iNumLay,resultsT(ind,:)');  xlabel('latitude'); title('GMT line T');  colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(19); pcolor(rlat,1:iNumLay,resultsO3(ind,:)'); xlabel('latitude'); title('GMT line O3'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');

xresultsWV = resultsWV'; xresultsWV = reshape(xresultsWV,iNumLay,72,64);
xresultsT  = resultsT';   xresultsT = reshape(xresultsT,iNumLay,72,64);
xresultsO3 = resultsO3'; xresultsO3 = reshape(xresultsO3,iNumLay,72,64);
figure(17); pcolor(rlat,1:iNumLay,squeeze(nanmean(xresultsWV,2))); xlabel('latitude'); title('zonal mean WV'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse'); 
figure(18); pcolor(rlat,1:iNumLay,squeeze(nanmean(xresultsT,2)));  xlabel('latitude'); title('zonal mean T');  colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(19); pcolor(rlat,1:iNumLay,squeeze(nanmean(xresultsO3,2))); xlabel('latitude'); title('zonal mean O3'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse');
figure(20); pcolor(rlat,pavg,squeeze(nanmean(xresultsWV,2))); xlabel('latitude'); title('zonal mean WV'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log')
figure(21); pcolor(rlat,pavg,squeeze(nanmean(xresultsT,2)));  xlabel('latitude'); title('zonal mean T');  colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log')
figure(22); pcolor(rlat,pavg,squeeze(nanmean(xresultsO3,2))); xlabel('latitude'); title('zonal mean O3'); colorbar; shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log')
for ii = 17:22; figure(ii); ylim([10 1000]); xlim([-60 +60]); caxis([-0.01 +0.01]); end; figure(18); caxis([-0.05 +0.05]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSave = -1;
if iSave > 0
  %% figs 17,18,5
  stm.WVrates = squeeze(nanmean(xresultsWV,2));
  stm.Trates  = squeeze(nanmean(xresultsT,2));
  stm.O3rates = squeeze(nanmean(xresultsO3,2));
  stm.rlat    = rlat;
  stm.pavg    = pavg;

  %% fig 7
  stm.f = f;
  stm.meanrates = nanmean(rates');
  stm.meanfits  = nanmean(fits');

  %% fig 15
  stm.landNocean_BT900 = nanmean(reshape(junkBT,72,64),1);
  stm.landNocean_sst   = nanmean(reshape(results(:,6),72,64),1);
  stm.landNocean_co2   = nanmean(reshape(results(:,1),72,64),1);
  stm.landNocean_n2o   = nanmean(reshape(results(:,2),72,64),1);
  stm.landNocean_ch4   = nanmean(reshape(results(:,3),72,64),1);

  stm.ocean_BT900 = nanmean(reshape(junkBTNAN,72,64),1);
  stm.ocean_sst   = nanmean(reshape(resultsNAN(:,6),72,64),1);
  stm.ocean_co2   = nanmean(reshape(resultsNAN(:,1),72,64),1);
  stm.ocean_n2o   = nanmean(reshape(resultsNAN(:,2),72,64),1);
  stm.ocean_ch4   = nanmean(reshape(resultsNAN(:,3),72,64),1);

  %% fig 17
  stm.chisqr      = nanmean(reshape(chisqr,72,64),1);

  fnameout = 'gather_results_oct2_2015.mat';
  save(fnameout,'-struct','stm');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPrint = -1;
if iPrint > 0

  lls4 = load('llsmap4');
  lls4.llsmap4 = lls4.llsmap4(2:end,:);

  figure(17); pcolor(rlat,pavg,squeeze(nanmean(xresultsWV,2))); xlabel('latitude'); ylabel('Pressure (mb)'); shading interp; colormap(lls4.llsmap4); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  figure(18); pcolor(rlat,pavg,squeeze(nanmean(xresultsT,2)));  xlabel('latitude'); ylabel('Pressure (mb)'); shading interp; colormap(lls4.llsmap4); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  figure(19); pcolor(rlat,pavg,squeeze(nanmean(xresultsO3,2))); xlabel('latitude'); ylabel('Pressure (mb)'); shading interp; colormap(lls4.llsmap4); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  for ii = 17:19; figure(ii); ylim([10 1000]); xlim([-90 +90]); caxis([-0.01 +0.01]); colorbar; end; figure(18); caxis([-0.05 +0.05]); 

  rrlon = X(:); rrlat = Y(:); airs_stemprate = results(:,6);
  figure(6); scatter_coast(X(:),Y(:),50,results(:,6)); shading interp; colorbar; caxis([-0.2 +0.2])
  colormap(lls4.llsmap4); ylabel('Longitude'); xlabel('Latitude'); caxis([-0.2 +0.2]);
  % save strow_retrievals_quicklook_oct4 X Y results

  addpath /asl/matlib/plotutils
  figure(18); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Oct2020/GriddedRetrievals/Figs/airsL1c_trends_Tz_vs_lat.pdf');
  figure(17); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Oct2020/GriddedRetrievals/Figs/airsL1c_trends_WVz_vs_lat.pdf');
  figure(19); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Oct2020/GriddedRetrievals/Figs/airsL1c_trends_O3z_vs_lat.pdf');
  figure(06); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Oct2020/GriddedRetrievals/Figs/airsL1c_trends_stemp_vs_lat_vs_lon.pdf');

  figure(15); plot(rlat,nanmean(reshape(results(:,6),72,64),1),'r','linewidth',2);    plotaxis2; ylim([-0.04 +0.04]); xlim([-90 +90]); %title('zonal rates LAND+OCEAN'); 
  figure(16); plot(rlat,nanmean(reshape(resultsNAN(:,6),72,64),1),'r','linewidth',2); plotaxis2; ylim([-0.04 +0.04]); xlim([-90 +90]); %title('zonal rates OCEAN'); 
  for ii = 15 : 16; xlabel('Latitude'); ylabel('dST/dt (K/yr)'); end
  figure(12); aslprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Oct2020/GriddedRetrievals/Figs/airsL1c_stemp_ocean_vs_lat.pdf');  

  addpath /home/sergio/MATLABCODE/NANROUTINES
  r_1_x(1) = nanlinearcorrelation(results(:,1),results(:,1));
  r_1_x(2) = nanlinearcorrelation(results(:,1),results(:,2));
  r_1_x(3) = nanlinearcorrelation(results(:,1),results(:,3));
  r_1_x(4) = nanlinearcorrelation(results(:,1),results(:,4));
  r_1_x(5) = nanlinearcorrelation(results(:,1),results(:,5));
  r_1_x(6) = nanlinearcorrelation(results(:,1),results(:,6));

  addpath /asl/matlib/h4tools
  [h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');
  r_prof(1) = nanlinearcorrelation(results(:,1),p.cngwat');
  r_prof(2) = nanlinearcorrelation(results(:,1),p.cprtop');
  r_prof(3) = nanlinearcorrelation(results(:,1),p.cfrac');
  r_prof(4) = nanlinearcorrelation(results(:,1),p.cpsize');
  r_prof(5) = nanlinearcorrelation(results(:,1),p.cngwat2');
  r_prof(6) = nanlinearcorrelation(results(:,1),p.cprtop2');
  r_prof(7) = nanlinearcorrelation(results(:,1),p.cfrac2');
  r_prof(8) = nanlinearcorrelation(results(:,1),p.cpsize2');

  plot(r_prof,'o','Markersize',5); ylabel('Correlation woth CO2 rate'); grid
  names = {'Ice Amt'; 'Ice Top'; 'Ice Frac'; 'Ice Size'; 'Water Amt'; 'Water Top'; 'Water Frac'; 'Water Size';};
  set(gca,'xtick',[1:8],'xticklabel',names,'fontsize',10); xtickangle(45); plotaxis2;
end
