maskLFchan = ones(2645,1) * maskLF;
mask = find(maskLF == 1);

redo_fig8_spectralrates_with_mask

i730 = find(f >= 730,1);   figure(9);  clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i730,:)-fits(i730,:))); title('BT730 bias')
i900 = find(f >= 900,1);   figure(10); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i900,:)-fits(i900,:))); title('BT900 bias')
i1419 = find(f >= 1419,1); figure(11); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i1419,:)-fits(i1419,:))); title('BT1419 bias')
i667 = find(f >= 667,1);   figure(12); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i667,:)-fits(i667,:))); title('BT667 bias')
for ii = 9:12; figure(ii); caxis([-1 +1]*0.01); colormap(llsmap5); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% have added 4 new figures, so these subsequent figures displaced by 4 eg fig 11 --> fig 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chisqr = rates(jacobian.chanset,:)-fits(jacobian.chanset,:); chisqr = sum(chisqr.*chisqr,1)/length(jacobian.chanset); figure(13); clf; scatter_coast(X(:),Y(:),50,log10(chisqr.*maskLF)); title('log10(\chi^2) oem chans')

junkBT  = rates(i900,:).*maskLF;
junkT20 = resultsT(:,iNumLay)'.*maskLF;
junkW20 = resultsWV(:,iNumLay)'.*maskLF;
figure(14); plot(rlat,nanmean(reshape(junkBT,72,64),1),'g',rlat,10*nanmean(reshape(junkW20,72,64),1),'bx-',...
                rlat,nanmean(reshape(junkT20,72,64),1),'k',rlat,nanmean(reshape(results(:,6).*maskLF',72,64),1),'r','linewidth',2); title('zonal rates LAND+OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

land = (landfrac > 0); 
fprintf(1,'Out of %5i grid points, the fraction over land is %8.6f \n',length(landfrac),sum(land)/length(landfrac))

resultsNAN = results;               resultsNAN(land,:) = NaN;
junkBTNAN  = rates(i900,:);         junkBTNAN(land) = NaN;
junkT20NAN = resultsT(:,iNumLay)';  junkT20NAN(land) = NaN;
junkW20NAN = resultsWV(:,iNumLay)'; junkT20NAN(land) = NaN;
figure(15); plot(rlat,nanmean(reshape(junkBTNAN,72,64),1),'g',rlat,10*nanmean(reshape(junkW20NAN,72,64),1),'bx-',...
                rlat,nanmean(reshape(junkT20NAN,72,64),1),'k',rlat,nanmean(reshape(resultsNAN(:,6),72,64),1),'r','linewidth',2); title('zonal rates OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% these are NOT masked for L/O
figure(16); pcolor(rlon,pavg,resultsWV((1:72)+30*72,:)'); xlabel('longitude'); title('Equator WV'); caxis([-1e-2 +1e-2]);      
figure(17); pcolor(rlon,pavg,resultsT((1:72)+30*72,:)');  xlabel('longitude'); title('Equator T');  caxis([-0.1 +0.1]);        
figure(18); pcolor(rlon,pavg,resultsO3((1:72)+30*72,:)'); xlabel('longitude'); title('Equator O3'); caxis([-0.01 +0.01]*0.25); 
for ii = 16:18; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

%% these are NOT masked for L/O
ind = 1:72:4608;
ind = ind + 31;
figure(19); pcolor(rlat,pavg,resultsWV(ind,:)'); xlabel('latitude'); title('GMT line WV'); caxis([-1e-2 +1e-2]);      
figure(20); pcolor(rlat,pavg,resultsT(ind,:)');  xlabel('latitude'); title('GMT line T');  caxis([-0.1 +0.1]);        
figure(21); pcolor(rlat,pavg,resultsO3(ind,:)'); xlabel('latitude'); title('GMT line O3'); caxis([-0.01 +0.01]*0.25); 
for ii = 19:21; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

boo = zeros(iNumLay,72,64); for ijunk = 1 : iNumLay; boo(ijunk,:,:) = maskLFmatr'; end
xresultsWV = resultsWV'; xresultsWV = boo .* reshape(xresultsWV,iNumLay,72,64);
xresultsT  = resultsT';   xresultsT = boo .* reshape(xresultsT,iNumLay,72,64);
xresultsO3 = resultsO3'; xresultsO3 = boo .* reshape(xresultsO3,iNumLay,72,64);
figure(22); pcolor(rlat,pavg,squeeze(nanmean(xresultsWV,2))); xlabel('latitude'); title('zonal mean WV'); caxis([-1e-2 +1e-2]);      
figure(23); pcolor(rlat,pavg,squeeze(nanmean(xresultsT,2)));  xlabel('latitude'); title('zonal mean T');  caxis([-0.1 +0.1]);        
figure(24); pcolor(rlat,pavg,squeeze(nanmean(xresultsO3,2))); xlabel('latitude'); title('zonal mean O3'); caxis([-0.01 +0.01]*0.25);
for ii = 22:24; figure(ii); ylim([10 1000]); xlim([-90 +90]); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  figure(6); clf; scatter_coast(X(:),Y(:),50,results(:,6)); shading interp; colorbar; caxis([-0.2 +0.2])
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
