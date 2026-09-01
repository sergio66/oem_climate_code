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
