%% see plot_driver_gather_gridded_retrieval_results.m

WV_Ylim = 100;
T_Ylim  = 10;

if iOCBset == 0
  ocbstr = 'UMBC';
elseif iOCBset == 1
  ocbstr = 'ERA5 synthetic';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aslmap(4,rlat65,rlon73,smoothn((reshape(results(:,1),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt CO2');  caxis([1.5 2.5])
clf;; scatter_coast(Xlon,Ylat,50,results(:,1)); title('d/dt CO2');  caxis([1.5 2.5]); caxis([2.0 2.5])

aslmap(6,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,6),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt ST');    caxis([-0.15 +0.15])

figure(9); scatter_coast(Xlon,Ylat,50,thedofs); jett = jet(64); jett(1,:) = 1; colormap(jett); title('ALL DOFS'); caxis([0 max(thedofs)]); caxis([0 30])
figure(10); boo = find(lencdofs == 66); sumc = mean(cdofs(boo,:),1); 
  plot(cumsum(sumc)); fl = ceil(sum(sumc)); line([6 6],[0 fl],'color','k'); line([26 26],[0 fl],'color','k'); line([46 46],[0 fl],'color','k');
  text(2,10,'TG','fontsize',10);   text(16,10,'WV','fontsize',10);   text(36,10,'T','fontsize',10);   text(56,10,'O3','fontsize',10); xlim([0 66])

figure(11); semilogy(wvsumcflip,pflip20,tsumcflip,pflip20,o3sumcflip,pflip20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([0.050 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Enter (-7) polar    L/O')
disp('      (-6) midlat   L/O')
disp('      (-5) tropical L/O')
disp('      (-4) polar land    (+4) polar ocean')
disp('      (-3) midlat land   (+3) midlat ocean')
disp('      (-2) tropical land (+2) tropical ocean')
disp('      (-1) land          (+1) ocean');
disp('      [0,default] ALL trends : ');
ixAorOorL = input('Enter region : ');
if length(ixAorOorL) == 0 
  ixAorOorL = 0;
end

clear xmaskLF ymaskLF zmaskLF
xmaskLF = zeros(1,4608);
ymaskLF = zeros(1,4608);
zmaskLF = zeros(1,4608);

co2minrate = 2.15;

if ixAorOorL == -7
  ymaskLF(abs(Ylat) > 60) = 1;
  zmaskLF(results(:,1) > co2minrate & abs(Ylat) > 60) = 1;
  xstrmaskLF = '_polar_';
elseif ixAorOorL == -6
  ymaskLF(abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  zmaskLF(results(:,1) > co2minrate & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  xstrmaskLF = '_midlat_';
elseif ixAorOorL == -5
  ymaskLF(abs(Ylat) < 30) = 1;
  zmaskLF(results(:,1) > co2minrate & abs(Ylat) < 30) = 1;
  xstrmaskLF = '_tropical_';
elseif ixAorOorL == 0
  ymaskLF = ones(1,4608);
  zmaskLF(results(:,1) > co2minrate) = 1;
  xstrmaskLF = '_all_';
elseif ixAorOorL == -1
  ymaskLF(landfrac == 1) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 1) = 1;
  xstrmaskLF = '_land_';
elseif ixAorOorL == +1
  ymaskLF(landfrac == 0) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 0) = 1;
  xstrmaskLF = '_ocean_';
elseif ixAorOorL == -2
  ymaskLF(landfrac == 1 & abs(Ylat) <= 30) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 1 & abs(Ylat) <= 30) = 1;
  xstrmaskLF = '_tropical_land_';  
elseif ixAorOorL == +2
  ymaskLF(landfrac == 0 & abs(Ylat) <= 30) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 0 & abs(Ylat) <= 30) = 1;
  xstrmaskLF = '_tropical_ocean_';  
elseif ixAorOorL == -3
  ymaskLF(landfrac == 1 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 1 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  xstrmaskLF = '_midlat_land_';  
elseif ixAorOorL == +3
  ymaskLF(landfrac == 0 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 0 & abs(Ylat) > 30 & abs(Ylat) <= 60) = 1;
  xstrmaskLF = '_midlat_ocean_';  
elseif ixAorOorL == -4
  ymaskLF(landfrac == 1 & abs(Ylat) > 60) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 1 & abs(Ylat) > 60) = 1;
  xstrmaskLF = '_polar_land_';  
elseif ixAorOorL == +4
  ymaskLF(landfrac == 0 & abs(Ylat) > 60) = 1;
  zmaskLF(results(:,1) > co2minrate & landfrac == 0 & abs(Ylat) > 60) = 1;
  xstrmaskLF = '_polar_ocean_';  
end

if iNorD == 1
  xstrmaskLF = [xstrmaskLF 'N_'];
else
  xstrmaskLF = [xstrmaskLF 'D_'];
end

fprintf(1,'using only landfrac/rlat               finds %5i tiles \n',length(find(ymaskLF == 1)))
fprintf(1,'using both landfrac/rlat and CO2 > 2.1 finds %5i tiles \n',length(find(zmaskLF == 1)))
junk = input('Enter (-1) landfrac/rlat or (+1) landfrac/rlat and CO2 > 2.1 (default(+1)) : ');
if length(junk) == 0
  junk = 1;
end
if junk == -1
  xmaskLF = ymaskLF;
  xmaskLFmatr = reshape(ymaskLF,72,64)';  %% simple behavior, just based on LF and rlat
  xmask = find(ymaskLF == 1);
else
  xmaskLF = zmaskLF;
  xmaskLFmatr = reshape(zmaskLF,72,64)';  %% simple behavior, just based on LF and rlat and CO2
  xmask = find(zmaskLF == 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see find_T_RH_trends.m
clear ta taHandle tafov; zonal_trends_umbc_airsL3_models_ind_figures
clear ta taHandle tafov; zonal_trends_umbc_airsL3_models_tiles

fprintf(1,'<cos(p.rlat(xmask))> = %8.6f \n',mean(cos(p.rlat(xmask)*pi/180)))
coswgt101 = ones(101,1) * cos(p.rlat(xmask)*pi/180);
coswgt100 = ones(100,1) * cos(p.rlat(xmask)*pi/180);
coswgt097 = ones(097,1) * cos(p.rlat(xmask)*pi/180);
coswgt012 = ones(012,1) * cos(p.rlat(xmask)*pi/180);
coswgt024 = ones(024,1) * cos(p.rlat(xmask)*pi/180);

mncos101 = 1./nanmean(coswgt101,2);
mncos100 = 1./nanmean(coswgt100,2);
mncos097 = 1./nanmean(coswgt097,2);
mncos012 = 1./nanmean(coswgt012,2);
mncos024 = 1./nanmean(coswgt024,2);

clear ta taHandle tafov; mean_wgt_trends_umbc_airsL3_models_figures
clear ta taHandle tafov; mean_wgt_trends_umbc_airsL3_models_tiles
clear ta taHandle tafov; mean_wgt_trends_umbc_airsL3_models_tiles_unc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFig = 39;
figure(iFig); clf
  plot(f,nanmean(rates(:,xmask),2),'k',f,nanmean(fits(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,xmask),2),'linewidth',2)
  plot(f,nanmean(rates(:,xmask),2),'k',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,xmask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,xmask),2),'linewidth',2)
  plotaxis2; hl = legend('AIRS Obs',ocbstr,'AIRS L3','CMIP6','ERA5','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

  plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,xmask),2),'x-',...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,xmask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,xmask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,xmask),2),f,nanmean(rates(:,xmask),2),'linewidth',1)
  plotaxis2; hl = legend(ocbstr,'AIRS L3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',8); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

cos2645 = ones(2645,1)*cos(p.rlat(xmask)*pi/180);
cos2645 = ones(size(cos2645));
mncos2645 = 1./nanmean(cos2645,2);
  plot(f,mncos2645.*nanmean(cos2645.*nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,xmask),2),'x-',...
       f,mncos2645.*nanmean(cos2645.*nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,xmask),2),...
       f,mncos2645.*nanmean(cos2645.*nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,xmask),2),...
       f,mncos2645.*nanmean(cos2645.*nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,xmask),2),...
       f,mncos2645.*nanmean(cos2645.*rates(:,xmask),2),'linewidth',1)
  plotaxis2; hl = legend(ocbstr,'AIRS L3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',8); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = input('print the figs +1/-1   (default = -1) : ');
if length(iPrint) == 0
  iPrint = -1;
end
addpath /asl/matlib/plotutils
if iPrint > 0
  printdir = 'Figs_STMOct2021/';

%{

  figure(28); figname = [printdir '/umbc' xstrmaskLF 'T_trend.pdf'];   aslprint(figname)
  figure(46); figname = [printdir '/era' xstrmaskLF 'T_trend.pdf'];    aslprint(figname)
  figure(47); figname = [printdir '/era5' xstrmaskLF 'T_trend.pdf'];   aslprint(figname)
  figure(48); figname = [printdir '/airsL3' xstrmaskLF 'T_trend.pdf']; aslprint(figname)
  figure(49); figname = [printdir '/cmip6' xstrmaskLF 'T_trend.pdf'];  aslprint(figname)

  figure(29); figname = [printdir '/umbc' xstrmaskLF 'RH_trend.pdf'];   aslprint(figname)
  figure(50); figname = [printdir '/era' xstrmaskLF 'RH_trend.pdf'];    aslprint(figname)
  figure(51); figname = [printdir '/era5' xstrmaskLF 'RH_trend.pdf'];   aslprint(figname)
  figure(52); figname = [printdir '/airsL3' xstrmaskLF 'RH_trend.pdf']; aslprint(figname)
  figure(53); figname = [printdir '/cmip6' xstrmaskLF 'RH_trend.pdf'];  aslprint(figname)

  figure(30); figname = [printdir '/umbc' xstrmaskLF 'WVfrac_trend.pdf'];   aslprint(figname)
  figure(54); figname = [printdir '/era' xstrmaskLF 'WVfrac_trend.pdf'];    aslprint(figname)
  figure(55); figname = [printdir '/era5' xstrmaskLF 'WVfrac_trend.pdf'];   aslprint(figname)
  figure(56); figname = [printdir '/airsL3' xstrmaskLF 'WVfrac_trend.pdf']; aslprint(figname)
  figure(57); figname = [printdir '/cmip6' xstrmaskLF 'WVfrac_trend.pdf'];  aslprint(figname)
%}


  figure(33); figname = [printdir '/tiled' xstrmaskLF 'WVfrac_trend.pdf'];     aslprint(figname)
  figure(34); figname = [printdir '/tiled' xstrmaskLF 'RH_trend.pdf'];         aslprint(figname)

%{
  figure(32); figname = [printdir '/tiled' xstrmaskLF 'T_trend.pdf'];          aslprint(figname)
  figure(33); figname = [printdir '/tiled' xstrmaskLF 'WVfrac_trend.pdf'];     aslprint(figname)
  figure(34); figname = [printdir '/tiled' xstrmaskLF 'RH_trend.pdf'];         aslprint(figname)
  figure(39); figname = [printdir '/spectral' xstrmaskLF 'rates.pdf'];         aslprint(figname)
  figure(40); figname = [printdir '/tiled' xstrmaskLF 'T_RH_WVfracrates.pdf']; aslprint(figname)
%}

end

