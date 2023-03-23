figure(49);  %% 3 panel dT/dt       UMBC, ERA5xAK,ERA5 
figure(50);  %% 3 panel dWVfrac/dt  UMBC, ERA5xAK,ERA5 
figure(51);  %% 2 panel dRH/dt      UMBC, ERA5 

% plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
% plotoptions.yLimits = [100 1000];
% plotoptions.yLinearOrLog = +1;
% plotoptions.str2 = 'ERA5';   
% if isfield(plotoptions,'str3')
%   plotoptions = rmfield(plotoptions,'str3');
% end
% z1 = deltaRHlat'; 
% z2 = era5.trend_RH; z2 = reshape(z2,100,72,64); z2 = squeeze(nanmean(z2,2));
% iFig = 51; figure(iFig); clf; profile_plots_2tiledlayout(rlat,plays,z1,z2,iFig,plotoptions);

%% can estimate max deltaRH expected from Isaac Held blog  https://www.gfdl.noaa.gov/blog_held/47-relative-humidity-over-the-oceans/
%%% see my notes, BK 45
averageTlat = squeeze(nanmean(reshape(p.ptemp,[101 72 64]),2))';
      Lo = 2.5e6;  %%% J/kg
      Rv = 461.52; %%% J/kg/K
      moo = exp(Lo/Rv * deltaTlat ./ averageTlat ./ averageTlat) - 1;
figure(52); clf; 
  pcolor(rlat,pavgLAY(1:97,3000),smoothn(moo(:,1:97)',1)); shading interp; colormap(usa2); set(gca,'ydir','reverse')
  ylim([100 1000]); caxis([-1 +1]*0.015); colorbar('horizontal'); colormap(cmap); title(['if RH0 were 100 \newline Zonal MAX dWVfrac/dt expected from dT/dt'])

figure(8); figure(28); figure(29); figure(30); 
aslmap(6,rlat65,rlon73,smoothn((reshape(maskLF.*results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap(31,rlat65,rlon73,smoothn((reshape(maskLF.*mmwPert - maskLF.*mmw0,72,64)') ,1), [-90 +90],[-180 +180]); title('dmmw/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap_polar(32,rlat65,rlon73,smoothn((reshape(maskLF.*results(:,6)',72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap_polar(33,rlat65,rlon73,smoothn((reshape(maskLF.*mmwPert - maskLF.*mmw0,72,64)') ,1), [-90 +90],[-180 +180]); title('dmmw/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see How closely do changes in surface and column water vapor follow Clausius-Clapeyron scaling in climate-change simulations?
%% P A O’Gorman, C J Muller, https://core.ac.uk/download/pdf/4426849.pdf, or saved in PDF/change_of_RH_with_stemp_GCM_PGorman.pdf
%% P A O'Gorman and C J Muller 2010 Environ. Res. Lett. 5 025207 DOI 10.1088/1748-9326/5/2/025207

booU = mmwPert - mmw0; 
if iAK > 0
  booE = era5.trend_mmw;
else
  booE = boo;
end

figure(36); clf
booUU = booU./mmw0 * 100;
booEE = booE./mmw0 * 100;

iPG = 2;
if iPG == 0
  %% do nothing
elseif iPG == 1
  booUU = booUU ./results(:,6)'; 
  booEE = booEE ./era5.trend_stemp; 
  booUU = nanmean(reshape(booUU,72,64),1);
  booEE = nanmean(reshape(booEE,72,64),1);
elseif iPG == 2
  zonal_umbc_st = nanmean(reshape(results(:,6),72,64),1);
  zonal_era5_st = nanmean(reshape(era5.trend_stemp,72,64),1);
  booUU = nanmean(reshape(booUU,72,64),1);
  booEE = nanmean(reshape(booEE,72,64),1);
  booUU = booUU ./ zonal_umbc_st;
  booEE = booEE ./ zonal_era5_st;
end

plot(rlat,booUU,rlat,booEE,'linewidth',2); 
plot(rlat,smooth(booUU,10),rlat,smooth(booEE,10),'linewidth',2);
if iPG == 0 | iPG == 1
  figure(36); ylim([-1 +1]*20)
elseif iPG == 2
  figure(36); ylim([0 +1]*20)
end
plotaxis2; ylabel('d% mmw/dST  %/K');                       xlabel('Latitude'); legend('UMBC','ERA5');

%%%%%%%%%%%%%%%%%%%%%%%%%

if std(pMean17years.landfrac) < eps
  [salti, landfrac] = usgs_deg10_dem(pMean17years.rlat, pMean17years.rlon);
  pMean17years.landfrac = landfrac;
else
  landfrac = pMean17years.landfrac;
end  
lfmaskA = ones(1,4608);
lfmaskL = (landfrac > 0);
lfmaskL = (landfrac == 1);
lfmaskO = (landfrac == 0);

figure(36); booUU = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,mmwPert-mmw0,    results(:,6)');
figure(36); booEE = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,mmwPertERA5-mmw0,era5.trend_stemp);
plot(rlat,smooth(booUU.all,10),'r',rlat,smooth(booEE.all,10),'r--',rlat,smooth(booUU.ocean,10),'b',rlat,smooth(booEE.ocean,10),'b--',rlat,smooth(booUU.land,10),'g',rlat,smooth(booEE.land,10),'g--','linewidth',2);
if iPG == 0 | iPG == 1
  figure(36); ylim([-1 +1]*20)
elseif iPG == 2
  figure(36); ylim([0 +1]*20)
end
plotaxis2; ylabel('d% mmw/dST  %/K');                       xlabel('Latitude'); legend('UMBC all','ERA5 all','UMBC ocean','ERA5 ocean','UMBC land','ERA5 land','location','best','fontsize',10);
  
%%%%%%%%%%%%%%%%%%%%%%%%%

clear junk;
junk.color = 'k'; junk.iNorS = +1; aslmap_polar(34,rlat65,rlon73,10*smoothn((reshape(booU,72,64)') ,1), [-90 +90],[-180 +180],junk); title('UMBC dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)
  junk.color = 'k'; junk.iNorS = +1; aslmap_polar(35,rlat65,rlon73,10*smoothn((reshape(booE,72,64)') ,1), [-90 +90],[-180 +180],junk); title('ERA5 dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)
junk.color = 'k'; junk.iNorS = -1; aslmap_polar(34,rlat65,rlon73,10*smoothn((reshape(booU,72,64)') ,1), [-90 +90],[-180 +180],junk); title('UMBC dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)
  junk.color = 'k'; junk.iNorS = -1; aslmap_polar(35,rlat65,rlon73,10*smoothn((reshape(booE,72,64)') ,1), [-90 +90],[-180 +180],junk); title('ERA5 dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)
aslmap(34,rlat65,rlon73,10*smoothn((reshape(booU,72,64)') ,1), [-90 +90],[-180 +180]); title('UMBC dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)
  aslmap(35,rlat65,rlon73,10*smoothn((reshape(booE,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 dmmw/dt mm/decade');  caxis([-1 +1]*2); colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
