% https://www.nature.com/articles/s43247-024-01620-3#data-availability
% Strong persistent cooling of the stratosphere after the Hunga eruption
% Matthias Stocker, Andrea K. Steiner, Florian Ladstädter, Ulrich Foelsche & William J. Randel 
% Communications Earth & Environment volume 5, Article number: 450 (2024) 

iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualT(:,i19km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*4); 
  title('T(t,lat) at 19 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualT(:,i27km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*4); 
  title('T(t,lat) at 27 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualT(:,i32km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*4); 
  title('T(t,lat) at 32 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualWV(:,i19km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 19 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualWV(:,i27km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 27 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),squeeze(annualWV(:,i32km,:))); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 32 km'); xlim([2021.75 2024]); ylim([-1 +1]*50); line([timeHT timeHT],[-1 +1]*50,'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Latitude [deg]');

iFig = iFig + 1; figure(iFig); clf; junk = squeeze(nanmean(annualT(tropicsHT,:,:),1)); pcolor(iaYears,havg,junk); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*4);
  title('T(t,z) -30S to +10 N'); xlim([2021.75 2024]); ylim([16 32]); line([timeHT timeHT],[16 32],'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Height [km]')
iFig = iFig + 1; figure(iFig); clf; junk = squeeze(nanmean(annualWV(tropicsHT,:,:),1)); pcolor(iaYears,havg,junk); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25);
  title('WV(t,z) -30S to +10 N'); xlim([2021.75 2024]); ylim([16 32]); line([timeHT timeHT],[16 32],'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Height [km]')
  axis([2020 2025 0 25])

iFig = iFig + 1; figure(iFig); clf; junk = squeeze(nanmean(annualT(tropicsHT2,:,:),1)); pcolor(iaYears,havg,junk); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*4);
  title('T(t,z) -30S to +30 N'); xlim([2021.75 2024]); ylim([16 32]); line([timeHT timeHT],[16 32],'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Height [km]')
iFig = iFig + 1; figure(iFig); clf; junk = squeeze(nanmean(annualWV(tropicsHT2,:,:),1)); pcolor(iaYears,havg,junk); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25);
  title('WV(t,z) -30S to +30 N'); xlim([2021.75 2024]); ylim([16 32]); line([timeHT timeHT],[16 32],'color','k','linewidth',2); 
  xlabel('Time'); ylabel('Height [km]')
  axis([2020 2025 0 25])

%%%%%%%%%%%%%%%%%%%%%%%%%

make_T_WV_O3_ST_trends_from_anoms

return

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

smn = 1;
smn = 0.1;
iFig = iFig - 8;
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualT(:,i19km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*2); 
  title('T(t,lat) at 19 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualT(:,i27km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*2); 
  title('T(t,lat) at 27 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualT(:,i32km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*2); 
  title('T(t,lat) at 32 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualWV(:,i19km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 19 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualWV(:,i27km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 27 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
iFig = iFig + 1; figure(iFig); clf; pcolor(iaYears,rlat(2:end),smoothn(squeeze(annualWV(:,i32km,:)),smn)); shading interp;  colorbar; colormap(cmap); caxis([-1 +1]*0.25); 
  title('WV(t,lat) at 32 km'); xlim([2021.75 2024]); ylim([-1 +1]*50)
