figure(3); clf
figure(4); clf

disp('REMEMBER THESE ARE NOT ZONAL AVG PLOTS but one lat bin, meridional!!!!')
disp('REMEMBER THESE ARE NOT ZONAL AVG PLOTS but one lat bin, meridional!!!!')
disp('REMEMBER THESE ARE NOT ZONAL AVG PLOTS but one lat bin, meridional!!!!')

if iType ~= 4
  figure(3); pcolor(rlon,plevsnwp,thesave.oz2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
    xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
    title('OZ frac rates straight from \newline CMIP6 or ERA5 zonal levels')

  figure(4); pcolor(rlon,plevsx,thesave.oz2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
    xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('OZ rates from rtp after klayers')
end

figure(7); pcolor(rlon,plevsnwp,thesave.t2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  title('T rates straight from \newline CMIP6 or ERA5 zonal levels')

figure(8); pcolor(rlon,plevsx,thesave.t2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

figure(9); pcolor(rlon,plevsnwp,thesave.wv2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  title('WV frac rates straight from \newline CMIP6 or ERA5 zonal levels')

figure(10); pcolor(rlon,plevsx,thesave.wv2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

figure(11); pcolor(rlon,plevsnwp,thesave.rh2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); 
  title('RH rates straight from \newline CMIP6 or ERA5 zonal levels')

figure(12); pcolor(rlon,plevsx,thesave.rh2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates from rtp after klayers')

if (YMEnd(1) - YMStart(1)) <= 5
  figure(11); caxis([-1 +1]*0.25*2*10)
  figure(12); caxis([-1 +1]*0.25*2*10)

  figure(7); caxis([-1 +1]*0.15*2*5)
  figure(8); caxis([-1 +1]*0.15*2*5)

  figure(9);  caxis([-1 +1]*0.15/10*2*5)
  figure(10); caxis([-1 +1]*0.15/10*2*5)

  figure(2); ylim([-1 +1]*0.2)

end
