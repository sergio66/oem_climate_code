ind = 1:72:4608; ind = ind + 31; %% this is GMT line from -90 S to + 90 N
ind = (1:72) + (iCompare-1)*72;         %% this is S. Pole
  
fprintf(1,'  LatBin %2i = %8.6f iiBin = [%4i,%4i] \n',iCompare,rlat(iCompare),ind(1),ind(end));

figure(19); subplot(311); pcolor(rlon,pjunk20,resultsWV(ind,:)');          shading interp; ylabel('UMBC'); title(['Latbin ' num2str(iCompare,'%02d') ' dWVfrac/dt /yr']); caxis([-1e-2 +1e-2]);      
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
            subplot(312); pcolor(rlon,pjunk20,waterrate_akF_era5(ind,:)'); shading interp; ylabel('ERA5*AK'); caxis([-1e-2 +1e-2]);      
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
            subplot(313); pcolor(rlon,pjunk20,waterrate_ak0_era5(ind,:)'); shading interp; ylabel('ERA5'); xlabel('longitude'); caxis([-1e-2 +1e-2]);      
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])

figure(20); subplot(311); pcolor(rlon,pjunk20,resultsT(ind,:)');          shading interp; ylabel('UMBC'); title(['Latbin ' num2str(iCompare,'%02d') ' dT/dt K/yr']); caxis([-1 +1]*0.1)
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
            subplot(312); pcolor(rlon,pjunk20,temprate_akF_era5(ind,:)'); shading interp; ylabel('ERA5*AK'); caxis([-1 +1]*0.1)
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
            subplot(313); pcolor(rlon,pjunk20,temprate_ak0_era5(ind,:)'); shading interp; ylabel('ERA5'); xlabel('longitude'); caxis([-1 +1]*0.1)
              colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])

gray = 0.8;
figure(21); plot(rlon,smooth(rates(1520,ind),3),'gx-',rlon,smooth(results(ind,6),3),'b',rlon,smooth(era5.trend_stemp(ind),3),'r','linewidth',2); 
hold on; plot(rlon,pMean17years.landfrac(ind)/10,'color',[gray gray gray],'linewidth',2); hold off
  plotaxis2; legend('d(BT1231)/dt','UMBC dSurfT/dt','ERA5 dSurfT/dt','landfrac/10','location','best','fontsize',10);
  
%figure(20); pcolor(rlat,pavg,resultsT(ind,:)');  xlabel('latitude'); title('GMT line T');  caxis([-0.1 +0.1]);        
%figure(21); pcolor(rlat,pavg,resultsO3(ind,:)'); xlabel('latitude'); title('GMT line O3'); caxis([-0.01 +0.01]*0.25); 
%for ii = 19:21; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

iCompare = input('Enter latbin over which to compare ERA5 vs UMBC trends (1:64, -1 to stop) : ');
