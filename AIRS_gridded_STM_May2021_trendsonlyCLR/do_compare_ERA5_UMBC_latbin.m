if ~exist('phmm')
  [hmm,~,phmm,~] = rtpread('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/2012/FixedNAN/all4608_era5_full12months_Qcumulative09.rtp');
end

ind = 1:72:4608; ind = ind + 36; %% this is GMT line from -90 S to + 90 N
ind = (1:72) + (iCompare-1)*72;         %% this is one latitude bin eg S. Pole or Equator
  
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

%moo = load(['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(iCompare,'%02d') '.mat']);
moo = load(['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(iCompare,'%02d') '.mat']);

gray = 0.8;
figure(21); plot(rlon,smooth(rates(1520,ind),3),'gx-',rlon,smooth(results(ind,6),3),'b',rlon,smooth(era5.trend_stemp(ind),3),'r','linewidth',2); 
hold on; plot(rlon,pMean17years.landfrac(ind)/10,'color',[gray gray gray],'linewidth',2); hold off
  plotaxis2; legend('d(BT1231)/dt','UMBC dSurfT/dt','ERA5 dSurfT/dt','landfrac/10','location','best','fontsize',10);
  fprintf(1,'midpoimt has spres = %10.4f and nlevs = %3i \n',phmm.spres(ind(36)),phmm.nlevs(ind(36)))

figure(22); 
moo = squeeze(moo.jacT(:,:,36)); pcolor(f(1:520),plays(1:100),moo(:,1:520)); shading interp; set(gca,'ydir','reverse'); title('T jac for lonbin = 36'); line([650 800],[phmm.spres(ind(36)) phmm.spres(ind(36))],'color','k','linewidth',2);
jett = jet; jett(1,:) = 1; colorbar; colormap(jett);

%figure(20); pcolor(rlat,pavg,resultsT(ind,:)');  xlabel('latitude'); title('GMT line T');  caxis([-0.1 +0.1]);        
%figure(21); pcolor(rlat,pavg,resultsO3(ind,:)'); xlabel('latitude'); title('GMT line O3'); caxis([-0.01 +0.01]*0.25); 
%for ii = 19:21; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

figure(23); 
  plot(f,nanmean(rates(:,ind),2),'b',f,nanmean(fits(:,ind),2),'r',f,nanmean(rates(:,ind),2)-nanmean(fits(:,ind),2),'k'); plotaxis2;
  title(['Latbin ' num2str(iCompare,'%02d') ' SPECTRA dBT/dt K/yr']); 
  hl = legend('input','fit','diff','location','best'); xlim([650 1650])

iCompare = input('Enter latbin over which to compare ERA5 vs UMBC trends (1:64, -1 to stop) : ');
