moo = load('/asl/s1/sergio/CESM3/cesm3_64x72_rates_stats_Sept2002_Aug2021_19yr_desc.mat');

figure(15); pcolor(meanvaluebin(moo.rlat),moo.Tlevs/100,squeeze(nanmean(moo.thestats64x72.ptemprate,1))'); colorbar; caxis([-1 +1]*0.15); 
            shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('CESM3 19 years dT/dt'); colormap(llsmap5)
figure(16); pcolor(meanvaluebin(moo.rlat),moo.Qlevs/100,squeeze(nanmean(moo.thestats64x72.waterrate,1))'); colorbar; caxis([-1 +1]*0.015); 
            shading interp; set(gca,'ydir','reverse'); ylim([100 1000]); title('CESM3 19 years dWVfrac/dt'); colormap(llsmap5)

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

figure(17)
aslmap(17,rlat65,rlon73,smoothn((reshape(moo.thestats64x72.stemprate,72,64)') ,1), [-90 +90],[-180 +180]); title('CESM3 19 years dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
