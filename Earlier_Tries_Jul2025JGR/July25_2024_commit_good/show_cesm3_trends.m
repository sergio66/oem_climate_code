addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/maps/

load llsmap5

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

[Y,X] = meshgrid(rlat,rlon);
do_XX_YY_from_X_Y

figure(17)
aslmap(17,rlat65,rlon73,smoothn((reshape(moo.thestats64x72.stemprate,72,64)') ,1), [-90 +90],[-180 +180]); title('CESM3 19 years dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap(17,rlat65,rlon73,smoothn(moo.thestats64x72.stemprate' ,1), [-90 +90],[-180 +180]); title('CESM3 19 years dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5)

figure(18); 
boo = nanmean(reshape(moo.thestats64x72.stemprate,72,64),1);
plot(rlat,boo); plotaxis2; title('zonal skt trend')

figure(19); clf; pcolor(reshape(cos(YY*pi/180),72,64)'); colormap jet; colorbar
plot(YY,cos(YY*pi/180))

zoo = moo.thestats64x72.stemprate; zoo = zoo(:); scatter_coast(XX,YY,50,zoo)

cosY    = cos(Y*pi/180);
coslat  = cos(YY*pi/180);
indSST  = moo.thestats64x72.stemprate;    

indSST = moo.thestats64x72.stemprate; 
aha = indSST.*cosY; 
  fprintf(1,'mean CESM3 SKT trend = %8.6f \n',sum(aha(:))/sum(cosY(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iNumYears = 20;
%  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CESM3_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr','stemptrend');
%  cesm3_spectral_olr = junk.climcapsL3_spectral_olr;
%  cesm3_spectral_olr.stemptrend = junk.stemptrend.cesm;
%{
addpath /asl/matlib/plotutils

dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(15); aslprint([dir0 'cesm19year_Trates.pdf']);
figure(16); aslprint([dir0 'cesm19year_WVrates.pdf']);
figure(17); aslprint([dir0 'cesm19year_sktrates.pdf']);
%}
