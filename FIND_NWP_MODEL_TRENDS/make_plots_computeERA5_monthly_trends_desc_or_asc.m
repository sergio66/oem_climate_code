if exist('llsmap5.mat')
  colorstr = 'llsmap5';
  load llsmap5
else
  colorstr = 'usa2';
end
colorstr = 'usa2';

figure(1); clf; scatter_coast(pall.rlon,pall.rlat,40,trend_stemp); title('ERA5 trend  stemp K/yr');    caxis([-1 +1]*0.15); colormap(colorstr);
figure(2); clf; scatter_coast(pall.rlon,pall.rlat,40,trend_RHSurf); title('ERA5 trend  UGH RHsurf pc/yr'); caxis([-1 +1]*0.4); colormap(colorstr);
aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 dST/dt');      caxis([-1 +1]*0.15); colormap(colorstr)
aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 UGH dRHSurf/dt'); caxis([-1 +1]*0.25); colormap(colorstr)

figure(3); clf; junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend ptemp K/yr');     caxis([-1 +1]*0.15); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log');    shading interp; ylim([10 1000]); colorbar
figure(4); clf; junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); shading interp; ylim([100 1000]); colorbar
figure(5); clf; junk = reshape(trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend WVfrac /yr');     caxis([-1 +1]*0.01); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); shading interp; ylim([100 1000]); colorbar

figure(6); clf; junk = squeeze(nanmean(pall.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean ptemp K');    caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(7);clf; junk = squeeze(nanmean(pall.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean RH percent'); caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%{
figure(8); clf; junk = squeeze(nanstd(pall.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');    caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(9); clf; junk = squeeze(nanstd(pall.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent'); caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(10); clf; junk = squeeze(max(pall.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');    caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(11); clf; junk = squeeze(max(pall.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent'); caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); clf; junk = reshape(trend_nwp_ptemp,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend ptemp K/yr');  caxis([-1 +1]*0.15); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(9); clf; junk = reshape(trend_nwp_rh,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(10); clf; junk = reshape(trend_nwp_gg,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend SH g/g/yr');  caxis([0 +2.5]*1e-5); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(11); clf; junk = reshape(trend_nwp_ppmv,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend PPMV /yr');  caxis([0 40]); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); shading interp; ylim([100 1000]); colorbar
figure(12); clf; junk = reshape(trend_nwp_frac,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  frac /yr');  caxis([-10 +10]*1e-3); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); shading interp; ylim([100 1000]); colorbar

pause(0.1);

junk1 = squeeze(nanmean(pall.nwp_gas_1,1)); junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2)); 
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(13); clf; pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  SH g/g/yr VERS2');  caxis([0 +2.5]*1e-5); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = toppmv(pall.nwp_plevs,pall.nwp_ptemp,pall.nwp_gas_1,18,21); junk1 = squeeze(nanmean(junk1,1));
junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2));
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(14); clf; pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1);
figure(14); clf; loglog(nanmean(junk1,2),trend_nwp_plevs_mean);  set(gca,'ydir','reverse'); xlim([1 1e4]); grid
figure(14); clf; pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  PPMV/yr VERS2');  caxis([0 40]); colormap(colorstr); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% this is MLS comparisons, Frank Werner JPL
figure(11); ylim([1 300]); colormap(colorstr)
figure(07); ylim([1 300]);
figure(09); ylim([1 300]); caxis([-1 1]*1e-7); colormap(colorstr)
figure(12); ylim([1 300]); caxis([-1 1]*1e-7); colormap(colorstr)
figure(10); ylim([1 300]); caxis([-1 1]*1e-1); colormap(colorstr)
figure(13); ylim([1 300]); caxis([-1 1]*1e-1); colormap(colorstr)
%}
