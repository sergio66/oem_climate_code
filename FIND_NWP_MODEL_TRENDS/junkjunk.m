
figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.stemp,1)); colormap(jet); title('ERA5 mean stemp')
figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.RHSurf,1)); colormap(jet); title('ERA5 mean RHsurf DO NOT BELIEVE')
figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.TwSurf,1)); colormap(jet); title('ERA5 mean TWSurf DO BELIEVE')
figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.mmw,1)); colormap(jet); title('ERA5 mean mmw')

figure(5); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.stemp); colormap(jet); title('ERA5 mean stemp')
figure(6); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.RHSurf); colormap(jet); title('ERA5 mean RHsurf DO NOT BELIEVE')
figure(7); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.TwSurf); colormap(jet); title('ERA5 mean TWSurf DO BELIEVE')
figure(8); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.mmw); colormap(jet); title('ERA5 mean mmw')

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = flipud(pN./pD);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

figure(9); clf; junk = reshape(a.pnew_op.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
  title('Mean T')
figure(10); clf; junk = reshape(a.pnew_op.RH,100,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([00 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
  title('Mean RH')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

if iTrendsOrAnoms > 0
  computeERA5_surface_trends

  figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');        caxis([-0.2 +0.2]); colormap(usa2);
  figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf UGH pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
  figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf UGH K/yr');  caxis([-0.2 +0.2]); colormap(usa2);
  figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');      caxis([-0.2 +0.2]); colormap(usa2);

  figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');        caxis([-1 +1]*0.15); colormap(llsmap5);
  figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf UGH pc/yr'); caxis([-1 +1]*0.4);  colormap(llsmap5);
  figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf UGH K/yr');  caxis([-1 +1]*0.1);  colormap(llsmap5);
  figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');      caxis([-1 +1]*0.2);  colormap(llsmap5);

  figure(1); clf; aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]);, title('ERA5 trend  stemp K/yr');        caxis([-1 +1]*0.15); colormap(llsmap5);
  figure(2); clf; aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 trend  RHsurf UGH pc/yr');  caxis([-1 +1]*0.4);  colormap(llsmap5);
  figure(3); clf; aslmap(3,rlat65,rlon73,smoothn((reshape(trend_TwSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 trend  TWSurf UGH K/yr');   caxis([-1 +1]*0.1);  colormap(llsmap5);
  figure(4); clf; aslmap(4,rlat65,rlon73,smoothn((reshape(trend_mmw,72,64)') ,1), [-90 +90],[-180 +180]);    title('ERA5 trend  colwater mm/yr');    caxis([-1 +1]*0.2);  colormap(llsmap5);
  pause(0.1)

  if ~exist(fout_trendjunk_surface)
    fprintf(1,'saving surface/scalar trend file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk_surface)
    saver = ['save ' fout_trendjunk_surface ' comment trend*'];
    eval(saver);
  end

elseif iTrendsOrAnoms < 0
  computeERA5_surface_anoms
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iTrendsOrAnoms == 1 
  computeERA5_atmos_trends

  if isfield(all,'nwp_plevs')
    trend_nwp_plevs_mean = mean(squeeze(nanmean(all.nwp_plevs,1)),2);
  end
  
  trend_plays = flipud(pN./pD);
  
  trend_rlat = all.rlat;
  trend_rlon = all.rlon;
  trend_rlat64 = rlat; trend_rlon72 = rlon;
  %trend_plevs37 = permute(all.nwp_plevs,[2 1 3]); trend_plevs37 = reshape(trend_plevs37,37,227*4608); trend_plevs37 = mean(trend_plevs37,2);
  
  %find_computeERA5_monthly_trends_foutname
  if ~exist(fout_trendjunk)
    fprintf(1,'saving trend file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk)
    saver = ['save ' fout_trendjunk ' comment trend*'];
    eval(saver);
  end
  
  figure(5); clf; pcolor_sin(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_ptemp,100,72,64),2))); caxis([-1 +1]*0.150); set(gca,'ydir','reverse'); colormap(llsmap5); title('ERA5 dT/dt'); shading interp
    set(gca,'yscale','log'); ylim([10 1000]); colorbar('horizontal')
  figure(6); clf; pcolor_sin(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_gas_1,100,72,64),2))); caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colormap(llsmap5); title('ERA5 dWVfrac/dt'); shading interp
    set(gca,'yscale','linear'); ylim([100 1000]); colorbar('horizontal')

elseif iTrendsOrAnoms == -1 
  computeERA5_atmos_anoms

  if isfield(all,'nwp_plevs')
    anom_nwp_plevs_mean = mean(squeeze(nanmean(all.nwp_plevs,1)),2);
  end
  
  anom_plays = flipud(pN./pD);
  
  anom_rlat = all.rlat;
  anom_rlon = all.rlon;
  anom_rlat64 = rlat; anom_rlon72 = rlon;
  %anom_plevs37 = permute(all.nwp_plevs,[2 1 3]); anom_plevs37 = reshape(anom_plevs37,37,227*4608); anom_plevs37 = mean(anom_plevs37,2);
  
  %find_computeERA5_monthly_anoms_foutname
  if ~exist(fout_anomjunk)
    fprintf(1,'saving anom file : can type in a separate window         watch "ls -lt %s " \n',fout_anomjunk)
    saver = ['save ' fout_anomjunk ' comment anom*'];
    eval(saver);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')

figure(1); clf; scatter_coast(all.rlon,all.rlat,40,trend_stemp); title('ERA5 trend  stemp K/yr');    caxis([-1 +1]*0.15); colormap(llsmap5);
figure(2); clf; scatter_coast(all.rlon,all.rlat,40,trend_RHSurf); title('ERA5 trend  UGH RHsurf pc/yr'); caxis([-1 +1]*0.4); colormap(llsmap5);
addpath /asl/matlib/maps/
aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 dST/dt');      caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 UGH dRHSurf/dt'); caxis([-1 +1]*0.25); colormap(llsmap5)

figure(3); junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend ptemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(5); junk = reshape(trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend WVfrac /yr');  caxis([-1 +1]*0.01); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(6); junk = squeeze(nanmean(all.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean ptemp K');  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(7); junk = squeeze(nanmean(all.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean RH percent');  caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%{
figure(8); junk = squeeze(nanstd(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(9); junk = squeeze(nanstd(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(10); junk = squeeze(max(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(11); junk = squeeze(max(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); junk = reshape(trend_nwp_ptemp,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend ptemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(9); junk = reshape(trend_nwp_rh,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(10); junk = reshape(trend_nwp_gg,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend SH g/g/yr');  caxis([0 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(11); junk = reshape(trend_nwp_ppmv,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend PPMV /yr');  caxis([0 40]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(12); junk = reshape(trend_nwp_frac,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  frac /yr');  caxis([-10 +10]*1e-3); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

pause(0.1);

junk1 = squeeze(nanmean(all.nwp_gas_1,1)); junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2)); 
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(13); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  SH g/g/yr VERS2');  caxis([0 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = toppmv(all.nwp_plevs,all.nwp_ptemp,all.nwp_gas_1,18,21); junk1 = squeeze(nanmean(junk1,1));
junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2));
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1);
figure(14); loglog(nanmean(junk1,2),trend_nwp_plevs_mean);  set(gca,'ydir','reverse'); xlim([1 1e4]); grid
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  PPMV/yr VERS2');  caxis([0 40]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% this is MLS comparisons, Frank Werner JPL
figure(11); ylim([1 300]); colormap(llsmap5)
figure(07); ylim([1 300]);
figure(09); ylim([1 300]); caxis([-1 1]*1e-7); colormap(llsmap5)
figure(12); ylim([1 300]); caxis([-1 1]*1e-7); colormap(llsmap5)
figure(10); ylim([1 300]); caxis([-1 1]*1e-1); colormap(llsmap5)
figure(13); ylim([1 300]); caxis([-1 1]*1e-1); colormap(llsmap5)
%}
