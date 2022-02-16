%{
compute_feedbacks_umbc   ; pause(0.1)
compute_feedbacks_airsL3 ; pause(0.1)
compute_feedbacks_era5   ; pause(0.1)
compute_feedbacks_cmip6  ; pause(0.1)
figure(75); scatter_coast(p.rlon,p.rlat,50,era5_spectral_olr.feedback.planck-airsL3_spectral_olr.feedback.planck); caxis([-1 +1]/100);  colormap(usa2);  title('ERA5-AIRSL3 \lambda_{Planck}')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')

if ~exist('umbc_spectral_olr')
  compute_feedbacks_umbc_ecRad   ; pause(0.1)
end

if ~exist('airsL3_spectral_olr')
  junkx = input('load in airsL3, era5, cmip6 flux calcs from earlier (-1/+1 default) : ? ');
  if length(junk) == 0
    junkx = 1;
  end
  if junkx > 0
    disp('loading in flux calcs from earlier');
    junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','airsL3_spectral_olr'); airsL3_spectral_olr = junk.airsL3_spectral_olr;
    junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','era5_spectral_olr'); era5_spectral_olr = junk.era5_spectral_olr;
    junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','cmip6_spectral_olr'); cmip6_spectral_olr = junk.cmip6_spectral_olr;
  end
end

if ~exist('airsL3_spectral_olr')
  compute_feedbacks_airsL3_ecRad ; pause(0.1)
  compute_feedbacks_era5_ecRad   ; pause(0.1)
  compute_feedbacks_cmip6_ecRad  ; pause(0.1)
end

ns0 = 500;
ns0 = 50;
ns0 = 1;
wonk = umbc_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(75,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('UMBC \lambda_{Planck}')
wonk = umbc_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(76,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('UMBC \lambda_{WV}')

wonk = airsL3_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(77,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('AIRSL3 \lambda_{Planck}')
wonk = airsL3_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(78,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('AIRSL3 \lambda_{WV}')

wonk = era5_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(79,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('ERA5 \lambda_{Planck}')
wonk = era5_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(80,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('ERA5 \lambda_{WV}')

wonk = cmip6_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(81,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('CMIP6 \lambda_{Planck}')
wonk = cmip6_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(82,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('CMIP6 \lambda_{WV}')

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(79); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(80); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(81); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(82); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(79); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(80); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(81); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(82); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_avg_feedback1     %% crude attempt at zonal avg
do_avg_feedback1cos  %% crude attempt at zonal avg with cosine(lat) wgt
do_avg_feedback2     %% better attempt at zonal avg
do_avg_feedback2cos  %% better attempt at zonal avg with cosine(rlat) wgt  BEST

%%% Ryan suggested normalizing using dERASST for all, instead of the individual dXSST X=ERA or CMIP6 or UMBC or AIRSL3 
iERAnorm = input('Do you wish to redo the feedback by using only dERA SKT instead of individual d SKT? (-1) no, default (+1) yes : ');
if iERAnorm > 0
  redo_feedbacks_dERA5ST_dt
  do_avg_feedback2cos_dERA5ST_dt  %% better attempt at zonal avg with cosine(rlat) wgt  BEST
end
