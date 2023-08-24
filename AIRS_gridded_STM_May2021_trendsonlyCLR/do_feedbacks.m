disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')

disp('WARNING : remember my "compute_olr" code simply uses SARTA and multiplies by pi ie has BAND GAPS, is NOT doing FIR calcs, iit is NOT true OLR calc')
disp('          you really need RRTM or ecRad for that')
disp('          so in eg compute_feedbacks_X_ecRad, X = umbc,era5,airsl3,cmip6 .. the "spectral_olr" calc is bunk/rubbish/nonsense')
disp('  ')
disp('WARNING : remember my "compute_olr" code simply uses SARTA and multiplies by pi ie has BAND GAPS, is NOT doing FIR calcs, iit is NOT true OLR calc')
disp('          you really need RRTM or ecRad for that')
disp('          so in eg compute_feedbacks_X_ecRad, X = umbc,era5,airsl3,cmip6 .. the "spectral_olr" calc is bunk/rubbish/nonsense')
disp('  ')
disp('WARNING : remember my "compute_olr" code simply uses SARTA and multiplies by pi ie has BAND GAPS, is NOT doing FIR calcs, iit is NOT true OLR calc')
disp('          you really need RRTM or ecRad for that')
disp('          so in eg compute_feedbacks_X_ecRad, X = umbc,era5,airsl3,cmip6 .. the "spectral_olr" calc is bunk/rubbish/nonsense')
disp('  ')

if ~exist('umbc_spectral_olr')
  compute_feedbacks_umbc_ecRad   ; pause(0.1)
end

fprintf(1,'models string = %s \n',strMODELS);

if ~exist('airsL3_spectral_olr')
  junkx = input('load in airsL3, era5, cmip6 flux calcs from earlier (-1/+1 default) : ? ');
  if length(junkx) == 0
    junkx = 1;
  end
  if junkx > 0
    disp('loading in flux calcs from earlier');
    %junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','airsL3_spectral_olr'); airsL3_spectral_olr = junk.airsL3_spectral_olr;
    %junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','era5_spectral_olr'); era5_spectral_olr = junk.era5_spectral_olr;
    %junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5','cmip6_spectral_olr'); cmip6_spectral_olr = junk.cmip6_spectral_olr;

    %savename1 = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';      %% oops this is wrong
    %savename2 = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat';  %% oops this is wrong

    savename1 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_AIRSL3_ERA5_CMIP6_save.mat';     %% this is base lambda for AIRSL3_ERA5_CMIP6
    savename2 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat'; %% this is base lambda for CLIMCAPS_MERRA2_AMIP6

    if ~exist('iJorC')
      iJorC = +1;  %% do Joel L3
    end
    if ~exist('iEorM')
      iEorM = +5;  %% do ERA5
    end
    if ~exist('iAorC')
      iAorC = -1;  %% do CMIP6
    end

    iLoadJorC = -1;
    iLoadEorM = -1;
    iLoadAorC = -1;

    if iJorC == +1
      iLoadJorC = +1;
      junk = load(savename1,'airsL3_spectral_olr'); airsL3_spectral_olr = junk.airsL3_spectral_olr;
    end
    if iEorM == +5
      iLoadEorM = +1;
      junk = load(savename1,'era5_spectral_olr');   era5_spectral_olr = junk.era5_spectral_olr;
    end
    if iAorC == -1
      iLoadAorC = +1;
      junk = load(savename1,'cmip6_spectral_olr');  cmip6_spectral_olr = junk.cmip6_spectral_olr;
    end

    if iJorC == -1
      iLoadJorC = +1;
      junk = load(savename2,'airsL3_spectral_olr'); airsL3_spectral_olr = junk.airsL3_spectral_olr;
    end
    if iEorM == +2
      iLoadEorM = +1;
      junk = load(savename2,'era5_spectral_olr');   era5_spectral_olr = junk.era5_spectral_olr;
    end
    if iAorC == +1
      iLoadAorC = +1;
      junk = load(savename2,'cmip6_spectral_olr');  cmip6_spectral_olr = junk.cmip6_spectral_olr;
    end

  end
end

if ~exist('airsL3_spectral_olr')
  disp('WARNING : computing feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
  junk = input('(+1/default) Go ahead with the calcs, you have run make_profile_spectral_trends or (-1) oops, re-run it here : ');
  if length(junk) == 0
    junk == 1;
  end
  if junk < 0
    clear nwp_spectral_trends_cmip6_era5_airsL3_umbc
    nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends(cmip6,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
  end
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

%% make sure you use _ecRad because that uses all 10-3000 cm-1, rather than _spectral_olr which only uses SARTA (645-2830 cm-1)
do_avg_feedback1     %% crude attempt at zonal avg
do_avg_feedback1cos  %% crude attempt at zonal avg with cosine(lat) wgt
do_avg_feedback2     %% better attempt at zonal avg
do_avg_feedback2cos  %% better attempt at zonal avg with cosine(rlat) wgt  BEST

figure(4); 
%%% Ryan suggested normalizing using dERASST for all, instead of the individual dXSST X=ERA or CMIP6 or UMBC or AIRSL3 
%%% iERAnorm = input('Do you wish to redo the feedback by using only dERA SKT instead of individual d SKT? (-1/default) no (+1) yes : ');
%%% if length(iERAnorm) == 0
%%%   iERAnorm = -1;
%%% end
%%% if iERAnorm > 0
%%%   redo_feedbacks_dERA5ST_dt
%%%   do_avg_feedback2cos_dERA5ST_dt  %% better attempt at zonal avg with cosine(rlat) wgt  BEST
%%% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ifield(umbc_spectral_olr,'allperts_no_tracegas_ecRad')
  quick_regress_OLR_SST_for_feedbacks
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = input('save the OLR feedbacks??? (-1/+1) : ');
if junk > 0
  %{
  %%%% this is how I made  olr_feedbacks_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  feedbackname = ['olr_feedbacks_' strMODELS '_save.mat'];
  saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays   airsL3_spectral_olr era5_spectral_olr cmip6_spectral_olr airsL3 era5 cmip6'];
  %%%% this is how I made  olr_feedbacks_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  %}

  iSaveAll = input('save ERA5, CMIP,airsL3 as well? (-1 default/+1) : ');
  if length(iSaveAll) == 0
    iSaveAll = -1;
  end
  if iSaveAll == -1
    feedbackname = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELS '_numyears_' num2str(iNumYears) '.mat'];
    saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];
  else
    feedbackname = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELS '_ALL_numyears_' num2str(iNumYears) '.mat'];
    saver = ['save ' feedbackname ' era5_spectral_olr cmip6_spectral_olr airsL3_spectral_olr umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];
  end
  fprintf(1,'saving to %s \n',feedbackname);
  eval(saver);
end

