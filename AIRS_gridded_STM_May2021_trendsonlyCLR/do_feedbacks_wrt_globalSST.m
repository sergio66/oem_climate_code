% me thinks these are incorrect because have to perturb by global SST
% redo_feedbacks_dERA5ST_dt
% do_avg_feedback2cos_dERA5ST_dt  %% better attempt at zonal avg with cosine(rlat) wgt  BEST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  ')
disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!!')
disp('module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')

disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')
disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')
disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')
disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('  ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear feedbackname*
clear *spectral_olr
iLambda_UseGlobalSST = +1;

iUMBCexist = -1;
if iSwap_ERA_2012_08_15 < 0
  feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];
else
  feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_swap_profile.mat'];
end
if iRaw_or_Unc == -1
  feedbacknameUMBC = feedbacknameUMBC(1:end-4);
  feedbacknameUMBC = [feedbacknameUMBC '_unc_factor' num2str(maxratio,'%0.2f') '.mat'];
end
if iAllorSeasonal == -1
  feedbacknameUMBC = feedbacknameUMBC(1:end-4);
  feedbacknameUMBC = [feedbacknameUMBC '_DJF.mat'];
elseif iAllorSeasonal == -2
  feedbacknameUMBC = feedbacknameUMBC(1:end-4);
  feedbacknameUMBC = [feedbacknameUMBC '_MAM.mat'];
elseif iAllorSeasonal == -3
  feedbacknameUMBC = feedbacknameUMBC(1:end-4);
  feedbacknameUMBC = [feedbacknameUMBC '_JJA.mat'];
elseif iAllorSeasonal == -4
  feedbacknameUMBC = feedbacknameUMBC(1:end-4);
  feedbacknameUMBC = [feedbacknameUMBC '_SON.mat'];
end

do_compute_save_UMBC_feedbacks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FRESH COMPUTE  NWP/AIRSL3/XMIP6 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'models string = %s \n',strMODELS);

if iSwap_ERA_2012_08_15 < 0
  feedbacknameNWP_ERA5 = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELS '_numyears_' num2str(iNumYears,'%02d') '.mat'];
else
  feedbacknameNWP_ERA5 = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELS '_numyears_' num2str(iNumYears,'%02d') '_swap_profile.mat'];
end
if iRaw_or_Unc == -1
  feedbacknameNWP_ERA5 = feedbacknameNWP_ERA5(1:end-4);
  feedbacknameNWP_ERA5 = [feedbacknameNWP_ERA5 '_unc_factor' num2str(maxratio,'%0.2f') '.mat'];
end
if exist(feedbacknameNWP_ERA5)
  lser = ['!ls -lt ' feedbacknameNWP_ERA5];
  eval(lser)
  fprintf(1,'for %2i years ERA/AIRSL3/CMIP6 file %s already exists \n',iNumYears,feedbacknameNWP_ERA5)
end

if exist('strMODELSX')
  if iSwap_ERA_2012_08_15 < 0
    feedbacknameNWP_MERRA2 = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELSX '_numyears_' num2str(iNumYears,'%02d') '.mat'];
  else
    feedbacknameNWP_MERRA2 = ['/asl/s1/sergio/JUNK/olr_feedbacks_' strMODELSX '_numyears_' num2str(iNumYears,'%02d') '_swap_profile.mat'];
  end
  if iRaw_or_Unc == -1
    feedbacknameNWP_MERRA2 = feedbacknameNWP_MERRA2(1:end-4);
    feedbacknameNWP_MERRA2 = [feedbacknameNWP_MERRA2 '_unc_factor' num2str(maxratio,'%0.2f') '.mat'];
  end
  if exist(feedbacknameNWP_MERRA2)
    lser = ['!ls -lt ' feedbacknameNWP_MERRA2];
    eval(lser)
    fprintf(1,'for %2i years MERRA2/CLIMCAPSL3/AMIP6 file %s already exists \n',iNumYears,feedbacknameNWP_MERRA2)
  end
else
  disp('hmm, guess you did not read in MERRA2/CLIMCAPSL3/AMIP6 trends ....')
end

iFreshComputeNWP_L3_feedbacks = input('(-1 [default]) read in old files or (+1) freshly brew compute one or both of ERA5/AIRSL3/CMIP6 or MERRA2/CLIMCAPSL3/AMIP6 feedbacks : ');
if length(iFreshComputeNWP_L3_feedbacks) == 0
  iFreshComputeNWP_L3_feedbacks = -1;
  iXFreshComputeNWP_L3_feedbacks = -10;
end

if iFreshComputeNWP_L3_feedbacks > 0
  iXFreshComputeNWP_L3_feedbacks = input(' freshly brew compute (0) ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6 feedbacks (+1) ERA5/AIRSL3/CMIP6 only (+2) MERRA2/CLIMCAPSL3/AMIP6 only : ');
end

if iFreshComputeNWP_L3_feedbacks < 0
  fprintf(1,'  loading in %s \n',feedbacknameNWP_ERA5);
  clear stemptrend
  loader = ['load ' feedbacknameNWP_ERA5];
  eval(loader)
  airsL3_spectral_olr.stemptrend = stemptrend.airsL3;
  cmip6_spectral_olr.stemptrend  = stemptrend.cmip6;
  era5_spectral_olr.stemptrend   = stemptrend.era5;
  clear stemptrend
  if exist('strMODELSX')
    clear stemptrend
    fprintf(1,'  loading in %s \n',feedbacknameNWP_MERRA2);
    loader = ['load ' feedbacknameNWP_MERRA2];
    eval(loader)
    climcapsL3_spectral_olr.stemptrend = stemptrend.airsL3;
    amip6_spectral_olr.stemptrend  = stemptrend.cmip6;
    merra2_spectral_olr.stemptrend   = stemptrend.era5;    
    clear stemptrend
  end
elseif iFreshComputeNWP_L3_feedbacks > 0 & iXFreshComputeNWP_L3_feedbacks == 2
  fprintf(1,'  loading in %s \n',feedbacknameNWP_ERA5);
  loader = ['load ' feedbacknameNWP_ERA5];
  eval(loader)
elseif iFreshComputeNWP_L3_feedbacks > 0 & iXFreshComputeNWP_L3_feedbacks == 1 & exist('strMODELSX')
  fprintf(1,'  loading in %s \n',feedbacknameNWP_MERRA2);
  loader = ['load ' feedbacknameNWP_MERRA2];
  eval(loader)
end

if iFreshComputeNWP_L3_feedbacks >= 0
  if ~exist('airsL3_spectral_olr')
    do_compute_save_AIRSL3_ERA5_CMIP6_feedbacks
  end
  if ~exist('climcapsL3_spectral_olr') & exist('amip6')
    do_compute_save_CLIMCAPSL3_MERRA2_AMIP6_feedbacks
  end  
end

whos *spectral_olr
iOhOh = -1;
if ~isfield(umbc_spectral_olr,'stemptrend')
  iOhOh = +1;
  disp('warning .. umbc_spectral_olr does not have field stemptrend')
end
if ~isfield(era5_spectral_olr,'stemptrend')
  iOhOh = +1;
  disp('warning .. era5_spectral_olr does not have field stemptrend')
end
if iOhOh > 0
  error('someone does not have field stemptrend')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LOAD IN OLDER NWP/AIRSL3/XMIP6 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('airsL3_spectral_olr') & iFreshComputeNWP_L3_feedbacks < 0
  junkx = input('need to load in airsL3, era5, cmip6 flux calcs from earlier (-1/+1 [default]) : ? ');
  if length(junkx) == 0
    junkx = 1;
  end
  if junkx > 0

    disp('loading in flux calcs from earlier');
    savename1 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat';     %% this is base lambda for AIRSL3_ERA5_CMIP6, globalSST
    savename2 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_globalSST_CLIMCAPS_MERRA2_AMIP6_save.mat'; %% this is base lambda for CLIMCAPS_MERRA2_AMIP6, globalSST

    savename1 = ['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_'   num2str(iNumYears,'%02d') '.mat'];
    savename2 = ['/asl/s1/sergio/JUNK/olr_feedbacks_CLIML3_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'];
    fprintf(1,'AIRSL3_ERA5_CMIP6   : will use %s \n',savename1);
    fprintf(1,'CLIML3_MERRA2_AMIP6 : will use %s \n',savename2);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%

    iJorC = +1;  %% do Joel L3
    iEorM = +5;  %% do ERA5
    iAorC = -1;  %% do CMIP6

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

    %%%%%%%%%%%%%%%%%%%%%%%%%

    iJorC = -1;  %% do Barnet L3
    iEorM = +2;  %% do MERRA2
    iAorC = +1;  %% do AMIP6

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns0 = 500;
ns0 = 50;
ns0 = 1;
wonk = umbc_spectral_olr.feedback_ecRad.planck.individual; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(75,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('UMBC \lambda_{Planck}')
wonk = umbc_spectral_olr.feedback_ecRad.wv.individual; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(76,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('UMBC \lambda_{WV}')

wonk = airsL3_spectral_olr.feedback_ecRad.planck.individual; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(77,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('AIRSL3 \lambda_{Planck}')
wonk = airsL3_spectral_olr.feedback_ecRad.wv.individual; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(78,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('AIRSL3 \lambda_{WV}')

wonk = era5_spectral_olr.feedback_ecRad.planck.individual; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(79,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('ERA5 \lambda_{Planck}')
wonk = era5_spectral_olr.feedback_ecRad.wv.individual; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(80,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('ERA5 \lambda_{WV}')

wonk = cmip6_spectral_olr.feedback_ecRad.planck.individual; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(81,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('CMIP6 \lambda_{Planck}')
wonk = cmip6_spectral_olr.feedback_ecRad.wv.individual; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
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

show_avg_feedbacks_plethora

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('now you can eg run "show_olr_ecRad_feedback" ')
disp('now you can eg run "show_olr_ecRad_feedback" ')
disp('now you can eg run "show_olr_ecRad_feedback" ')

%disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
%disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
%disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
