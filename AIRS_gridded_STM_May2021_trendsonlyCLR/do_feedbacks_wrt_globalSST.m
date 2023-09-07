% me thinks these are incorrect because have to perturb by global SST
% redo_feedbacks_dERA5ST_dt
% do_avg_feedback2cos_dERA5ST_dt  %% better attempt at zonal avg with cosine(rlat) wgt  BEST

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

disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')
disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')
disp('each model takes about 15 minutes for SARTA and ecRad to run completely ie 15 min for UMBC, 45 min for ERA5/AIRSL3/CMIP6, 45 min for MERRA2/CLIMCAPSL3/AMIP6')
disp('so much better if you can read in pre-computed ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6')
disp('  ')

clear feedbackname*
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
  
if ~exist('umbc_spectral_olr') 
  if exist(feedbacknameUMBC)
    iUMBCexist = +1;
    fprintf(1,'for %2i years the UMBC OLR ecRad file %s exists \n',iNumYears,feedbacknameUMBC)
    lser = ['!ls -lt ' feedbacknameUMBC];
    eval(lser)
    junk = input('re-do the OLR feedbacks??? (-1 [default] /+1) : ');    
    if junk > 0
      umbc_spectral_olr = struct;    %% so it has no fields
      umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,'UMBC');
      %%% compute_feedbacks_umbc_ecRad   ; pause(0.1)
    else
      fprintf(1,'  reading in %s \n',feedbacknameUMBC)
      loader = ['load ' feedbacknameUMBC];
      eval(loader);
    end
  else
    umbc_spectral_olr = struct;    %% so it has no fields
    umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,'UMBC');
    %%% compute_feedbacks_umbc_ecRad   ; pause(0.1)
  end
end

junk = input('save the OLR feedbacks??? (-1 [default] /+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk > 0
  %{
  %%%% this is how I made  olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  if iSwap_ERA_2012_08_15 < 0
    feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];
  else
    feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_swap_profile.mat'];
  end
  if iRaw_or_Unc == -1
    feedbacknameUMBC = feedbacknameUMBC(1:end-4);
    feedbacknameUMBC = [feedbacknameUMBC '_unc_factor' num2str(maxratio,'%0.2f') '.mat'];
  end
  saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays   airsL3_spectral_olr era5_spectral_olr cmip6_spectral_olr airsL3 era5 cmip6'];
  saver = ['save ' feedbackname ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays   airsL3_spectral_olr era5_spectral_olr cmip6_spectral_olr airsL3 era5 cmip6'];
  %%%% this is how I made  olr_feedbacks_globalSST_AIRSL3_ERA5_CMIP6_save.mat olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_save.mat %%%%
  %}
  
  if junk > 0
    junk2 = +1;
    if exist(feedbacknameUMBC)
      lser = ['!ls -lth ' feedbacknameUMBC];
      eval(lser);
      junk2 = input('file already exists, overwrite (-1 default/+1) : ');
      if length(junk2) == 0
        junk2 = -1;
      end
    end       
    saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
    if junk2 > 0
      fprintf(1,'saving to %s \n',feedbacknameUMBC);
      eval(saver);
    else
      fprintf(1,'this already exists %s not saving \n',feedbacknameUMBC);
    end
  end
end

fprintf(1,'models string = %s \n',strMODELS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FRESH COMPUTE  NWP/AIRSL3/XMIP6 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

iFreshComputeNWP_L3_feedbacks = input('(-1/default) read in old files or (+1) freshly brew compute ERA5/AIRSL3/CMIP6 or MERRA2/CLIMCAPSL3/AMIP6 feedbacks : ');
if length(iFreshComputeNWP_L3_feedbacks) == 0
  iFreshComputeNWP_L3_feedbacks = -1;
end
if iFreshComputeNWP_L3_feedbacks < 0
  fprintf(1,'  loading in %s \n',feedbacknameNWP_ERA5);
  loader = ['load ' feedbacknameNWP_ERA5];
  eval(loader)
  if exist('strMODELSX')
    fprintf(1,'  loading in %s \n',feedbacknameNWP_MERRA2);
    loader = ['load ' feedbacknameNWP_MERRA2];
    eval(loader)
  end
elseif iFreshComputeNWP_L3_feedbacks > 0
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
  
    % compute_feedbacks_airsL3_ecRad  ; pause(0.1)
    % compute_feedbacks_era5_ecRad    ; pause(0.1)
    % compute_feedbacks_cmip6_ecRad   ; pause(0.1)
  
    era5_spectral_olr = struct;    %% so it has no fields
    era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,-1,'ERA5');
    
    aL3trend.stemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
    aL3trend.ptemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
    aL3trend.gas_1 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1;
    aL3trend.gas_3 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3;
    airsL3_spectral_olr = struct;    %% so it has no fields
    airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,-1,'AIRS L3');
    
    c6trend.stemp = cmip6.trend_stemp;
    c6trend.ptemp = cmip6.trend_ptemp;
    c6trend.gas_1 = cmip6.trend_gas_1;
    c6trend.gas_3 = cmip6.trend_gas_3;
    cmip6_spectral_olr = struct;    %% so it has no fields
    cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,-1,'CMIP6');
  
    junk = input('save ERA3/AIRSL3/CMIP6 only (-1 [default]/+1) : ');
    if length(junk) == 0
      junk = -1;
    end
    if junk > 0
      junk2 = +1;
      if exist(feedbacknameNWP_ERA5)
        lser = ['!ls -lth ' feedbacknameNWP_ERA5];
        eval(lser);
        junk2 = input('file already exists, overwrite (-1 default/+1) : ');
        if length(junk2) == 0
          junk2 = -1;
        end
      end       
      stemptrend.era5   = era5.trend_stemp;
      stemptrend.airsL3 = aL3trend.stemp;
      stemptrend.cmip6  = c6trend.stemp;
      saver = ['save ' feedbacknameNWP_ERA5 ' era5_spectral_olr cmip6_spectral_olr airsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
      if junk2 > 0
        fprintf(1,'saving to %s \n',feedbacknameNWP_ERA5);
        eval(saver);
      else
        fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_ERA5);
      end
    end  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~exist('climcapsL3_spectral_olr') & exist('amip6')
    disp('WARNING : computing feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
    junk = input('(+1/default) Go ahead with the calcs, you have run make_profile_spectral_trends or (-1) oops, re-run it here : ');
    if length(junk) == 0
      junk == 1;
    end
    if junk < 0
      clear nwp_spectral_trends_amip6_merra2_climcapsL3_umbc
      nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = make_profile_spectral_trends(amip6,merra2,climcapsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
    end
  
    % compute_feedbacks_airsL3_ecRad  ; pause(0.1)
    % compute_feedbacks_era5_ecRad    ; pause(0.1)
    % compute_feedbacks_cmip6_ecRad   ; pause(0.1)
  
    merra2_spectral_olr = struct;    %% so it has no fields
    merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,-1,'MERRA2');

    if ~exist('nwp_spectral_trends_amip6_merra2_climcapsL3_umbc')
     nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = make_profile_spectral_trends(amip6,merra2,climcapsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
    end
    cL3trend.stemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.stemp;
    cL3trend.ptemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.ptemp;
    cL3trend.gas_1 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_1;
    cL3trend.gas_3 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_3;
    climcapsL3_spectral_olr = struct;    %% so it has no fields
    climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,cL3trend.ptemp,cL3trend.gas_1,cL3trend.gas_3,climcapsL3_spectral_olr,-1,'CLIMCAPS L3');
    
    a6trend.stemp = amip6.trend_stemp;
    a6trend.ptemp = amip6.trend_ptemp;
    a6trend.gas_1 = amip6.trend_gas_1;
    a6trend.gas_3 = amip6.trend_gas_3;
    amip6_spectral_olr = struct;    %% so it has no fields
    amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,a6trend.ptemp,a6trend.gas_1,a6trend.gas_3,amip6_spectral_olr,-1,'AMIP6');
  
    junk = input('save MERRA2/CLIMCAPSL3/AMIP6 only (-1 [default]/+1) : ');
    if length(junk) == 0
      junk = -1;
    end
    if junk > 0
      junk2 = +1;
      if exist(feedbacknameNWP_MERRA2)
        lser = ['!ls -lth ' feedbacknameNWP_MERRA2];
        eval(lser);
        junk2 = input('file already exists, overwrite (-1 default/+1) : ');
        if length(junk2) == 0
          junk2 = -1;
        end
      end       
      stemptrend.era5   = merra2.trend_stemp;
      stemptrend.airsL3 = cL3trend.stemp;
      stemptrend.cmip6  = a6trend.stemp;
      saver = ['save ' feedbacknameNWP_MERRA2 ' merra2_spectral_olr amip6_spectral_olr climcapsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
      if junk2 > 0
        fprintf(1,'saving to %s \n',feedbacknameNWP_MERRA2);
        eval(saver);
      else
        fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_MERRA2);
      end
    end
  
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LOAD IN OLDER NWP/AIRSL3/XMIP6 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('airsL3_spectral_olr') & iFreshComputeNWP_L3_feedbacks < 0
  junkx = input('need to load in airsL3, era5, cmip6 flux calcs from earlier (-1/+1 default) : ? ');
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

disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
disp('now you can clear the memory and run driver_olr_fluxchanges.m to make comparisons against CERES trends')
