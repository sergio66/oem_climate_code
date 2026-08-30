%% this is /home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m

%{
cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES
  %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
  set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m : iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1; 
  
  %% clear, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
  rad0 : 
    set_rtp.m : use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';
    rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p cpu2021 --array=1-4608 sergio_matlab_jobB.sbatch 6
    watch 'ls -lt /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof* | wc -l'

    when done
      mv JUNK/allkcbands_prof* JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Raw/.

  %% clear, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
  pert : 
    set_rtp.m : use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp';            
    rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p cpu2021 --array=1-4608 sergio_matlab_jobB.sbatch 6
    watch 'ls -lt /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof* | wc -l'

    when done
      mv JUNK/allkcbands_prof* JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Pert/.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/PLOTMISC
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlib/plotutils

%% set control params

%% a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 05;  %% use CarbonTracker CO2 trends
%% a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 15;  %% use CarbonTracker CO2 trends
%% a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2v2.mat';    iNumYears = 10;  %% use CarbonTracker CO2 trends
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 20;  %% use CarbonTracker CO2 trends ****

iMakePertProf = -1;
iRawOrPert = +1;
iRawOrPert = -1;
iReadSave_Or_MakePlots = +1;
iReadSave_Or_MakePlots = -1;

if iReadSave_Or_MakePlots > 0
  if iRawOrPert > 0
    disp('reading RAW calc')
  else
    disp('reading PERT calc')
  end
else
  disp('already have read in kCARTA calcs, made convolutions ... now look at results')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do we need to make perturbed (trends) rtp profile???

%%%
%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_rtp.m
%% use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';            %% clear, raw,  see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
%% use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp';       %% clear, pert, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m
%%%

if iMakePertProf > 0
  pert_kcarta_fluxes_for_paper
  rtpwrite(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_' num2str(2002 + iNumYears) '_monthlyERA5_pert.rp.rtp'],h,ha,pert,pa);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% do we need to read in and convolve kcarta runs? or just use the convolved results

if iReadSave_Or_MakePlots > 0
  iaFound = [];
  allrads.dall = zeros(2420000,4608);
  
  disp('reading 4608 files, progress marked by dots = 100, crosses = 1000 ')
  for ii = 1 : 4608
    if mod(ii,1000) == 0
      fprintf(1,'+')
    elseif mod(ii,100) == 0
      fprintf(1,'.')
    end
  
    fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof' num2str(ii) '.mat'];
    if iRawOrPert > 0
      fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Raw_' num2str(iNumYears,'%02d') 'years/allkcbands_prof' num2str(ii) '.mat'];
    elseif iRawOrPert < 0
      fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Pert_' num2str(iNumYears,'%02d') 'years/allkcbands_prof' num2str(ii) '.mat'];
    end
    if exist(fin)
      iaFound(ii) = +1;
      loader = ['a = load(''' fin ''');'];
      eval(loader);
      allrads.dall(:,ii) = a.dall;
    else
      iaFound(ii) = 0;
    end
  end
  fprintf(1,' done reading\n');
    
  allrads.wall = a.wall;
  
  moo = find(iaFound > 0);
  fprintf(1,'found %4i; max = %4i \n',sum(iaFound),max(moo))
  
  disp('convolving')

  moo = find(allrads.wall >= 605.00-1e-05);
  [fc,qc] = convolve_airs(allrads.wall(moo),allrads.dall(moo,:),1:2834); 
  moo = find(allrads.wall <= fc(1));
  [fgc,qgc] = quickconvolve(allrads.wall(moo),allrads.dall(moo,:),0.5,0.5); 
  
  oo = find(iaFound==1); oo = oo(end); 
  plot(a.wall,a.dall,[fgc; fc],[qgc(:,oo); qc(:,oo)])
  
  if oo >= 1000
    moo = [qgc; qc];
    plot([fgc; fc],[qgc(:,1000); qc(:,1000)],[fgc; fc],nanmean(moo,2))
  end

  comment = 'see /home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m'
  if iRawOrPert > 0
    fname_save_kcrads_FarIR_ThermalIR = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/radsraw_' num2str(iNumYears,'%02d') 'years.mat'];
  else
    fname_save_kcrads_FarIR_ThermalIR = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/radspert_' num2str(iNumYears,'%02d') 'years.mat'];
  end
  saver = ['save ' fname_save_kcrads_FarIR_ThermalIR ' fc qc fgc qgc comment'];
  eval(saver)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  read_fileMean17years
  h = hMean17years;
  p = pMean17years;

  a0 = load(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/radsraw_' num2str(iNumYears,'%02d') 'years.mat']);
  aF = load(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/radspert_' num2str(iNumYears,'%02d') 'years.mat']);

  loader = ['load h2645structure.mat'];
  eval(loader)

  %% see driver_contrast_obsrates_05_10_15_20.m
  dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/';
  if iNumYears == 05  
    aNumYears = load([dir0 'iType_10_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  elseif iNumYears == 10
    aNumYears = load([dir0 'iType_11_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  elseif iNumYears == 15
    aNumYears = load([dir0 'iType_12_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  elseif iNumYears == 20
    aNumYears = load([dir0 'iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat']);
  end

  %% now see /home/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/check_new_AIRS2645.m, airs_convolve_file_numchans.m for 2834 --> 2645
  fc   = [a0.fgc; a0.fc(h.ichan)];
  rad0 = [a0.qgc; a0.qc(h.ichan,:)];
  radF = [aF.qgc; aF.qc(h.ichan,:)];

  cosrlat = cos(p.rlat * pi/180);
  cosFIR_ThIR = ones(length(fc),1) * cosrlat;
  cosThIR     = ones(h.nchan,1) * cosrlat;

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(1); clf; plot(fc,nanmean(radF-rad0,2)); ylabel('\delta Rad')

  %% see eg ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/clust_tile_fits_quantiles_loop72lons.m  which calls
  %% see eg ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/tile_fits_quantiles.m for radiance trend ... d.rad_quantile_desc(iT,iQ,iV)   iT=1:457 timesteps iQ = 1:5 quants, iV = 1:2645 chans
  %% 
  %%    %% Strows stuff (run from his dir)
  %%    fdirpre      = '/home/strow/Work/Airs/Tiles/Data/Quantv1';        %% symbolic link to /home/strow/Work/Airs/Tiles/Data/Quantv1 -> /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon
  %%    fdirpre_out  = '/home/strow/Work/Airs/Tiles/Data/Quantv1_fits';
  %%    
  %%    %% Sergio stuff (run from my dir)
  %%    fdirpre      = '../DATAObsStats_StartSept2002_CORRECT_LatLon/';   %% symbolic link to ./DATAObsStats_StartSept2002_CORRECT_LatLon -> /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon
  %%    fdirpre_out  = '../DATAObsStats_StartSept2002_CORRECT_LatLon/';
  %% 
  %%    fdirpre_out  = '/home/strow/Work/Airs/Tiles/Data/Quantv1_fits';
  %%    if iQAX == 1
  %%      fn_summary = sprintf('LatBin%1$02d/LonBin%2$02d/summarystats_LatBin%1$02d_LonBin%2$02d_timesetps_001_%3$03d_V1.mat',lati,loni,i16daysSteps);
  %%    elseif iQAX == 3
  %%      fn_summary = sprintf('LatBin%1$02d/LonBin%2$02d/iQAX_3_summarystats_LatBin%1$02d_LonBin%2$02d_timesetps_001_%3$03d_V1.mat',lati,loni,i16daysSteps);
  %%    end
  %%    fn_summary = fullfile(fdirpre,fn_summary);  
  %% 
  %%    fnout = ['LatBin' num2str(lati,'%02d') '/LonBin' num2str(loni,'%02d') '/iQAX_3_fits_LonBin' num2str(loni,'%02d') '_LatBin' num2str(lati,'%02d') '_V1_'];
  %%    fnout = [fnout    num2str(startdate,'%04d') '_' num2str(stopdate,'%04d')  '_TimeStepsX' num2str(i16daysStepsX,'%03d')];
  %%    fnout = fullfile(fdirpre_out,fnout);
  %%
  %% so eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin32/LonBin36/iQAX_3_summarystats_LatBin32_LonBin36_timesetps_001_457_V1.mat 
  %%    has  rad_quantile_desc(iT,iQ,iV) ie the time series
  %% while eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin32/LonBin36/iQAX_3_fits_LonBin36_LatBin32_V1_TimeSteps457.mat
  %%    has  b_desc(iV,iQ,1:10) which are the radiance fits (1=const,2=trend,3-10 are sin/cos)
  %%         dbt_desc(iV,iQ)    which are the BT       trends
  %% doo = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin32/LonBin36/iQAX_3_fits_LonBin36_LatBin32_V1_TimeSteps457.mat');
  %% figure(1); clf; plot(1:2645,doo.b_desc(:,3,2),1:2645,doo.dbt_desc(:,3)); hl = legend('rad','bt'); ylabel('dX/dt')

  get_or_make_radiance_trends

  moo = direct_rad_trend; moo = nanmean(moo,2);
  figure(1); clf; plot(fc,nanmean(radF-rad0,2),h.vchan,moo); ylabel('\delta Rad')

  moo = direct_rad_trend;  moo = nanmean(moo.*cosThIR,2);
  figure(2); clf; plot(fc,nanmean((radF-rad0).*cosFIR_ThIR,2),h.vchan,moo); ylabel('\delta Rad coswgt')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(3); clf; plot(fc,nanmean(rad2bt(fc,radF)-rad2bt(fc,rad0),2));  ylabel('\delta BT')

  moo = (reshape(aNumYears.b_desc,4608,2645))'; moo = nanmean(moo,2);
  moo = (reshape(aNumYears.b_desc,4608,2645))'; moo = nanmean(moo.*cosThIR,2);

  figure(4); clf; plot(fc,nanmean(cosFIR_ThIR.*(rad2bt(fc,radF)-rad2bt(fc,rad0)),2),h.vchan,moo);  ylabel('\delta BT coswgt')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  gray = [1 1 1]*0.5;
  
  figure(5); clf
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

  tafov(1) = nexttile;
  plot(fc,nanmean(radF,2),'linewidth',2); hold on
  plot(fc,ttorad(fc,284),'color',gray); hold off
  ylabel('Radiance \newline units');
  xlim([0 2800]); plotaxis2;
  shade2(5,0,0,650,150,gray,0.3);

  tafov(2) = nexttile;
  moo = direct_rad_trend; moo = nanmean(moo,2);
  plot(fc,nanmean(radF-rad0,2),'b',h.vchan,moo,'r'); hold on;
  plot(fc,ttorad(fc,284+0.02)-ttorad(fc,284),'color',gray); hold off
  ylabel('d(Radiance)/dt \newline [yr-1]')
  xlim([0 2800]); ylim([-0.1 +0.051])
  xlabel('Wavenumber [cm-1]'); plotaxis2; 
  shade2(5,0,-0.1,650,0.15,gray,0.3);
  hl = legend('Simulated','Observed','location','best','fontsize',10);
  fairs = h.vchan;
  trendairs = moo;

  akleidon.fairs     = fairs;
  akleidon.trendairs = moo;
  akleidon.kcarta_chan       = fc;
  akleidon.kcartaradFinal    = nanmean(radF,2);
  akleidon.kcartarad284      = ttorad(fc,284);
  akleidon.kcartaradInit     = nanmean(rad0,2);  
  akleidon.kcartarad284delta = ttorad(fc,284+0.02);

  figure(6); plot(akleidon.kcarta_chan,akleidon.kcartaradFinal-akleidon.kcartaradInit,'b',akleidon.fairs,akleidon.trendairs,'r',akleidon.kcarta_chan,akleidon.kcartarad284delta-akleidon.kcartarad284,'k')
  fid = fopen('akleidon_airs.txt','w');
    fprintf(fid,'%8.2f %12.8f \n',akleidon.fairs,akleidon.trendairs);
    fclose(fid);
  figure(6); plot(akleidon.kcarta_chan,akleidon.kcartaradFinal,'b',akleidon.kcarta_chan,akleidon.kcartarad284delta,'r',akleidon.kcarta_chan,akleidon.kcartaradInit,'c',akleidon.kcarta_chan,akleidon.kcartarad284,'m')
  fid = fopen('akleidon_kcarta.txt','w');
    fprintf(fid,'%8.2f %12.8f %12.8f %12.8f %12.8f \n',akleidon.kcarta_chan,akleidon.kcartaradInit,akleidon.kcartaradFinal,akleidon.kcartarad284,akleidon.kcartarad284delta);
    fclose(fid);

  % Get rid of all extra space I can
  ta.Padding = 'none';
  ta.TileSpacing = 'none';
  ta.TileSpacing = 'compact';
  tafov(1).XTickLabel = '';

  %% aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/radiance_trends_20years.pdf')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
