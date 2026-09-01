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

file_trend_airsobs = ['direct_rad_trends_Q03_' num2str(iNumYears,'%02d') 'years.mat'];
if ~exist(file_trend_airsobs)
  iCnt = 0;
  for JJ = 1 : 64
    fprintf(1,'reading rad and bt trends from latbin %2i of 64 \n',II)
    for II = 1 : 72
      iCnt = iCnt + 1;
      finx = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/'];
      finx = [finx 'DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' num2str(JJ,'%02i') '/LonBin' num2str(II,'%02i') '/iQAX_3_fits_LonBin' num2str(II,'%02i') '_LatBin' num2str(JJ,'%02i') '_V1_TimeSteps457.mat'];
      junk = load(finx,'b_desc','dbt_desc');
      direct_dbt_trend(:,iCnt) = junk.dbt_desc(:,3);
      direct_rad_trend(:,iCnt) = squeeze(junk.b_desc(:,3,2));
    end
  end
  comment = 'see driver_read_kcarta_fluxes_for_paper.m and get_or_make_radiance_trends.m';
  saver = ['save ' file_trend_airsobs ' direct_dbt_trend direct_rad_trend comment'];
  eval(saver);  
else
  loader = ['load ' file_trend_airsobs];
  eval(loader);
end
