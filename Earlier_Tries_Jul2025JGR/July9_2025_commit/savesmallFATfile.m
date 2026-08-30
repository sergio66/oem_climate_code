if ~exist('strMODELS')
  strMODELS = 'NoMODELS';
end
disp('suggested final save names ...')
if topts.set_era5_cmip6_airsL3 == 5
  start_apriori_str = 'ERA5';
elseif topts.set_era5_cmip6_airsL3 == 2
  start_apriori_str = 'MERRA2';
elseif topts.set_era5_cmip6_airsL3 == 1
  start_apriori_str = 'ERA-I';
elseif topts.set_era5_cmip6_airsL3 == 6
  start_apriori_str = 'CMIP6';
elseif topts.set_era5_cmip6_airsL3 == -6
  start_apriori_str = 'AMIP6';
elseif topts.set_era5_cmip6_airsL3 == 3
  start_apriori_str = 'AIRSL3';
elseif topts.set_era5_cmip6_airsL3 == -3
  start_apriori_str = 'CLIMCAPS';
elseif topts.set_era5_cmip6_airsL3 == 8
  start_apriori_str = 'MLSL3';
else
  start_apriori_str = '0';
end

if topts.ocb_set == 1
  start_apriori_str = ['_cal_' start_apriori_str];
elseif topts.ocb_set == 2
  start_apriori_str = ['_bias_' start_apriori_str];
end

numretlayers_str = [num2str(topts.iNlays_retrieve) 'fatlayers'];

if length(iaSequential) == 1 & iaSequential(1) == -1
  cGorS = 'GULP';
else
  cGorS = 'SEQN';
end

genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_'         numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX1   %s \n',genericoutname);
genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_uncX3_'   numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX3   %s \n',genericoutname);
genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_uncX100_' numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX100 %s \n',genericoutname);
genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_'         numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX1   %s \n',genericoutname);
genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_uncX3_'   numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX3   %s \n',genericoutname);
genericoutname = ['/asl/s1/sergio/JUNK/smallgather_tileCLRnight_' cGorS '_dataset' num2str(dataset,'%02d') '_Q' num2str(iQuantile,'%02d') '_newERA5_2021jacs_startwith' start_apriori_str '_uncX100_' numretlayers_str '_' strMODELS '_feedback.mat']; 
  fprintf(1,'suggested name uncX100 %s \n',genericoutname);

junk = input('save savesmallFATfile???? (-1 default/+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk == 1
  junk2 = +1;
  genericoutname = input('Enter name of savesmallFATfile : ');
  if exist(genericoutname)
    lser = ['!ls -lth ' genericoutname];
    eval(lser);
    junk2 = input('file already exists, overwrite (-1 default/+1) : ');
    if length(junk2) == 0
      junk2 = -1;
    end
  end

  if junk2 > 0
    saver = ['save -v7.3 ' genericoutname];
    saver = [saver ' deltaRH deltaRHlat deltaRHunc results resultsunc resultsO3 resultsO3unc resultsT resultsTunc resultsWV resultsWVunc '];
    saver = [saver ' topts fracWV fracWVunc fracO3 fracO3unc deltaT deltaTunc save_cov_set chisqrR pert_unc'];
    saver = [saver ' rates fits nedt'];
    saver = [saver ' pavg iNumLay xb* pert pjunk20 maskLF rlat rlat65 rlon rlon73 nlays_straight_from_results '];
    if iAK > 0
      figure(29); figure(30); 
      figure(45); figure(46); 
      deltaRH = real(deltaRH);
      figure(54); clf; waha = squeeze(nanmean(reshape(deltaRH,100,72,64),2));        pcolor(waha); shading interp; colorbar; set(gca,'ydi','reverse'); title('UMBC dRH/dt'); colormap(llsmap5); caxis([-1 +1]*0.5)
      figure(55); clf; waha = squeeze(nanmean(reshape(era5.trend_RH,100,72,64),2));  pcolor(waha); shading interp; colorbar; set(gca,'ydi','reverse'); title('ERA5 dRH/dt'); colormap(llsmap5); caxis([-1 +1]*0.5)
      figure(56); clf; waha = save_cov_set.xb_wvz(2,:);                              aslmap(56,rlat65,rlon73,smoothn(reshape(waha,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.015); title('xb(WVfrac(gnd))')
      era5_stemprate = era5.trend_stemp;
      era5_wvrate    = era5.trend_gas_1;
      era5_rhrate    = era5.trend_RH;
      era5_ozrate    = era5.trend_gas_3;
      era5_tzrate    = era5.trend_ptemp;
      era5_plays     = era5.trend_plays;
      saver = [saver ' *_ak*_era5* era5_plays era5_*rate pjunk20 h p'];
    end
    eval(saver);
    
    quick_compare_era5_retrieval

  end 
end
