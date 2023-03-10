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
    saver = [saver ' deltaRH deltaRHlat deltaRHunc results resultsO3 resultsO3unc resultsT resultsTunc resultsWV resultsWVunc rlat rlat65 rlon rlon73 topts fracO3 fracWV fracWVunc deltaT deltaTunc'];
    if iAK > 0
      figure(54); waha = squeeze(nanmean(reshape(deltaRH,100,72,64),2));        pcolor(waha); shading interp; colorbar; set(gca,'ydi','reverse'); title('UMBC dRH/dt'); colormap(llsmap5); caxis([-1 +1]*0.5)
      figure(55); waha = squeeze(nanmean(reshape(era5.trend_RH,100,72,64),2));  pcolor(waha); shading interp; colorbar; set(gca,'ydi','reverse'); title('ERA5 dRH/dt'); colormap(llsmap5); caxis([-1 +1]*0.5)

      era5_stemprate = era5.trend_stemp;
      saver = [saver ' *_ak*_era5* era5_stemprate pjunk20'];
    end
    eval(saver);
  end 
end
