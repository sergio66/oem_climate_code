fprintf(1,'set_driver_jacfile.m : anomalies!!! latbin = %2i timestep = %3i \n',driver.iLat,driver.i16daytimestep)
junk = num2str(driver.i16daytimestep,'%03d');

if iXJac == 1
  %% sarta time vary jacs
  driver.jacobian.filename = [];
  fprintf(1,'iXJac == 1 reading in anomaly timestep sarta jac file %s \n',driver.jacobian.filename)

elseif iXJac == 2  
  rlat = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat'); 
  rlat65 = rlat.latB2; rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
  rlon73 = (1:73); rlon73 = -180 + (rlon73-1)*5;  rlon = (1:72); rlon = -177.5 + (rlon-1)*5;
  [Y,X] = meshgrid(rlat,rlon);
  %[salti, landfrac] = usgs_deg10_dem(Y,X);    
  %lf = landfrac; lf = lf(:);
  
  do_XX_YY_from_X_Y
  [mmjunk,nnjunk] = size(XX);
  if mmjunk == 1
    %% need these reshaped for the pcolor plots below
    XX = XX';
    YY = YY';
  end

  %% kcarta time vary jac
  junk = load(driver.rateset.datafile);
  if strfind(driver.anomalyinfo.datafile,'_tile_')
    junk = load(driver.anomalyinfo.datafile,'LatBin','LonBin');
    %driver.iLon = junk.LonBin;
    junk = junk.LatBin;
  else
    %junk = junk.usethese{driver.anomalyinfo.latbin};
    junk = junk.usethese{driver.anomalyinfo.fatprofilebin};

    %YYmean = nanmean(YY(junk));
    %junk = find(rlat65 >= YYmean,1) - 1;
    junk = round(nanmean(junk));
  end
  driver.jacobian.filename = [];
  AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/';
  AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/';

  AHA = [AHA '/kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(junk,'%02i') '.mat']; %% ERA5,  2002-2022 20 year <avg cld = Q09> and NOT Q05
  driver.jacobian.filename = AHA;
  fprintf(1,'iXJac == 2 reading in timestep kcarta jac file %s \n',driver.jacobian.filename)
  iVersJac = 2021;   %% ERA5 CLR  from 2002-2021

elseif iXJac == 0
  %% constant kcarta jacs
  driver.jacobian.filename = [];
  fprintf(1,'iXJac == 0 reading in constant kcarta jac file %s \n',driver.jacobian.filename)
end

topts.jacobian.filename = driver.jacobian.filename;
driver.jacobian.filename = driver.jacobian.filename;

iVersJac = iXJac;
topts.iVersJac = iVersJac;
