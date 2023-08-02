%% note this is a SCRIPT

if length(intersect(settings.dataset,[6 8])) == 1
  [co2x,n2ox,ch4x] = get_co2_n2o_ch4_for_strow_override0(driver,iVersJac); %% sets co2x,n2ox,ch4x using special CRIS or OCO years
else
  [co2x,n2ox,ch4x] = get_co2_n2o_ch4_for_strow_override(driver,settings,iVersJac); %% sets co2x,n2ox,ch4x using 2-19 years :: AIRS which starts in 2002/09
end

iAdjCo2 = +1;
co2x = co2x * iAdjCo2;

if iAdjCo2 >= 10
  boo = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2022_08.mat');
  JOBJOBJOB = (driver.iLat-1)*72 + driver.iLon;
  co2x = boo.trend(JOBJOBJOB);
  disp('set_CO2_CH4_N2O_ESRL.m : co2 --> co2 CarbonTracker')
elseif iAdjCo2 > -1 & iAdjCo2 < +1
  %co2x = co2x * 1.1  %%%% SERGIO PUT THIS JULY 23, 2023  makes WV at surface pretty good across all lats, dT/dt < dT(ERA5)/dt in tropical troposphere
  %co2x = co2x * 0.9  %%%% SERGIO PUT THIS JULY 23, 2023   makes dWVfrac/dt < 0 at southern latitudes!!!s  dT/dt > dT(ERA5)/dt in tropical troposphere
  %co2x = co2x * 0.95  %%%% SERGIO PUT THIS JULY 23, 2023   makes dWVfrac/dt < 0 at southern latitudes!!!!!
  %co2x = co2x * 1.05  %%%% SERGIO PUT THIS JULY 23, 2023   makes dWVfrac/dt < 0 at southern latitudes!!!!!
  co2x = co2x * (1 + iAdjCo2);
  %disp('set_CO2_CH4_N2O_ESRL.m : co2 --> co2 x 1.05')  
end
driver.co2adj_ESRL = iAdjCo2;

if settings.set_tracegas == +1 & driver.i16daytimestep < 0 & settings.ocb_set ~= 1
  fprintf(1,'setting constant rates for tracegas apriori : CO2 = %8.6f  N2O = %8.6f   CH4 = %8.6f \n',co2x,n2ox,ch4x)
  if settings.co2lays == 1
    xb(1) = co2x;  % Set CO2 apriori
    xb(2) = 1;
    xb(3) = ch4x;
    xb(4) = 0.0; %% clouds, so dunno value
    xb(5) = 0.0; %% clouds, so dunno value

    xb(1) = co2x * 1;    % Set CO2 apriori
    xb(2) = n2ox * 1;    % set N2O 
    xb(3) = ch4x * 1;    % set CH4
    xb(4) = 0.0; %% clouds, so dunno value
    xb(5) = 0.0; %% clouds, so dunno value

    %xb(1:5) = xb(1:5)*10;
    %[hhh,~,ppp,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2019_monthlyERA5.rp.rtp');
    [hhh,~,ppp,~] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019.rtp');
    JOBJOBJOB = (driver.iLat-1)*72 + driver.iLon;
    if ppp.landfrac(JOBJOBJOB) < eps
      %% ocean, so use Isaac Held blog ideas  https://www.gfdl.noaa.gov/blog_held/47-relative-humidity-over-the-oceans/
      dBT1231 = 0.0025/2;
      dBT1231 = driver.rateset.rates(1520); 
      dBT1231_WV = dBT1231 * 0.07; %% remember saturation vapor pressure changes at 0.07/K and we want dRH = 0 BUT THIS 0.07 is for 285 K !!!!!!!!!!!!

      %%% see my notes, BK 45
      Lo = 2.5e6;  %%% J/kg
      Rv = 461.52; %%% J/kg/K
      moo = exp(Lo/Rv * dBT1231/ppp.stemp(JOBJOBJOB)/ppp.stemp(JOBJOBJOB))-1;
      moo = Lo/Rv * dBT1231/ppp.stemp(JOBJOBJOB)/ppp.stemp(JOBJOBJOB);
      dBT1231_WV = moo; %% remember saturation vapor pressure changes at 0.07/K and we want dRH = 0 BUT THIS 0.07 is for 285 K !!!!!!!!!!!!

    else
      %% land, so whos knows
      dBT1231_WV = 0.0;
      dBT1231    = 0.0;
    end
    
    iAdjLowerAtmWVfrac = topts.iAdjLowerAtmWVfrac;
    iAdjLowerAtmWVfracX = 0.0;
 
    %% code before 10 July 2023
    %% if abs(iAdjLowerAtmWVfrac-1) < eps
    %%   %% if iAdjLowerAtmWVfrac == 1, iAdjLowerAtmWVfracX = 1 so use dBT1231_WV as is from above
    %%   iAdjLowerAtmWVfracX = 1.0;
    %% elseif iAdjLowerAtmWVfrac > 0
    %%   %% see Bk 46 : dWVfrac/dt =  1/RH dRH/dt + Lo/Rv 
    %%   iAdjLowerAtmWVfracX = (1 + iAdjLowerAtmWVfrac/0.8);
    %% end

    %% code after 10 July 2023
    if iAdjLowerAtmWVfrac > 0
      %% see Bk 46 : dWVfrac/dt =  1/RH dRH/dt + Lo/Rv 
      iAdjLowerAtmWVfracX = iAdjLowerAtmWVfrac;
    else
      iAdjLowerAtmWVfracX = 0.0;
    end

    fprintf(1,'d/dt BT1231 from rates = %8.6f K/year --> WVfractional change needed for constant RH over ocean %8.6f frac/yr with iAdjLowerAtmWVfracX = %8.6f\n',driver.rateset.rates(1520),dBT1231_WV,iAdjLowerAtmWVfracX);
    if iAdjLowerAtmWVfracX > eps
      iVers = 0;  %% before noon Mar 11,2023, only sets WV   LOWEST 2 layers
      iVers = 1;  %% after  noon Mar 11,2023,      sets WV   LOWEST 6 layers
      iVers = 2;  %% after       Mar 27,2023,      sets WV+T LOWEST 6 layers

      TfacAdjAtmosphericAmplification = 1/4; %% factor of 1/4 is pretty good *****
      TfacAdjAtmosphericAmplification = 1/5; %% factor of 1/5, nah not so good, dCOWV/dt << dcolWV/dt for ERA5
      TfacAdjAtmosphericAmplification = 1/3; %% factor of 1/3 is pretty good *****
      TfacAdjAtmosphericAmplification = 1/2; %% factor of 1/2 is probably too large, messes up dRH/dt by making that too big at surface

      topts.iRHdelta0adjVers = iVers;
      topts.TfacAdjAtmosphericAmplification = TfacAdjAtmosphericAmplification;

      fprintf(1,'d/dt BT1231 from rates = %8.6f K/year --> iVers,TfacAdjAtmosphericAmplification = %2i %8.6f \n',iVers,TfacAdjAtmosphericAmplification)

      if iVers == 0
        %% tested and  savesmallFATfile --> /asl/s1/sergio/JUNK/test7_guessstartWV_Vers0_march11_2023.mat, commit on Sat Mar 11 19:17:59 2023 -0500
        xb(6+length(driver.jacobian.water_i)-0) = dBT1231_WV * iAdjLowerAtmWVfracX;
        xb(6+length(driver.jacobian.water_i)-1) = dBT1231_WV * iAdjLowerAtmWVfracX;
      elseif iVers == 1
        %% tested and savesmallFATfile --> /asl/s1/sergio/JUNK/test7_guessstartWV_Vers1_march11_2023.mat, commits on Sun Mar 12 10:01:58 2023 -0400 and Sun Mar 12 09:51:40 2023 -0400
        xb(6+length(driver.jacobian.water_i)-0) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 0/6);
        xb(6+length(driver.jacobian.water_i)-1) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 1/6);
        xb(6+length(driver.jacobian.water_i)-2) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 2/6);
        xb(6+length(driver.jacobian.water_i)-3) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 3/6);
        xb(6+length(driver.jacobian.water_i)-4) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 4/6);
        xb(6+length(driver.jacobian.water_i)-5) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 5/6);      
        xb(6+length(driver.jacobian.water_i)-6) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 6/6);
      elseif iVers == 2
        %% tested and savesmallFATfile --> /asl/s1/sergio/JUNK/test7_guessstartWV_Vers1_march11_2023.mat, commits on Sun Mar 12 10:01:58 2023 -0400 and Sun Mar 12 09:51:40 2023 -0400
        xb(6+length(driver.jacobian.water_i)-0) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 0/6);
        xb(6+length(driver.jacobian.water_i)-1) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 1/6);
        xb(6+length(driver.jacobian.water_i)-2) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 2/6);
        xb(6+length(driver.jacobian.water_i)-3) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 3/6);
        xb(6+length(driver.jacobian.water_i)-4) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 4/6);
        xb(6+length(driver.jacobian.water_i)-5) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 5/6);      
        xb(6+length(driver.jacobian.water_i)-6) = dBT1231_WV * iAdjLowerAtmWVfracX * (1 - 6/6);
       
        xb(6+length(driver.jacobian.temp_i)-0) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 0/6);
        xb(6+length(driver.jacobian.temp_i)-1) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 1/6);
        xb(6+length(driver.jacobian.temp_i)-2) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 2/6);
        xb(6+length(driver.jacobian.temp_i)-3) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 3/6);
        xb(6+length(driver.jacobian.temp_i)-4) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 4/6);
        xb(6+length(driver.jacobian.temp_i)-5) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 5/6);      
        xb(6+length(driver.jacobian.temp_i)-6) = dBT1231 * TfacAdjAtmosphericAmplification * iAdjLowerAtmWVfracX * (1 - 6/6);
      end

    end

  elseif settings.co2lays == 3
    xb(1) = co2x * 1;        % Set CO2 apriori lower trop
    xb(2) = co2x * 1;        % Set CO2 apriori mid trop
    xb(3) = co2x * 1;        % Set CO2 apriori strat

    xb(4) = n2ox * 1;        % set N2O 
    xb(5) = ch4x * 1;        % set CH4
    xb(6) = 0.0; %% clouds, so dunno value
    xb(7) = 0.0; %% clouds, so dunno value

    xb(8+length(driver.jacobian.water_i)) = 0.01;

  end

elseif settings.set_tracegas == +1 & driver.i16daytimestep > 0 & settings.ocb_set ~= 1
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1)/junk;
  str = ['setting time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 for ' num2str(junk) ' years']; 
  disp(str);
  if settings.co2lays == 1
    deltaT = 365/16; %% days per timestep

    %% default all this while getting good results
    xb(1) = co2x * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori
    xb(2) = n2ox * (driver.i16daytimestep-1)/deltaT * 1.0;    % set N2O 
    xb(3) = ch4x * (driver.i16daytimestep-1)/deltaT * 1.0;    % set CH4
    xb(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC11, before Aug 23 the mult was 1
    xb(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC12, before Aug 23 the mult was 1
   
  elseif settings.co2lays == 3
    deltaT = 365/16; %% days per timestep

    xb(1) = co2x * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori lower trop
    xb(2) = co2x * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori mid trop
    xb(3) = co2x * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori strat

    xb(4) = n2ox * (driver.i16daytimestep-1)/deltaT * 1.0;        % set N2O 
    xb(5) = ch4x * (driver.i16daytimestep-1)/deltaT * 1.0;        % set CH4
    xb(6) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC12, before Aug 23 the mult was 1
  end

elseif settings.set_tracegas == +2 & driver.i16daytimestep > 1 & settings.ocb_set ~= 1
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1);
  str = ['setting bootstrap time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 based on previous timestep'];
  disp(str);
  fminus = ['OutputAnomaly_OBS/' num2str(driver.iibin,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep-1) '.mat'];
  fprintf(1,'looking for %s to fill in co2/n2o/ch4/cfc11/cfc12 rates ... \n',fminus)
  if exist(fminus)
    prev = load(fminus);
    if settings.co2lays == 1
      deltaT = 365/16; %% days per timestep

      %% default all this while getting good results
      xb0(1) = co2x * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori
      xb0(2) = n2ox * (driver.i16daytimestep-1)/deltaT * 1.0;    % set N2O 
      xb0(3) = ch4x * (driver.i16daytimestep-1)/deltaT * 1.0;    % set CH4
      xb0(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC11, before Aug 23 the mult was 1
      xb0(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC12, before Aug 23 the mult was 1

      %% overwrite!!!
      xb(1) = prev.oem.finalrates(1);
      xb(2) = prev.oem.finalrates(2);
      xb(3) = prev.oem.finalrates(3);
      xb(4) = prev.oem.finalrates(4);
      xb(5) = prev.oem.finalrates(5);

      fprintf(1,'changed CO2   apriori from %8.6f to %8.6f ppv \n',xb0(1),xb(1))
      fprintf(1,'changed N2O   apriori from %8.6f to %8.6f ppb \n',xb0(2),xb(2))
      fprintf(1,'changed CH4   apriori from %8.6f to %8.6f ppb \n',xb0(3),xb(3))
      fprintf(1,'changed CRF11 apriori from %8.6f to %8.6f ppt \n',xb0(4),xb(4))
      fprintf(1,'changed CFC12 apriori from %8.6f to %8.6f ppt \n',xb0(5),xb(5))

    else
      error('not doing this')
    end
  else
    error('oops previous file DNE');
  end
end
