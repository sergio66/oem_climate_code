[co2x,n2ox,ch4x] = get_co2_n2o_ch4_for_strow_override(driver,iVersJac); %% sets co2x,n2ox,ch4x

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
      %% land, so whos know
      dBT1231_WV = 0.0;
    end
    
    fprintf(1,'d/dt BT1231 from rates = %8.6f K/year --> WVfractional change needed for constant RH over ocean %8.6f frac/yr with iAdjLowerAtmWVfrac = %8.6f\n',driver.rateset.rates(1520),dBT1231_WV,iAdjLowerAtmWVfrac);
    if iAdjLowerAtmWVfrac > eps
      iVers = 0;  %% before noon Mar 11,2023
      iVers = 1;  %% after  noon Mar 11,2023
      
      if iVers == 0
        xb(6+length(driver.jacobian.water_i)-0) = dBT1231_WV * iAdjLowerAtmWVfrac;
        xb(6+length(driver.jacobian.water_i)-1) = dBT1231_WV * iAdjLowerAtmWVfrac;
      elseif iVers == 1
        xb(6+length(driver.jacobian.water_i)-0) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 0/6);
        xb(6+length(driver.jacobian.water_i)-1) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 1/6);
        xb(6+length(driver.jacobian.water_i)-2) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 2/6);
        xb(6+length(driver.jacobian.water_i)-3) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 3/6);
        xb(6+length(driver.jacobian.water_i)-4) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 4/6);
        xb(6+length(driver.jacobian.water_i)-5) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 5/6);      
        xb(6+length(driver.jacobian.water_i)-6) = dBT1231_WV * iAdjLowerAtmWVfrac * (1 - 6/6);
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
