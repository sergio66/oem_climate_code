%% trends

set_iMidPoint_TimeStepUse

if settings.descORasc == -1
  AHAJAC = 'NOTHING .. should this be same as settings.descORasc';
  driver.jacobian.filename = [AHAJAC];  
  fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  error('??????')
  
elseif settings.descORasc == +1
  %% we can "fool" the code by using midpoint anomaly jac
  %junk = num2str(iMidPoint_TimeStepUse,'%03d');
  %driver.jacobian.filename = ['AHA']; 

  %% for now assume same jacs
  AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/';
  AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/';
  AHA = '/home/sergio/nogit/TILEJACS/AVG/KCARTA/T_WV_O3/LatBin65/SubsetJacLatbin/';
  
  %% figure out which latbin
  %% latbin 1 has jacs/points for indices 1:72 + (1-1)*72
  %% latbin 2 has jacs/points for indices 1:72 + (2-1)*72
  driver.jac_latbin         = floor((driver.iibin-1)/72)+1;
  driver.jac_indexINSIDEbin = driver.iibin - (driver.jac_latbin-1)*72;     %% so this should be lonbin

  %% Default values
  iKCARTAorSARTA = +1;
  iVersJac = 2019;     %% ERA5 from 2002-2019
  iVersJac = 2021;     %% ERA5 from 2002-2021
  iVersJac = 2022;     %% ERA5 from 2002-2022
  iVersJac = 2025;     %% ERA5 from 2002-2025

  if settings.dataset == 30
    iVersJac = 2022;   %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)  AMSU AMSU AMSU
  elseif settings.dataset == 19 | settings.dataset == 20
    iVersJac = 2025;   %% ERA5 cldQ from 2002-2025, s
  elseif settings.dataset == 5
    iVersJac = 2014;   %% AMIP6/CMIp6 2002-2014, 12 years
  elseif settings.dataset == 6
    iVersJac = 2012;   %% CrIS NSR 2012-2019, 07 years
  elseif settings.dataset == 7
    iVersJac = 2022;   %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)
    iVersJac = 2021;   %% ERA5 CLR from 2002-2021
  elseif settings.dataset == 8
    iVersJac = 2015;   %% OCO2 2015-2021, 7 years
  elseif settings.dataset == 9

    %% Feb 9, 2023 commit
    if settings.ocb_set == 1
      disp(' settings.dataset == 9 but settings.ocb_set == 1 so set iVersJac = 2021')
      iVersJac = 2021;   %% ERA5 clr 2021
      iVersJac = 2025;   %% reset in Aug 2026 to do this
      iVersJac = 2022;   %% reset in Aug 2026 to do this                  
    elseif settings.ocb_set == 0
      iVersJac = 2022;   %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)
      iVersJac = 2021;   %% ERA5 CLR  from 2002-2021
      iVersJac = 2025;   %% reset in Aug 2026 to do this
      iVersJac = 2022;   %% reset in Aug 2026 to do this            
    end

    if settings.ocb_set == 1
      iVersJac = 2022; iOldORNew = +5;  %% ERA5 clr from 2002-2022
      iVersJac = 2025;   %% reset in Aug 2026 to do this
      iVersJac = 2022;   %% reset in Aug 2026 to do this                  
    elseif settings.ocb_set == 0
      if driver0.iQuantile < 4
        iVersJac = 2022;  iOldORNew = +9;  %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)
        iVersJac = 2022;  iOldORNew = +5;  %% ERA5 clr from 2002-2022   THIS IS NEW DEC 29, 2023 !!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!!!!
        iVersJac = 2025;   %% reset in Aug 2026 to do this
        iVersJac = 2022;   %% reset in Aug 2026 to do this                  
      elseif driver0.iQuantile >= 4
        iVersJac = 2021;                   %% ERA5 CLR  from 2002-2021
        iVersJac = 2022;  iOldORNew = +5;  %% ERA5 clr from 2002-2022
        iVersJac = 2025;   %% reset in Aug 2026 to do this
        iVersJac = 2022;   %% reset in Aug 2026 to do this                  	
      end
    end
  end
  disp(' ')
  fprintf(1,' set_driver_jacfile.m : settings.ocb_set = %2i     iVersJac = %4i    iOldORNew = %2i \n',settings.ocb_set,iVersJac,iOldORNew);
  disp(' ')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if settings.dataset == 30
    disp('set_trends_jacfile.m : reading in SARTA AMSU jacobians')  
    %% SARTA AMSU jacs
    AHA = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/AMSU_12channels_20years_Trends_Anomalies/JAC_20year_avg/'];
    AHA = [AHA '/amsu_jacT_ST_WV_' num2str(driver.iibin,'%04i') '.mat'];
    fprintf(1,'AMSU JAC %s \n',AHA)

  elseif iKCARTAorSARTA < 0
    disp('set_trends_jacfile.m : reading in SARTA jacobians')  
    %% AHA = [AHA '/subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
    AHA = [AHA '/clr_subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];

  else
    if topts.iXJac == 2
      disp('set_trends_jacfile.m : reading in trends jacs : kCARTA : DEFAULT')
    elseif topts.iXJac == 1
      disp('set_trends_jacfile.m : on the fly SARTA trends jacs')
    end

    % AHA = [AHA '/kcarta_subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];                                   %% 40 latbins
    if iVersJac == 2012 | iVersJac == 2015
      AHA = [AHA '/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_07yr_' num2str(driver.jac_latbin,'%02i') '.mat'];   %% ERA5, 2012-2019 year
    elseif iVersJac == 2014
      %% see ~/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Sept2022_startSept2002_endAug2014_trendsonly/clust_put_together_jacs_clrERA5.m
      AHA = [AHA '/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_12yr_' num2str(driver.jac_latbin,'%02i') '.mat']; %% ERA5,  2002-2014 12 year
    elseif iVersJac == 2019
      AHA = [AHA '/kcarta_clr_subjacLatBin_newSARTA_' num2str(driver.jac_latbin,'%02i') '.mat'];                     %% ERA-I, 2002-2019 17 year
    elseif iVersJac == 2021 
      AHA = [AHA '/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(driver.jac_latbin,'%02i') '.mat']; %% ERA5, 2002-2021 19 year
    elseif iVersJac == 2022
      if iOldORNew == 9
        %% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Nov2022_startSept2002_endAug2022_trendsonly_cldy_Q09/clust_put_together_jacs_cldERA5.m, but this has TONS of clouds
        AHA = [AHA '/kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(driver.jac_latbin,'%02i') '.mat']; %% ERA5,  2002-2022 20 year <avg cld = Q09> and NOT Q05
      elseif iOldORNew == 5
        %% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Mar2023_startSept2002_endAug2022_trendsonly/clust_put_together_jacs_clrERA5.m
        AHA = [AHA '/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_' num2str(driver.jac_latbin,'%02i') '.mat']; %% ERA5,  2002-2022 20 year <CLR>
      end
    elseif iVersJac == 2025
      AHA = [AHA '/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(driver.jac_latbin,'%02i') '.mat']; %% ERA5, 2002-2021 19 year	
    else
      iVersJac
      error('iVersJac = [2012,2015 = 2012/05-2019/04]  or 2014, 2019, 2021, 2022 and fake 2023 [2002/09-20XY/08] and 2025 only')
    end
  end

  topts.jacobian.filename = AHA;
  driver.jacobian.filename = AHA;
  topts.iVersJac = iVersJac;

  clear AHA
  if topts.iXJac == 2  
    fprintf(1,'reading in jac version %4i constant kcarta jac file %s \n',iVersJac,driver.jacobian.filename)
  elseif topts.iXJac == 1
    fprintf(1,'not reading in jac version %4i constant kcarta jac file %s \n',iVersJac,driver.jacobian.filename)
    disp('since we are running SARTA analytic jacs on the fly')
  end
end

