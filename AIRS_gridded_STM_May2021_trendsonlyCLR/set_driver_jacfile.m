function [driver,iVersJac,iXJac] = set_driver_jacfile(driver0,settings);

driver = driver0;

iXJac = settings.iXJac;
%if driver.i16daytimestep > 0
%  iXJac = 0; %% const geo kcarta jcs, default for trends
%  iXJac = 1; %% varying geo sarta jacs
%  iXJac = 2; %% varying geo kcarta jacs, default for anomaly
%  iXJac = 3; %% should really be Q(X --> 1) but hard to get ERA5 conditions for this!
%end

if driver.i16daytimestep < 0
  iMidPoint_TimeStepUse = 1;            %% put in constant Jacobian, at timestep 1 (2002/09)
  iMidPoint_TimeStepUse = 388;          %% put in constant Jacobian, at timestep 365 (2018/08)
  iMidPoint_TimeStepUse = floor(388/2); %% put in constant Jacobian, half way through (365/2 ==> 2009/09)

  iMidPoint_TimeStepUse = 270; 
  iMidPoint_TimeStepUse = 335; 
  iMidPoint_TimeStepUse = 194; %% default "fool" the code by using midpoint anomaly jac

  if settings.descORasc == -1
    driver.jacobian.filename = [AHAJAC];
    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)

  elseif settings.descORasc == +1
    %% we can "fool" the code by using midpoint anomaly jac
    %junk = num2str(iMidPoint_TimeStepUse,'%03d');
    %driver.jacobian.filename = ['AHA']; 

    %% for now assume same jacs
    AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/';
    AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/';
    %% figure out which latbin
    %% latbin 1 has jacs/points for indices 1:72 + (1-1)*72
    %% latbin 2 has jacs/points for indices 1:72 + (2-1)*72
    driver.jac_latbin         = floor((driver.iibin-1)/72)+1;
    driver.jac_indexINSIDEbin = driver.iibin - (driver.jac_latbin-1)*72;     %% so this should be lonbin

    iKCARTAorSARTA = +1;
    iVersJac = 2019;     %% ERA5 from 2002-2019
    iVersJac = 2021;     %% ERA5 from 2002-2021

    if settings.dataset == 5
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
      elseif settings.ocb_set == 0
        iVersJac = 2022;   %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)
        iVersJac = 2021;   %% ERA5 CLR  from 2002-2021
      end

      if settings.ocb_set == 1
        iVersJac = 2022; iOldORNew = +5;  %% ERA5 clr from 2002-2022
      elseif settings.ocb_set == 0
        iVersJac = 2021;   %% ERA5 CLR  from 2002-2021
        iVersJac = 2022;  iOldORNew = +5;  %% ERA5 clr from 2002-2022
        iVersJac = 2022;  iOldORNew = +9;  %% ERA5 cldQ from 2002-2022, so use for Q-8 etc (cloudy)
      end

    end

    if iKCARTAorSARTA < 0
      %% AHA = [AHA '/subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
      AHA = [AHA '/clr_subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
    else
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
        iOldORNew
      else
        iVersJac
        error('iVersJac = [2012,2015 = 2012/05-2019/04]  or 2014, 2019, 2021, 2022 and fake 2023 [2002/09-20XY/08] only')
      end
    end

    topts.jacobian.filename = AHA;
    driver.jacobian.filename = AHA;
    clear AHA
    fprintf(1,'reading in jac version %4i constant kcarta jac file %s \n',iVersJac,driver.jacobian.filename)
  end

elseif driver.i16daytimestep > 0
  junk = num2str(driver.i16daytimestep,'%03d');
  if iXJac == 1
    %% sarta time vary jacs
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 1 reading in timestep sarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 2  
    %% kcarta time vary jac
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 2 reading in timestep kcarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 0
    %% constant kcarta jacs
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 0 reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
end
