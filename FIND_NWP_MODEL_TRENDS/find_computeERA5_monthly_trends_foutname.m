%{
if iNumYears <= 069
  if iCldORClr == -1
    if iDorA == +1
      fout_trendjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
    elseif iDorA == -1
      fout_trendjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc.mat'];
  elseif iDorA == +10
      fout_trendjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc_randompt.mat'];
    elseif iDorA == -10
      fout_trendjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc_randompt.mat'];
    end
  else
    if iDorA == +1
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
    elseif iDorA == -1
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc.mat'];
    elseif iDorA == +10
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc_randompt.mat'];
    elseif iDorA == -10
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc_randompt.mat'];
    end
  end
else
  error('gugugug')
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iCldORClr == -1
  if iDorA == +1
    fout_trendjunk = ['ERA5_atm_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_desc.mat'];
  elseif iDorA == -1
    fout_trendjunk = ['ERA5_atm_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_asc.mat'];
  elseif iDorA == +10
    fout_trendjunk = ['ERA5_atm_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_desc_randompt.mat'];
  elseif iDorA == -10
    fout_trendjunk = ['ERA5_atm_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_asc_randompt.mat'];
  end
else
  if iDorA == +1
    fout_trendjunk = ['ERA5_atm_N_cld_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_desc.mat'];
  elseif iDorA == -1
    fout_trendjunk = ['ERA5_atm_N_cld_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_asc.mat'];
  elseif iDorA == +10
    fout_trendjunk = ['ERA5_atm_N_cld_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_desc_randompt.mat'];
  elseif iDorA == -10
    fout_trendjunk = ['ERA5_atm_N_cld_data_' num2str(yymmS(1)) '_' num2str(yymmS(2),'%02i') '_to_' num2str(yymmE(1)) '_' num2str(yymmE(2),'%02i') '_trends_asc_randompt.mat'];
  end
end  

fprintf(1,'will be saving to %s \n',fout_trendjunk)

if iAllorSeasonal == -1
  fout_trendjunk = fout_trendjunk(1:end-4);
  fout_trendjunk = [fout_trendjunk '_DJF.mat'];
elseif iAllorSeasonal == -2
  fout_trendjunk = fout_trendjunk(1:end-4);
  fout_trendjunk = [fout_trendjunk '_MAM.mat'];
elseif iAllorSeasonal == -3
  fout_trendjunk = fout_trendjunk(1:end-4);
  fout_trendjunk = [fout_trendjunk '_JJA.mat'];
elseif iAllorSeasonal == -4
  fout_trendjunk = fout_trendjunk(1:end-4);
  fout_trendjunk = [fout_trendjunk '_SON.mat'];
end

if iTrendsOrAnoms < 0
  fout_trendjunk = fout_trendjunk(1:end-4);
  fout_trendjunk = [fout_trendjunk '_anomaly.mat'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
if iCldORClr == -1
  %% CLR ONLY
  if iNumYears == 12
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2014_08_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2014_08_trends_asc.mat';
    end
  elseif iNumYears == 07
    %%%%%%%%%%%% SPECIAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_data_2012_05_to_2019_04_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_data_2012_05_to_2014_04_trends_asc.mat';
    end
    %%%%%%%%%%%% SPECIAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iNumYears == 18
    if iDorA > 0
      %fout_trendjunk = 'ERA5_atm_data_2002_09_to_2019_08_trends_desc.mat';
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2020_08_trends_desc.mat';
    else
      %fout_trendjunk = 'ERA5_atm_data_2002_09_to_2019_08_trends_asc.mat';
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2020_08_trends_asc.mat';
    end
  elseif iNumYears == 19
    if iDorA > 0
      %fout_trendjunk = 'ERA5_atm_data_2002_09_to_2021_07_trends_desc.mat';
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2021_08_trends_desc.mat';
    else
      %fout_trendjunk = 'ERA5_atm_data_2002_09_to_2021_07_trends_asc.mat';
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2021_08_trends_asc.mat';
    end
  elseif iNumYears == 20
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2022_08_trends_asc.mat';
    end
  elseif iNumYears <= 20
  else
    iNumYears
    error('unknown iNumYears .. accepting 07, 12,18,19,20')
  end
elseif iCldORClr == +1
  %% CLD fields
  if iNumYears == 12
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2014_08_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2014_08_trends_asc.mat';
    end
  elseif iNumYears == 07
    %%%%%%%%%%%% SPECIAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_N_cld_data_2012_05_to_2019_04_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_N_cld_data_2012_05_to_2014_04_trends_asc.mat';
    end
    %%%%%%%%%%%% SPECIAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iNumYears == 18
    if iDorA > 0
      %fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2019_08_trends_desc.mat';
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2020_08_trends_desc.mat';
    else
      %fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2019_08_trends_asc.mat';
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2020_08_trends_asc.mat';
    end
  elseif iNumYears == 19
    if iDorA > 0
      %fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_07_trends_desc.mat';
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_08_trends_desc.mat';
    else
      %fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_07_trends_asc.mat';
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_08_trends_asc.mat';
    end
  elseif iNumYears == 20
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat';
    end
  elseif iNumYears <= 20
    if iDorA > 0
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
    else
      fout_trendjunk = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc.mat'];
    end
  else
    iNumYears
    error('unknown iNumYears .. accepting 07, 12,18,19,20')
  end
end  
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
