if iCldORClr == -1
  %% CLR ONLY
  if iNumYears == 12
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2014_08_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_data_2002_09_to_2014_08_trends_asc.mat';
    end
  elseif iNumYears == 07
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_data_2012_05_to_2019_04_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_data_2012_05_to_2014_04_trends_asc.mat';
    end
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
    if iDorA > 0
      fout_trendjunk = 'ERA5_atm_N_cld_data_2012_05_to_2019_04_trends_desc.mat';
    else
      fout_trendjunk = 'ERA5_atm_N_cld_data_2012_05_to_2014_04_trends_asc.mat';
    end
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
  else
    iNumYears
    error('unknown iNumYears .. accepting 07, 12,18,19,20')
  end
end  
