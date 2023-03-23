
if iCldORClr == -1
  %% CLR ONLY
  if iNumYears == 12
    if iDorA > 0
      foutjunk = 'ERA5_atm_data_2002_09_to_2014_08_desc.mat';
    else
      foutjunk = 'ERA5_atm_data_2002_09_to_2014_08_asc.mat';
    end
  elseif iNumYears == 18
    if iDorA > 0
      %foutjunk = 'ERA5_atm_data_2002_09_to_2019_08_desc.mat';
      foutjunk = 'ERA5_atm_data_2002_09_to_2020_08_desc.mat';
    else
      %foutjunk = 'ERA5_atm_data_2002_09_to_2019_08_asc.mat';
      foutjunk = 'ERA5_atm_data_2002_09_to_2020_08_asc.mat';
    end
  elseif iNumYears == 19
    if iDorA > 0
  %    foutjunk = 'ERA5_atm_data_2002_09_to_2021_07_desc.mat';
      foutjunk = 'ERA5_atm_data_2002_09_to_2021_08_desc.mat';
    else
  %    foutjunk = 'ERA5_atm_data_2002_09_to_2021_07_asc.mat';
      foutjunk = 'ERA5_atm_data_2002_09_to_2021_08_asc.mat';
    end
  elseif iNumYears == 20
    if iDorA > 0
      foutjunk = 'ERA5_atm_data_2002_09_to_2022_08_desc.mat';
    else
      foutjunk = 'ERA5_atm_data_2002_09_to_2022_08_asc.mat';
    end
  else
    iNumYears
    error('unknown iNumYears .. accepting 12,18,19')
  end
elseif iCldORClr == +1
  %% CLD fields as well
  if iNumYears == 12
    if iDorA > 0
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2014_08_desc.mat';
    else
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2014_08_asc.mat';
    end
  elseif iNumYears == 18
    if iDorA > 0
      %foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2019_08_desc.mat';
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2020_08_desc.mat';
    else
      %foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2019_08_asc.mat';
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2020_08_asc.mat';
    end
  elseif iNumYears == 19
    if iDorA > 0
  %    foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_07_desc.mat';
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_08_desc.mat';
    else
  %    foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_07_asc.mat';
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2021_08_asc.mat';
    end
  elseif iNumYears == 20
    if iDorA > 0
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_desc.mat';
    else
      foutjunk = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_asc.mat';
    end
  else
    iNumYears
    error('unknown iNumYears .. accepting 12,18,19')
  end
end
