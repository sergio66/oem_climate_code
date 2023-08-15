if iNumYears <= 069
  if iDorA > 0
    foutjunk = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_desc.mat'];
  else
    foutjunk = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_asc.mat'];
  end
else
  error('gugugug')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
elseif iNumYears == 18
  if iDorA > 0
    %save -v7.3 MERRA2_atm_data_2002_09_to_2019_08_desc.mat comment all
    save -v7.3 MERRA2_atm_data_2002_09_to_2020_08_desc.mat comment all
  else
    %save -v7.3 MERRA2_atm_data_2002_09_to_2019_08_asc.mat comment all
    save -v7.3 MERRA2_atm_data_2002_09_to_2020_08_asc.mat comment all
  end
elseif iNumYears == 19
  if iDorA > 0
%    save -v7.3 MERRA2_atm_data_2002_09_to_2021_07_desc.mat comment all
    save -v7.3 MERRA2_atm_data_2002_09_to_2021_08_desc.mat comment all
  else
%    save -v7.3 MERRA2_atm_data_2002_09_to_2021_07_asc.mat comment all
    save -v7.3 MERRA2_atm_data_2002_09_to_2021_08_asc.mat comment all
  end
elseif iNumYears == 20
  if iDorA > 0
    save -v7.3 MERRA2_atm_data_2002_09_to_2022_08_desc.mat comment all
  else
    save -v7.3 MERRA2_atm_data_2002_09_to_2022_08_asc.mat comment all
  end
else
  iNumYears
  error('unknown iNumYears .. accepting 12,18,19,20')
end
%}
