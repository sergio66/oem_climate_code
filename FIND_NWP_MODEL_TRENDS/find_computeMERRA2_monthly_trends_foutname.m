if iNumYears <= 069
  if iDorA > 0
    fout_trendjunk = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
  else
    fout_trendjunk = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_asc.mat'];
  end
else
  error('gugugug')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if iNumYears == 12
  if iDorA > 0
    %save MERRA2_atm_data_2002_09_to_2021_07_trends_desc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2014_08_trends_desc.mat comment trend*
  else
    %save MERRA2_atm_data_2002_09_to_2021_07_trends_asc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2014_08_trends_asc.mat comment trend*
  end
elseif iNumYears == 18
  if iDorA > 0
    %save MERRA2_atm_data_2002_09_to_2019_08_trends_desc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2020_08_trends_desc.mat comment trend*
  else
    %save MERRA2_atm_data_2002_09_to_2019_08_trends_asc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2020_08_trends_asc.mat comment trend*
  end
elseif iNumYears == 19
  if iDorA > 0
    %save MERRA2_atm_data_2002_09_to_2021_07_trends_desc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2021_08_trends_desc.mat comment trend*
  else
    %save MERRA2_atm_data_2002_09_to_2021_07_trends_asc.mat comment trend*
    save MERRA2_atm_data_2002_09_to_2021_08_trends_asc.mat comment trend*
  end
elseif iNumYears == 20
  if iDorA > 0
    save MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat comment trend*
  else
    save MERRA2_atm_data_2002_09_to_2022_08_trends_asc.mat comment trend*
  end
end
%}
