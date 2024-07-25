if iForwardFrom2002_or_BackwardFrom2022  > 0
  if iNumYears <= 069
    if iCldORClr == -1
      if iDorA == +1
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_night.mat'];
      elseif iDorA == -1
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_day.mat'];
      elseif iDorA == 0
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_night.mat'];
      elseif iDorA == +10
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_night_randompt.mat'];
      elseif iDorA == -10
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_randompt.mat'];
      elseif iDorA == 100
        foutjunk = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_night_randompt.mat'];
      end
    elseif iCldORClr == +1
      if iDorA == +1
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_night.mat'];
      elseif iDorA == -1
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_day.mat'];
      elseif iDorA == 0
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_night.mat'];
      elseif iDorA == +10
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_night_randompt.mat'];
      elseif iDorA == -10
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_randompt.mat'];
      elseif iDorA == 100
        foutjunk = ['ERA5_atm_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_day_night_randompt.mat'];
      end
    end
  end
else
  error('gugugug')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
