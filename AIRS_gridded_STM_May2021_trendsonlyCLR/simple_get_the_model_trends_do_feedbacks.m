clear airsL3 era5 cmip6
clear *spectral_olr
clear nwp_spectral_trends_cmip6_era5_airsL3_umbc

if ~exist('iNumYears')
  if a.topts.dataset == 4
    iNumYears = 19; %% 2002 - 2021
  elseif a.topts.dataset == 5
    iNumYears = 12;  %% 2002 - 2014
  elseif a.topts.dataset == 9
    iNumYears = 20;  %% 2002 - 2022
  elseif a.topts.dataset == 10
    iNumYears = 05;  %% 2002 - 2022
  elseif a.topts.dataset == 11
    iNumYears = 10;  %% 2002 - 2022
  elseif a.topts.dataset == 12
    iNumYears = 15;  %% 2002 - 2022
  end
end

if ~exist('iNumYears')
  a.topts.dataset
  fprintf(1,'simple_get_the_model_trends_do_feedbacks.m : iNumYears DNE : so setting iNumYears = 20 for a.topts.dataset = %2i \n',a.topts.dataset)
  iNumYears = 20; %% 2002 - 2021
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../FIND_NWP_MODEL_TRENDS
disp(' ')
driver_get_the_model_trends
figure(8); plot_profile_trends3_cmip6    %% CMIP6_ERA5_AIRSL3,       with mask, but uses mean as mean(profile) and unc as mean(unc_profile),,,, I could call plot_profile_trends3_amip6 but I have cleared the variable

junkCLIMCAPS_MERRA_AMIP = input('if you choose eg AIRSL3/ERA5/CMIP6 .. can have another go to get eg CLIMCAPSL3/MERRA2/XMIP6 : 20 years only right now : (-1 [default] / +1) : ');
if length(junkCLIMCAPS_MERRA_AMIP) == 0
  junkCLIMCAPS_MERRA_AMIP = -1;
end
if junkCLIMCAPS_MERRA_AMIP > 0
  get_the_MERRA2_CLIMCAPSL3_AMIP6_model_trends
  figure(8); plot_profile_trends3_amip6    %% AMIP6_MERRA2_CLIMCAPSL3, with mask, but uses mean as mean(profile) and unc as mean(unc_profile),,,, I could call plot_profile_trends3_amip6 but I have cleared the variable
end

%figure(8); plot_profile_trends          %% ERA_ERA5_AIRSL3,         with mask
%figure(8); plot_profile_trends2         %% CMIP6_ERA5_AIRSL3,       with mask, but uses mean as mean(profile) and unc as std(profile)
  
junk = input('Do you want to compute OLR feedbacks using ecRad (-1 /+1 [all default]) : ');
if length(junk) == 0
  junk = +1;
end
if junk > 0
  if ~exist('lps0')
    addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
    [mmw0,lps0] = mmwater_rtp_pstop_lapse(h,p);
  end
  %do_feedbacks
  do_feedbacks_wrt_globalSST
end
