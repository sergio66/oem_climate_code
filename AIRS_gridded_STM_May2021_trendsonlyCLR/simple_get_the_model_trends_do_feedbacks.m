clear airsL3 era5 cmip6
clear *spectral_olr
clear nwp_spectral_trends_cmip6_era5_airsL3_umbc

addpath ../FIND_NWP_MODEL_TRENDS
driver_get_the_model_trends

%figure(8); plot_profile_trends          %% ERA_ERA5_AIRSL3,   with mask
%figure(8); plot_profile_trends2         %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as std(profile)
figure(8); plot_profile_trends3_cmip6    %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as mean(unc_profile),,,, I could call plot_profile_trends3_amip6 but I have cleared the variable

junk = input('Do you want to compute OLR feedbacks (-1 default /+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk > 0
  if ~exist('lps0')
    addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
    [mmw0,lps0] = mmwater_rtp(h,p);
  end
  do_feedbacks
end
