clear airsL3 era5 cmip6
clear *spectral_olr

addpath ../FIND_NWP_MODEL_TRENDS
driver_get_the_model_trends

%figure(8); plot_profile_trends          %% ERA_ERA5_AIRSL3,   with mask
%figure(8); plot_profile_trends2         %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as std(profile)
figure(8); plot_profile_trends3_cmip6    %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as mean(unc_profile),,,, I could call plot_profile_trends3_amip6 but I have cleared the variable

do_feedbacks
