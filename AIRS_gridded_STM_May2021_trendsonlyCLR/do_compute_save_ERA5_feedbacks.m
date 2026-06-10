disp(' ')
disp('WARNING : computing ERA5 feedbacks needs an UP-TO-DATE nwp_spectral_trends_era5 so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing ERA5 feedbacks needs an UP-TO-DATE nwp_spectral_trends_era5 so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing ERA5 feedbacks needs an UP-TO-DATE nwp_spectral_trends_era5 so make sure you have re-run make_profile_spectral_trends');
disp(' ')

% search for    ABC this may help in do_compute_save_ERA5_feedbacks.m    inside driver_gather_gridded_retrieval_results.m
% or     era5 = get_ERA5_trends_thermodynamic_and_spectral(min(iNumYears,23),iNorD); 
% junk = input('(+1/default) Go ahead with the calcs, you have run make_profile_spectral_trends or (-1) oops, re-run it here : ');
% if length(junk) == 0
%   junk == 1;
% end
% if junk < 0
%   clear nwp_spectral_trends_cmip6_era5_airsL3_umbc
%   nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends(cmip6,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
% end

% compute_feedbacks_airsL3_ecRad  ; pause(0.1)
% compute_feedbacks_era5_ecRad    ; pause(0.1)
% compute_feedbacks_cmip6_ecRad   ; pause(0.1)

%% results is the 4608 x 6 set of UMBC retrievals
era5_spectral_olr       = struct;    %% so it has no fields
era5_spectral_olr       = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,-1,rlat65,rlon73,-1,'ERA5');
era5_spectral_delta_olr = struct;
era5_spectral_delta_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_delta_olr,8888,rlat65,rlon73,-1,'ERA5');

if ~exist('feedbacknameNWP_ERA5')
  feedbacknameNWP_ERA5 = 'era5_only_OLR_feedbacks_23years.mat';
end  
junk = input('save ERA3 only (-1 [default]/+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk > 0
  junk2 = +1;
  if exist(feedbacknameNWP_ERA5)
    lser = ['!ls -lth ' feedbacknameNWP_ERA5];
    eval(lser);
    junk2 = input('file already exists, overwrite (-1 default/+1) : ');
    if length(junk2) == 0
      junk2 = -1;
    end
  end       
  stemptrend.era5   = era5.trend_stemp;
  saver = ['save ' feedbacknameNWP_ERA5 ' era5_spectral_olr era5_spectral_delta_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
  if junk2 > 0
    fprintf(1,'saving to %s \n',feedbacknameNWP_ERA5);
    eval(saver);
  else
    fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_ERA5);
  end
end  
