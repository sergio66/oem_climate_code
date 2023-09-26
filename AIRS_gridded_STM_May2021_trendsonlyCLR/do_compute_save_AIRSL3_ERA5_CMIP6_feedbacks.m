disp(' ')
disp('WARNING : computing AIRSL3_ERA5_CMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing AIRSL3_ERA5_CMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing AIRSL3_ERA5_CMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp(' ')

junk = input('(+1/default) Go ahead with the calcs, you have run make_profile_spectral_trends or (-1) oops, re-run it here : ');
if length(junk) == 0
  junk == 1;
end
if junk < 0
  clear nwp_spectral_trends_cmip6_era5_airsL3_umbc
  nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends(cmip6,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
end

% compute_feedbacks_airsL3_ecRad  ; pause(0.1)
% compute_feedbacks_era5_ecRad    ; pause(0.1)
% compute_feedbacks_cmip6_ecRad   ; pause(0.1)

era5_spectral_olr = struct;    %% so it has no fields
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,-1,'ERA5');

aL3trend.stemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
aL3trend.ptemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
aL3trend.gas_1 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1;
aL3trend.gas_3 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3;
airsL3_spectral_olr = struct;    %% so it has no fields
airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,-1,'AIRS L3');

c6trend.stemp = cmip6.trend_stemp;
c6trend.ptemp = cmip6.trend_ptemp;
c6trend.gas_1 = cmip6.trend_gas_1;
c6trend.gas_3 = cmip6.trend_gas_3;
cmip6_spectral_olr = struct;    %% so it has no fields
cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,-1,'CMIP6');

junk = input('save ERA3/AIRSL3/CMIP6 only (-1 [default]/+1) : ');
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
  stemptrend.airsL3 = aL3trend.stemp;
  stemptrend.cmip6  = c6trend.stemp;
  saver = ['save ' feedbacknameNWP_ERA5 ' era5_spectral_olr cmip6_spectral_olr airsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
  if junk2 > 0
    fprintf(1,'saving to %s \n',feedbacknameNWP_ERA5);
    eval(saver);
  else
    fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_ERA5);
  end
end  
