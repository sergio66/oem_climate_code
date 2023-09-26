disp(' ')
disp('WARNING : computing CLIMCAPSL3_MERRA2_AMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing CLIMCAPSL3_MERRA2_AMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp('WARNING : computing CLIMCAPSL3_MERRA2_AMIP6 feedbacks needs an UP-TO-DATE nwp_spectral_trends_cmip6_era5_airsL3_umbc so make sure you have re-run make_profile_spectral_trends');
disp(' ')

junk = input('(+1/default) Go ahead with the calcs, you have run make_profile_spectral_trends or (-1) oops, re-run it here : ');
if length(junk) == 0
  junk == 1;
end
if junk < 0
  clear nwp_spectral_trends_amip6_merra2_climcapsL3_umbc
  nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = make_profile_spectral_trends(amip6,merra2,climcapsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
end

% compute_feedbacks_airsL3_ecRad  ; pause(0.1)
% compute_feedbacks_era5_ecRad    ; pause(0.1)
% compute_feedbacks_cmip6_ecRad   ; pause(0.1)

merra2_spectral_olr = struct;    %% so it has no fields
merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,-1,'MERRA2');

if ~exist('nwp_spectral_trends_amip6_merra2_climcapsL3_umbc')
 nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = make_profile_spectral_trends(amip6,merra2,climcapsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
end
cL3trend.stemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.stemp;
cL3trend.ptemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.ptemp;
cL3trend.gas_1 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_1;
cL3trend.gas_3 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_3;
climcapsL3_spectral_olr = struct;    %% so it has no fields
climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,cL3trend.ptemp,cL3trend.gas_1,cL3trend.gas_3,climcapsL3_spectral_olr,-1,'CLIMCAPS L3');

a6trend.stemp = amip6.trend_stemp;
a6trend.ptemp = amip6.trend_ptemp;
a6trend.gas_1 = amip6.trend_gas_1;
a6trend.gas_3 = amip6.trend_gas_3;
amip6_spectral_olr = struct;    %% so it has no fields
amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,a6trend.ptemp,a6trend.gas_1,a6trend.gas_3,amip6_spectral_olr,-1,'AMIP6');

junk = input('save MERRA2/CLIMCAPSL3/AMIP6 only (-1 [default]/+1) : ');
if length(junk) == 0
  junk = -1;
end
if junk > 0
  junk2 = +1;
  if exist(feedbacknameNWP_MERRA2)
    lser = ['!ls -lth ' feedbacknameNWP_MERRA2];
    eval(lser);
    junk2 = input('file already exists, overwrite (-1 default/+1) : ');
    if length(junk2) == 0
      junk2 = -1;
    end
  end       
  stemptrend.era5   = merra2.trend_stemp;
  stemptrend.airsL3 = cL3trend.stemp;
  stemptrend.cmip6  = a6trend.stemp;
  saver = ['save ' feedbacknameNWP_MERRA2 ' merra2_spectral_olr amip6_spectral_olr climcapsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
  if junk2 > 0
    fprintf(1,'saving to %s \n',feedbacknameNWP_MERRA2);
    eval(saver);
  else
    fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_MERRA2);
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(feedbacknameNWP_MERRA2,'AIRSL3_MERRA2_AMIP6')
  disp('warning : saved or saving a file which has "AIRSL3_MERRA2_AMIP6" but you should change it to eg "CESM3_MERRA2_AMIP6" ')
end
