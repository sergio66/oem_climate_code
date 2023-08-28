strMODELScopy = strMODELS;

if exist('era5')
  era5copy = era5;
end
if exist('cmip6')
  cmip6copy = cmip6;
end
if exist('airsL3')
  airsL3copy = airsL3;
end

disp(' ')
driver_get_the_model_trends

strMODELSX = strMODELS;  
strMODELS  = strMODELScopy;

if exist('era5copy')
  merra2 = era5;
  era5   = era5copy;
end
if exist('cmip6copy')
  amip6 = cmip6;
  cmip6 = cmip6copy;
end
if exist('airsL3copy')
  if isfield(airsL3,'airsL3_spectral_rates')
    airsL3.climcapsL3_spectral_rates = airsL3.airsL3_spectral_rates;
    airsL3 = rmfield(airsL3,'airsL3_spectral_rates');
  end
  climcapsL3 = airsL3;
  climcapsl3 = airsL3;
  airsL3     = airsL3copy;
end

