load STS/NIGHTorAVG/ERA5//20yrs/reconstruct_era5_spectra_geo_rlat32_2002_09_2022_08.mat
plot(fchanx,nanmean(thesave.xtrendSpectral,2))
xlim([645 1645])

iNumYears = 20;
iNumYears = 23;
iNorD = 1; %% night

for jj = 1 : 64
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  
  ind = (1:72) + (jj-1)*72;
  if iNorD == +1
    loader = ['load STS/NIGHTorAVG/ERA5//' num2str(iNumYears) 'yrs/reconstruct_era5_spectra_geo_rlat' num2str(jj,'%02d') '_2002_09_' num2str(2002 + iNumYears) '_08.mat'];
    eval(loader)
    era5_rates_desc(:,ind)     = thesave.xtrendSpectral;
    era5_rates_desc_unc(:,ind) = thesave.xtrendSpectral_unc;    
  elseif iNorD == -1  
    loader = ['load STS/NIGHTorAVG/ERA5//' num2str(iNumYears) 'yrs/reconstruct_era5_spectra_geo_rlat' num2str(jj,'%02d') '_2002_09_' num2str(2002 + iNumYears) '_08_asc.mat'];
    eval(loader)    
    era5_rates_asc(:,ind)     = thesave.xtrendSpectral;
    era5_rates_asc_unc(:,ind) = thesave.xtrendSpectral_unc;    
  end
end

if iNorD > 0
  saver = ['save STS/NIGHTorAVG/ERA5//' num2str(iNumYears) 'yrs/ERA5_spectraltrends_2002_09_' num2str(2002 + iNumYears) '_08.mat era5_rates_desc era5_rates_desc_unc fchanx'];
else
  saver = ['save STS/NIGHTorAVG/ERA5//' num2str(iNumYears) 'yrs/ERA5_spectraltrends_2002_09_' num2str(2002 + iNumYears) '_08_asc.mat era5_rates_asc era5_rates_asc_unc fchanx'];
end
eval(saver)
