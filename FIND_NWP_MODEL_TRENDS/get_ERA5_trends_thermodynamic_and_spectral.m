function era5 = get_ERA5_trends_thermodynamic_and_spectral(iNumYears);

load('../AIRS_gridded_STM_May2021_trendsonlyCLR/h2645structure.mat');

if nargin == 0
  iNumYears = 20;
end

if iNumYears >= 0
  iNumYearsX = roundN(iNumYears,5);  %% only have spectra for 5,10,15,20 ie 2022-2007,2012,2017,2022
  if iNumYearsX == 0
    iNumYearsX = 5;
  end
else
  iNumYearsX = -4; %% only have spectra for -4   ie 2018-2022
end
if iNumYears ~= iNumYearsX
  fprintf(1,' WARNING : get_ERA5_trends_thermodynamic_and_spectral.m using ERA5 trends from %2i years in stead of %2i years \n',iNumYearsX,iNumYears);
end

if iNumYearsX < 0
  iStart2002 = -1;
  titlestr = ['ERA5 const tracegas ' num2str(2022 - abs(iNumYearsX)) '/08-2022/08'];
  iEorMstrSPECTRA = 'NIGHTorAVG';
  iEorMFstr = '/ERA5_ConstG/reconstruct_era5_spectra_geo_rlat';
else
  iStart2002 = +1;
  titlestr = ['ERA5 vary tracegas 2002/08-' num2str(2002 + abs(iNumYearsX)) '/08'];
  iEorMstrSPECTRA = 'NIGHTorAVG';
  iEorMFstr = '/ERA5/reconstruct_era5_spectra_geo_rlat';
end

era5 = getdata_NWP(5,1,0,iNumYearsX);

%% from get_the_model_spectral_trends.m
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    if iStart2002 == -1
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2018_09_2022_08.mat'];
    end
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    era5.era5_spectral_rates(:,ind) = thesave;
  end

era5.vchan = h.vchan;
figure(5); clf; plot(era5.vchan,nanmean(era5.era5_spectral_rates,2)); plotaxis2; title(titlestr);
axis([640 1640 -0.1 +0.1])
