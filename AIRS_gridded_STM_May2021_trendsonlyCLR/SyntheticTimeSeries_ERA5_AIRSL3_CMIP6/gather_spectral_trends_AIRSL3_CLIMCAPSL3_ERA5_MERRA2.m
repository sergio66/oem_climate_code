fprintf(1,'\n\n doing AIRS DESC ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['AIRSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save AIRSL3_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat const trend trenderr

fprintf(1,'\n\n doing AIRS ASC  ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['AIRSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_asc_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save AIRSL3_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat const trend trenderr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n\n doing CLIMCAPS DESC ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['CLIMCAPSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save CLIMCAPSL3_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat const trend trenderr

fprintf(1,'\n\n doing CLIMCAPS ASC  ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['CLIMCAPSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_asc_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save CLIMCAPSL3_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat const trend trenderr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n\n doing ERA5 DESC ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save ERA5_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat const trend trenderr

fprintf(1,'\n\n doing ERA5 ASC  ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_asc_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save ERA5_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat const trend trenderr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n\n doing MERRA2 AVG  ')
const    = nan(2645,4608);
trend    = nan(2645,4608);
trenderr = nan(2645,4608);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  x = load(['MERRA2_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02i') '//sarta_spectral_trends_latbin' num2str(ii,'%02i') '_2002_09_2022_08.mat']);
  ind = (ii-1)*72 + (1:72);
  const(:,ind)    = x.thesave.const;
  trend(:,ind)    = x.thesave.xtrend;
  trenderr(:,ind) = x.thesave.xtrendErr;
end
save MERRA2_SARTA_SPECTRAL_RATES/all_4608_2002_09_2022_08.mat const trend trenderr

fprintf(1,'\n\n DONE \n')
