addpath /home/sergio/MATLABCODE/TIME

dir0 = 'gml.noaa.gov/aftp/data/ozwv/WaterVapor/Hilo_New/';
thedir = dir([dir0 '/HIH*.txt']);

fprintf(1,'gong to read %3i files \n',length(thedir));

iCnt = 0;
for ii = 1 : length(thedir)
  if mod(ii,100) == 0
    fprintf(1,'+');
  elseif mod(ii,10) == 0
    fprintf(1,'.');
  end
  fname = thedir(ii).name;
  junkdate = str2num(fname(9:16));
  Ax = read_hilo(dir0,fname);
  if length(Ax) > 0 & junkdate > 20020000
    iaFound(ii) = +1;
    iCnt = iCnt + 1;
    A(iCnt,:,:) = Ax;
    thedate(iCnt,:) = str2num(fname(9:16));
  else
    iaFound(ii) = 0;
  end
end
fprintf(1,'\n');

addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
junk = num2str(thedate);
yy = str2num(junk(:,1:4));
mm = str2num(junk(:,5:6));
dd = str2num(junk(:,7:8));
doy = change2days(yy,mm,dd,2002);
hilo_pav = nanmean(squeeze(A(:,:,3)),1);

oo = find(doy >= 8*30); %% Aug 30,2002

for ii = 1 : 20
  [B stats] = Math_tsfit_lin_robust(double(doy),double(squeeze(A(oo,ii,4))),4);
  hilo_T_trend_trend(ii) = B(2);
  hilo_T_trend_err(ii) = stats.se(2);

  [B stats] = Math_tsfit_lin_robust(double(doy),double(squeeze(A(oo,ii,6))),4);
  hilo_RH_trend_trend(ii) = B(2);
  hilo_RH_trend_err(ii) = stats.se(2);

  junk = squeeze(A(oo,ii,8));
  [B stats] = Math_tsfit_lin_robust(double(doy),double(junk),4);
  hilo_ppmv_trend_trend(ii) = B(2);
  hilo_ppmv_trend_err(ii) = stats.se(2);

  hilo_ppmv_avg(ii) = nanmean(junk);
  junk = junk/nanmean(junk);
  [B stats] = Math_tsfit_lin_robust(double(doy),double(junk),4);
  hilo_wvfrac_trend_trend(ii) = B(2);
  hilo_wvfrac_trend_err(ii) = stats.se(2);
end

clf; plot(hilo_wvfrac_trend_trend.*hilo_ppmv_avg,hilo_pav,'o-',hilo_ppmv_trend_trend,hilo_pav); set(gca,'ydir','reverse'); plotaxis2;
subplot(131); plot(hilo_T_trend_trend,hilo_pav);      set(gca,'ydir','reverse'); plotaxis2; title('dT/dt K/yr')
subplot(132); plot(hilo_RH_trend_trend,hilo_pav);     set(gca,'ydir','reverse'); plotaxis2; title('dRH/dt %/yr')
subplot(133); plot(hilo_wvfrac_trend_trend,hilo_pav); set(gca,'ydir','reverse'); plotaxis2; title('dwvfrac/dt 1/yr')

save hilo.mat hilo_*

