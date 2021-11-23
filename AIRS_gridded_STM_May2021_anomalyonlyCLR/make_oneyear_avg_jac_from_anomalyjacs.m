%% see strow_override_defaults_latbins_AIRS_fewlays.m, assume iXJac = 2

num_per_year = ceil(365/16);
year_offset = 2010-2002+1;

iCnt = 0;
for ii = year_offset*num_per_year : (year_offset+1)*num_per_year - 1
  ii
  iCnt = iCnt + 1;
  %junk = num2str(driver.i16daytimestep,'%03d');
  junk = num2str(ii,'%03d');
  driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];
  loader = ['x = load(''' driver.jacobian.filename ''');'];
  eval(loader);
  xall(iCnt,:,:,:) = x.M_TS_jac_all;
end

M_TS_jac_all = squeeze(nanmean(xall,1));
f = x.f;
qrenorm = x.qrenorm;
str1 = x.str1;
str2 = x.str2;
comment = 'see make_oneyear_avg_jac_from_anomalyjacs.m';
save /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/avgjac2010_kcarta__M_TS_jac_all_5_97_97_97_2645.mat f M_TS_jac_all qrenorm str1 str2 comment
