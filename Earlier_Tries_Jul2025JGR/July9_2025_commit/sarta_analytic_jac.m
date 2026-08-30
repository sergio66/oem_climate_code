function [m_ts_jac0,nlays,qrenorm,freq2645,colo3,profilejunk] = sarta_analytic_jac(driver,info_about_time_lat,iRunSartaJac);

if nargin == 2
  iRunSartaJac = +1;
end

if iRunSartaJac < 0
  m_ts_jac0 = [];
  nlays = [];
  qrenorm = [];
  freq2645 = [];
  colo3 = [];
end
  
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/TIME
sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';

fip = mktempS('xxxjunk_ip');
fop = mktempS('xxxjunk_op');
frp = mktempS('xxxjunk_rp');

yyS = 2002; yyE = 2030;
iCnt = 0;
for yy = yyS : yyE
  mmS = 01; mmE = 12;
  if yy == yyS
    mmS = 09;
  elseif yy == yyE
    mmS = 08;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;
    ddsave(iCnt) = 15;  
  end
end
yearsSince2002 = yysave + (mmsave-1)/12 + (ddsave-1)/12/30;

rtimesave = utc2taiSergio(yysave,mmsave,ddsave,12*ones(size(mmsave)));

closest = abs(rtimesave - info_about_time_lat.rtime(driver.anomalyinfo.timestep16day));
closest = find(closest == min(closest),1);

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m
era5name = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(closest,'%03d') '.mat'];

era5 = load(era5name);
hnew = era5.hnew_op;
pnew = era5.pnew_op;
if driver.anomalyinfo.global == 1
  pavg = find_average_rtp(hnew,pnew,2,1:4608);
elseif driver.anomalyinfo.global == +2
  boo = find(pnew.rlat >= -30 & pnew.rlat <= +30);
  pavg = find_average_rtp(hnew,pnew,2,boo);
else
  do_XX_YY_from_X_Y
  moo = info_about_time_lat.usethese{driver.anomalyinfo.latbin};
  boo = find(pnew.rlat >= rlat65(min(moo)) & pnew.rlat <= rlat65(max(moo)+1));
  pavg = find_average_rtp(hnew,pnew,1,boo);
end

%% there seems to be a droop after June 2024
%% driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m:223:co2ppm = 370 + 2.2*((yy+mm/12)-2002);
ppmv2expect = 370 + 2.2 * (yearsSince2002(closest)-2002);
ppmv2 = layers2ppmv(hnew,pavg,1:length(pavg.stemp),2);
i500 = find(pavg.plevs >= 500,1);
ppmv2 = ppmv2(i500);
pavg.gas_2 = pavg.gas_2 * ppmv2expect/ppmv2;

iDoClr = +1;
if iDoClr > 0
  disp('  sarta_analytic_jac.m : setting clouds = 0')
  pavg.cfrac = 0;
  pavg.cngwat = 0;
  pavg.ctype = -9999;

  pavg.cfrac2 = 0;
  pavg.cngwat2 = 0;
  pavg.ctype2 = -9999;

  pavg.cfrac12 = 0;
end

profilejunk.pavg = pavg;
profilejunk.lps = compute_lapse_rate(hnew,pavg);

if iRunSartaJac < 0
  disp('just wanted tropoopause info, exiting sarta_analytic_jac.m')
  return
end

rtpwrite(fop,hnew,[],pavg,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' listj=100,1,2,3,4,5,6'];
eval(sartaer);

[w,jacT,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacTZ'],100);
[w,jac1,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG1'],1);
[w,jac2,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG2'],2);
[w,jac3,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG3'],3);
[w,jac4,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG4'],4);
[w,jac5,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG5'],5);
[w,jac6,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG6'],6);
rmer = ['!/bin/rm ' fip ' ' fop  ' ' frp ' ' frp '_jac*' ]; eval(rmer);

jac2 = nansum(jac2,1);
jac4 = nansum(jac4,1);
jac6 = nansum(jac6,1);
jacST = jacT(iaNumLay+1,:);
jacTZ = jacT(1:iaNumLay,:);
jacWV = jac1(1:iaNumLay,:);
jacOZ = jac3(1:iaNumLay,:);

%% output
m_ts_jac0 = [jac2; jac4; jac6; 0*jac2; 0*jac2; jacST; jacWV; jacTZ; jacOZ]';
nlays = iaNumLay;
[mm,nn] = size(m_ts_jac0); qrenorm = ones(1,mm);
freq2645 = w;
colo3 = nansum(jac3,1);
