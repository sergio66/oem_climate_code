function [m_ts_jac0,nlays,qrenorm,freq2645,colo3,profilejunk] = sarta_analytic_jac(driver,info_about_time_lat,iRunSartaJac);

%% see sarta_analytic_jac.m
%% see sarta_analytic_jac_trend.m

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

iH20XY = 20;
iH20XY = 24;

if iH20XY == 20
  %% to use "same" jacs as in the JGR 2025 paper, which is H2020, CKD 3.2
  sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
  sarta = '/home/sergio/git/sarta_scatter_rtp_klayers_sergio/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020';
elseif iH20XY == 24
  %% to use H2024, CKD 4.3
  sarta = '/home/sergio/git/sarta_scatter_rtp_klayers_sergio/JACvers/bin/jac_airs_l1c_2834_cloudy_apr26_H2024';
end

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

%% see /home/sergio/git/rtpmake/CLUST_RTPMAKE/CLUSTMAKE_ERA5_MONTHLY/clust_loop_make_monthly_tile_center_asc_or_desc.m
fdir0 = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/';
fdir0 = '/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/';
fdir0 = '/home/sergio/nogit/TILES/ERA5_profiles/';

era5name = [fdir0 '/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(closest,'%03d') '.mat'];

fprintf(1,'ERA5 name in sarta_analytic_jac.m for anomalies = %s \n',era5name)

era5 = load(era5name);
hnew = era5.hnew_op;
pnew = era5.pnew_op;
if driver.anomalyinfo.global == 1
  disp('doing global avg of ERA5 file, cosweighted')
  pavg = find_average_rtp(hnew,pnew,2,1:4608);
elseif driver.anomalyinfo.global == +2
  disp('doing TRP avg of ERA5 file, cosweighted')  
  boo = find(pnew.rlat >= -30 & pnew.rlat <= +30);
  pavg = find_average_rtp(hnew,pnew,2,boo);
else
  disp('doing zonal avg of ERA5 file, cosweighted')  
  do_XX_YY_from_X_Y
  moo = info_about_time_lat.usethese{driver.anomalyinfo.latbin};
  boo = find(pnew.rlat >= rlat65(min(moo)) & pnew.rlat <= rlat65(max(moo)+1));
  pavg = find_average_rtp(hnew,pnew,1,boo);
end

%% there seems to be a drop after June 2024
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

[w,xjacT,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacTZ'],100);
[w,xjac1,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG1'],1);
[w,xjac2,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG2'],2);
[w,xjac3,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG3'],3);
[w,xjac4,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG4'],4);
[w,xjac5,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG5'],5);
[w,xjac6,iaProf,iaNumLay] = readsarta_jacV2([frp '_jacG6'],6);
rmer = ['!/bin/rm ' fip ' ' fop  ' ' frp ' ' frp '_jac*' ]; eval(rmer);

jac2 = nansum(xjac2,1);
jac4 = nansum(xjac4,1);
jac6 = nansum(xjac6,1);
jacST = xjacT(iaNumLay+1,:);
jacTZ = xjacT(1:iaNumLay,:);
jacWV = xjac1(1:iaNumLay,:);
jacOZ = xjac3(1:iaNumLay,:);

%% see read_fileMean17years.m for h,p
%iNumYears = driver.iNumYears;
%read_fileMean17years
%[hMean17years,ha,pMean17years,pa] = rtpread(fileMean17years);
%[havg,pavg] = subset_rtp_allcloudfields(hMean17years,pMean17years,[],[],iiBin);
havg = hnew;
pavg = pavg;

%xjac3 = jacOZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% output after renormalizing
qrenorm_sarta_analyticjac

