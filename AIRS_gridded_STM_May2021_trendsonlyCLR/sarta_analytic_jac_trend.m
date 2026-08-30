function [m_ts_jac0,nlays,qrenorm,freq2645,colo3,profilejunk] = sarta_analytic_jac_trend(driver,hMean17years,ha,pMean17years,pa,iiBin);

%% see sarta_analytic_jac.m
%% see sarta_analytic_jac_trend.m

iRunSartaJac = +1;

if iRunSartaJac < 0
  m_ts_jac0 = [];
  nlays = [];
  qrenorm = [];
  freq2645 = [];
  colo3 = [];
end
  
sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
sarta = '/home/sergio/git/sarta_scatter_rtp_klayers_sergio/JACvers/bin/jac_airs_l1c_2834_cloudy_apr26_H2024';

fip = mktempS('xxxjunk_ip');
fop = mktempS('xxxjunk_op');
frp = mktempS('xxxjunk_rp');

[havg,pavg] = subset_rtp_allcloudfields(hMean17years,pMean17years,[],[],iiBin);

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

pavg.plays = plevs2plays(pavg.plevs);

lps = compute_lapse_rate(havg,pavg);

profilejunk = pavg;
profilejunk.nlays = pavg.nlevs - 1;
profilejunk.plays = plevs2plays(pavg.plevs);
profilejunk.ptemp = pavg.ptemp;
profilejunk.gas_1 = pavg.gas_1;
profilejunk.gas_3 = pavg.gas_3;
profilejunk.stemp = pavg.stemp;
profilejunk.spres = pavg.spres;
profilejunk.lps_tropoapauseP   = lps.trp_pHI;
profilejunk.lps_tropoapauseind = lps.trp_ind;

if iRunSartaJac < 0
  disp('just wanted tropoopause info, exiting sarta_analytic_jac.m')
  return
end

rtpwrite(fop,havg,[],pavg,[]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%

%% output
m_ts_jac0 = [jac2; jac4; jac6; 0*jac2; 0*jac2; jacST; jacWV; jacTZ; jacOZ]';
nlays = iaNumLay;

[mm,nn] = size(m_ts_jac0);
qrenorm = ones(1,mm);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% see [xm_ts_jac0,xnlays,xqrenorm,xfreq2645,~,xprofilejunk]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,iVersJac,iOldORNew,topts)
qrenorm(1) = 2.2;
qrenorm(2) = 1.0;
qrenorm(3) = 5.0;
qrenorm(4) = 1.0;
qrenorm(5) = 1.0;
qrenorm(6) = 0.1;

qrenorm(7:end) = 0.01;
scale = 0.01;

m_ts_jac0 = [jac2/100/qrenorm(1); jac4/100/qrenorm(2); jac6/100/qrenorm(3); 0*jac2; 0*jac2; jacST/100/qrenorm(6);];
m_ts_jac0 = [m_ts_jac0;  jacWV*scale; jacTZ*scale; jacOZ*scale;]';

%%%%%%%%%%%%%%%%%%%%%%%%%

freq2645 = w;
colo3 = nansum(jac3,1);

