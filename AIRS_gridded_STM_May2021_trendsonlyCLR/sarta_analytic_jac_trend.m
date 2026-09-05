function [m_ts_jac0,nlays,qrenorm,freq2645,colo3,profilejunk] = sarta_analytic_jac_trend(driver,hMean17years,ha,pMean17years,pa,iiBin);

%% see sarta_analytic_jac.m        %% for anomalies
%% see sarta_analytic_jac_trend.m  %% for trends

iRunSartaJac = +1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% output after renormalizing
qrenorm_sarta_analyticjac
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFiniteDiff = +1;
iFiniteDiff = -1;
if iFiniteDiff > 0
  boo = load(['/umbc/rs/pi_sergio/WorkDirDec2025/oem_climate_code/AIRS_gridded_STM_May2021_trendsonlyCLR//AllDemJacsClrCol_20yrs/individual_prof_convolved_kcarta_airs_' num2str(iiBin) '_coljac.mat']);
  
  pxavg = pavg;
  rtpwrite(fop,havg,[],pxavg,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
  eval(sartaer);
  [hhh,~,ppp,~] = rtpread(frp);
  t0 = rad2bt(hhh.vchan,ppp.rcalc);

  pxavg = pavg;
  pxavg.gas_2 = pxavg.gas_2 * 1.1;
  rtpwrite(fop,havg,[],pxavg,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
  eval(sartaer);
  [hhh,~,ppp,~] = rtpread(frp);
  txavg = rad2bt(hhh.vchan,ppp.rcalc);  
  finitediff2 = (txavg-t0)/log(1.1);

  pxavg = pavg;
  pxavg.gas_4 = pxavg.gas_4 * 1.1;
  rtpwrite(fop,havg,[],pxavg,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
  eval(sartaer);
  [hhh,~,ppp,~] = rtpread(frp);
  txavg = rad2bt(hhh.vchan,ppp.rcalc);    
  finitediff4 = (txavg-t0)/log(1.1);

  pxavg = pavg;
  pxavg.gas_6 = pxavg.gas_6 * 1.1;
  rtpwrite(fop,havg,[],pxavg,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
  eval(sartaer);
  [hhh,~,ppp,~] = rtpread(frp);
  txavg = rad2bt(hhh.vchan,ppp.rcalc);    
  finitediff6 = (txavg-t0)/log(1.1);

  pxavg = pavg;
  pxavg.stemp = pxavg.stemp + 0.1;
  rtpwrite(fop,havg,[],pxavg,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
  eval(sartaer);
  [hhh,~,ppp,~] = rtpread(frp);
  txavg = rad2bt(hhh.vchan,ppp.rcalc);    
  finitediffST = (txavg-t0)/(0.1);

  figure(1); clf; plot(hhh.vchan,jac2,hhh.vchan,finitediff2,boo.fKc,boo.rKc(:,1));   xlim([645 1645]); title('SARTA CO2 jac'); legend('analytic','finitediff','kcarta','location','best')
  figure(2); clf; plot(hhh.vchan,jac4,hhh.vchan,finitediff4,boo.fKc,boo.rKc(:,2));   xlim([645 1645]); title('SARTA N2O jac'); legend('analytic','finitediff','kcarta','location','best')
  figure(3); clf; plot(hhh.vchan,jac6,hhh.vchan,finitediff6,boo.fKc,boo.rKc(:,4));   xlim([645 1645]); title('SARTA CH4 jac'); legend('analytic','finitediff','kcarta','location','best')
  figure(4); clf; plot(hhh.vchan,jacST,hhh.vchan,finitediffST,boo.fKc,boo.rKc(:,8)); xlim([645 1645]); title('SARTA SKT jac'); legend('analytic','finitediff','kcarta','location','best')

  keyboard_nowindow
  figure(5); clf; plot(hhh.vchan,xjac2,'b',hhh.vchan,jac2,'r')
end
