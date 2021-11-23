addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/PLOTTER

latbins = equal_area_spherical_bands(20);
latbinsx = 0.5*(latbins(1:end-1)+latbins(2:end));

onlytropics = find(abs(latbinsx) <= 30);
extralats   = find(abs(latbinsx) <= 50);

tropics = extralats;
whos tropics

[hjunk,ha,pjunk,pa] = rtpread('/asl/rtp/rtp_airicrad_v6/clear/2002/era_airicrad_day244_clear.rtp');  

rtime0   = utc2taiSergio(2004,01,01,0.0);
rtime365 = utc2taiSergio(2004,12,31,23);

for iii = 1 : length(tropics)
  ii = tropics(iii);
  stats30 = load(['/home/strow/Work/Airs/Stability/Data/Desc/statlat' num2str(ii) '.mat']);

  [h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
  pbefore.solzen = 150 * ones(size(pbefore.stemp));
  boo = ii;
  [h,pbefore] = subset_rtp(h,pbefore,[],[],boo);
  p = pbefore;

  good = find(stats30.rtime_mean >= rtime0 & stats30.rtime_mean <= rtime365);
  p.gas_1 = nanmean(stats30.gas1_mean(good,:)',2);
  p.gas_3 = nanmean(stats30.gas3_mean(good,:)',2);
  p.ptemp = nanmean(stats30.ptemp_mean(good,:)',2);
  p.nlevs = nanmean(stats30.nlevs_mean(good)');
  p.spres = nanmean(stats30.spres_mean(good)');
  p.stemp = nanmean(stats30.stemp_mean(good)');
  p.satzen = nanmean(stats30.satzen_mean(good)');
  p.solzen = nanmean(stats30.solzen_mean(good)');
  p.robs1 = nanmean(stats30.robs(good,:)',2);
  p.rcalc = nanmean(stats30.rclr(good,:)',2);
  p.rtime = nanmean(stats30.rtime_mean(good)');

  boo = find(isnan(p.ptemp) | isnan(p.gas_1) | isnan(p.gas_3));
  p.ptemp(boo) = 300;
  p.gas_1(boo) = 0.0;
  p.gas_3(boo) = 0.0;

  h.nchan = hjunk.nchan;
  h.ichan = hjunk.ichan;
  h.vchan = hjunk.vchan;

  if iii == 1
    hall = h;
    pall = p;
  else
    [hall,pall] = cat_rtp(hall,pall,h,p);
  end

end

rtpwrite('loop_tropical_savestats_W_T_secondderiv.op.rtp',hall,ha,pall,pa);
sarta   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
sartaer =   ['!' sarta ' fin=loop_tropical_savestats_W_T_secondderiv.op.rtp fout=loop_tropical_savestats_W_T_secondderiv.rp.rtp'];
eval(sartaer)
[hall2,~,pall2,~] = rtpread('loop_tropical_savestats_W_T_secondderiv.rp.rtp');
tobs = rad2bt(hall.vchan,pall.robs1);
tcal0 = rad2bt(hall.vchan,pall.rcalc);
tcalF = rad2bt(hall.vchan,pall2.rcalc);

plot(hall2.vchan,mean(tobs'-tcal0'),'b',hall2.vchan,mean(tobs'-tcalF'),'r',...
    hall2.vchan,std(tobs'-tcal0'),'c--',hall2.vchan,std(tobs'-tcalF'),'m--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now this is like step_delta_T_WV_CO2_secondorder.m
%% and/or kcarta_2ndorderjacs_wv.m

pall3 = pall2;
pall3.gas_1 = pall3.gas_1*1.1;
rtpwrite('pert_testKC_wv.op.rtp',h,ha,pall3,pa);
sartaer = ['!' sarta ' fin=pert_testKC_wv.op.rtp fout=pert_testKC_wv.rp.rtp'];
eval(sartaer);
[hx,hax,px,pax] = rtpread('pert_testKC_wv.rp.rtp');
tx = rad2bt(hx.vchan,px.rcalc);
dq = pall3.gas_1(1:97,:)-pall2.gas_1(1:97,:);
dq = sum(dq,1);
deltaBT = tx-tcalF; %% wierd but true
plot(hx.vchan,deltaBT./(ones(2645,1)*dq));
plot(hx.vchan,mean(deltaBT./(ones(2645,1)*dq),2),'b',hx.vchan,std(deltaBT./(ones(2645,1)*dq),1,2),'c')

for ii = 1 : length(pall2.stemp)
  wv_second_deriv
  t_second_deriv
  co2_second_deriv
  stemp_second_deriv
end
savex.vchan = hall2.vchan;
savex.plevs = pnew0.plevs(:,1);

figure(1); clf
subplot(211);  plot(savex.vchan,nanmean(savex.water_stemp_true'),'b',savex.vchan,nanstd(double(savex.water_stemp_true)'),'c')
  title('true (d^2BT)/(dWV dstemp)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);
subplot(212);  plot(savex.vchan,nanmean(savex.water_stemp_log'),'b',savex.vchan,nanstd(savex.water_stemp_log'),'c')
  title(' wv*(d^2BT)/(dWV dstemp)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

figure(2); clf
subplot(211);  plot(savex.vchan,nanmean(savex.co2_ta_true'),'b',savex.vchan,nanstd(double(savex.co2_ta_true)'),'c')
  title('true (d^2BT)/(dCO2 dTa)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);
subplot(212);  plot(savex.vchan,nanmean(savex.co2_ta_log'),'b',savex.vchan,nanstd(savex.co2_ta_log'),'c')
  title(' co2*(d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

figure(3); clf
plot(savex.vchan,nanmean(savex.ta_stemp_true'),'b',savex.vchan,nanstd(double(savex.ta_stemp_true)'),'c')
  title('true (d^2BT)/(dstemp dTa)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

% save tropical_savex_savestats_W_T_secondderiv.mat savex
% save pm50_savex_savestats_W_T_secondderiv.mat savex
