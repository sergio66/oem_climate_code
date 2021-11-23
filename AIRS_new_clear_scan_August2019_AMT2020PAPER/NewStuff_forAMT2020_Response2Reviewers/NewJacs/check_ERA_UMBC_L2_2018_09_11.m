g175_L2 = load('/asl/s1/sergio/rtp/rtp_airsL2/2018/254/gran175.mat');

p175 = load('/asl/s1/sergio/rtp/singlefootprintretrievals_airibrad_v5/2018/09/11/ecm_retr175_-1_iDET_4_iStemp_ColWV_5151_iCenterFov_-1_iBadDCC_atrack_-200.mat');

[h,ha,p,pa] = rtpread('/asl/s1/sergio/rtp/rtp_airibrad_v5/2018/09/11/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.09.11.175.rtp'); %% levels
[h,ha,p,pa] = rtpread('/asl/s1/sergio/rtp/rtp_airibrad_v5/2018/09/11/cloudy_airs_l1b_ecm_sarta_baum_ice.2018.09.11.175.rp.rtp'); %% layers

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
if h.ptype == 0
  rtpwrite('junkyjunk.ip.rtp',h,ha,p,pa);
  klayerser = ['!' klayers ' fin=junkyjunk.ip.rtp fout=junkyjunk.op.rtp'];
  eval(klayerser);
  [h,ha,p,pa] = rtpread('junkyjunk.op.rtp');
  rmer = ['!/bin/rm junkyjunk.ip.rtp junkyjunk.op.rtp'];
  eval(rmer);
end

l2_850 = find(g175_L2.gprof.plevs >= 850,1);
l2_500 = find(g175_L2.gprof.plevs >= 500,1);
l2_250 = find(g175_L2.gprof.plevs >= 250,1);

era_850 = find(p.plevs(:,1) >= 850,1);
era_500 = find(p.plevs(:,1) >= 500,1);
era_250 = find(p.plevs(:,1) >= 250,1);

[l2_850 l2_500 l2_250;era_850 era_500 era_250]

plot(g175_L2.gprof.ptemp(l2_850,:)-p.ptemp(era_850,:))

inumL2 = length(l2_good);
inumP = length(latbins_results(8).stemp_era);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you are running driver_process_cutup_16years_clear.m
dbt = -5:0.25:5;
l2_good = find(g175_L2.gprof.Qual_Clear_OLR <= 1 & p.landfrac == 0);
l2_good = find(g175_L2.gprof.Qual_Clear_OLR == 0 & g175_L2.gprof.Qual_TSurf == 0 & p.landfrac == 0);

plot(dbt,histc(g175_L2.gprof.ptemp(l2_850,l2_good)-p.ptemp(era_850,l2_good),dbt),dbt,histc(latbins_results(8).tz_umbc(1,:)-latbins_results(8).tz_era(1,:),dbt))
  title('850 mb T diff between ERA and X'); hl = legend('L2-ERA','UMBC-ERA','location','best','fontsize',10); grid

plot(dbt,histc(g175_L2.gprof.ptemp(l2_500,l2_good)-p.ptemp(era_500,l2_good),dbt),dbt,histc(latbins_results(8).tz_umbc(2,:)-latbins_results(8).tz_era(2,:),dbt))
  title('500 mb T diff between ERA and X'); hl = legend('L2-ERA','UMBC-ERA','location','best','fontsize',10); grid

plot(dbt,histc(g175_L2.gprof.ptemp(l2_250,l2_good)-p.ptemp(era_250,l2_good),dbt),dbt,histc(latbins_results(8).tz_umbc(2,:)-latbins_results(8).tz_era(2,:),dbt))
  title('250 mb T diff between ERA and X'); hl = legend('L2-ERA','UMBC-ERA','location','best','fontsize',10); grid

plot(dbt,histc(g175_L2.gprof.stemp(l2_good)-p.stemp(l2_good),dbt),'g',dbt,histc(latbins_results(8).stemp_umbc-latbins_results(8).stemp_era,dbt),'g--',...
     dbt,histc(g175_L2.gprof.ptemp(l2_850,l2_good)-p.ptemp(era_850,l2_good),dbt),'k',dbt,histc(latbins_results(8).tz_umbc(1,:)-latbins_results(8).tz_era(1,:),dbt),'k--',...
     dbt,histc(g175_L2.gprof.ptemp(l2_500,l2_good)-p.ptemp(era_500,l2_good),dbt),'b',dbt,histc(latbins_results(8).tz_umbc(2,:)-latbins_results(8).tz_era(2,:),dbt),'b--',...
     dbt,histc(g175_L2.gprof.ptemp(l2_250,l2_good)-p.ptemp(era_250,l2_good),dbt),'r',dbt,histc(latbins_results(8).tz_umbc(3,:)-latbins_results(8).tz_era(3,:),dbt),'r--','linewidth',2)
title('T diff between X and ERA'); hl = legend('STEMP L2-ERA','STEMP UMBC-ERA','850 L2-ERA','850 UMBC-ERA','500 L2-ERA','500 UMBC-ERA','250 L2-ERA','250 UMBC-ERA','location','best','fontsize',10); grid

%plot(g175_L2.gprof.ptemp(:,l2_good),g175_L2.gprof.plevs(:,l2_good),'b',p.ptemp(:,l2_good),p.plevs(:,l2_good),'r');
plot(nanmean(g175_L2.gprof.ptemp(:,l2_good)'),nanmean(g175_L2.gprof.plevs(:,l2_good)'),'b',nanmean(p.ptemp(:,l2_good)'),nanmean(p.plevs(:,l2_good)'),'r');
axis([180 320 0 1000]); set(gca,'ydir','reverse')

plot(nanmean(g175_L2.gprof.plevs(1:97,l2_good)'-p.plevs(1:97,l2_good)'))
plot(nanmean(g175_L2.gprof.ptemp(1:97,l2_good)'-p.ptemp(1:97,l2_good)'),nanmean(g175_L2.gprof.plevs(1:97,l2_good)'),'b',...
     nanstd(g175_L2.gprof.ptemp(1:97,l2_good)'-p.ptemp(1:97,l2_good)'),nanmean(g175_L2.gprof.plevs(1:97,l2_good)'),'c--')
plotaxis2; set(gca,'ydir','reverse'); axis([-2 +2 0 1000])
%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 12150
  junk = find(g175_L2.gprof.plevs(:,ii) >= 850,1); g175_L2.gprof.ptemp3levs(1,ii) = g175_L2.gprof.ptemp(junk,ii); g175_L2.gprof.gas_13levs(1,ii) = g175_L2.gprof.gas_1(junk,ii); 
  junk = find(g175_L2.gprof.plevs(:,ii) >= 500,1); g175_L2.gprof.ptemp3levs(2,ii) = g175_L2.gprof.ptemp(junk,ii); g175_L2.gprof.gas_13levs(2,ii) = g175_L2.gprof.gas_1(junk,ii); 
  junk = find(g175_L2.gprof.plevs(:,ii) >= 250,1); g175_L2.gprof.ptemp3levs(3,ii) = g175_L2.gprof.ptemp(junk,ii); g175_L2.gprof.gas_13levs(3,ii) = g175_L2.gprof.gas_1(junk,ii); 

  junk = find(p.plevs(:,ii) >= 850,1); p.ptemp3levs(1,ii) = p.ptemp(junk,ii); p.gas_13levs(1,ii) = p.gas_1(junk,ii);
  junk = find(p.plevs(:,ii) >= 500,1); p.ptemp3levs(2,ii) = p.ptemp(junk,ii); p.gas_13levs(2,ii) = p.gas_1(junk,ii);
  junk = find(p.plevs(:,ii) >= 250,1); p.ptemp3levs(3,ii) = p.ptemp(junk,ii); p.gas_13levs(3,ii) = p.gas_1(junk,ii);
end

figure(1)
plot(dbt,histc(g175_L2.gprof.stemp(l2_good)-p.stemp(l2_good),dbt),'g',dbt,histc(latbins_results(8).stemp_umbc-latbins_results(8).stemp_era,dbt),'g--',...
     dbt,histc(g175_L2.gprof.ptemp3levs(1,l2_good)-p.ptemp3levs(1,l2_good),dbt),'k',dbt,histc(latbins_results(8).tz_umbc(1,:)-latbins_results(8).tz_era(1,:),dbt),'k--',...
     dbt,histc(g175_L2.gprof.ptemp3levs(2,l2_good)-p.ptemp3levs(2,l2_good),dbt),'b',dbt,histc(latbins_results(8).tz_umbc(2,:)-latbins_results(8).tz_era(2,:),dbt),'b--',...
     dbt,histc(g175_L2.gprof.ptemp3levs(3,l2_good)-p.ptemp3levs(3,l2_good),dbt),'r',dbt,histc(latbins_results(8).tz_umbc(3,:)-latbins_results(8).tz_era(3,:),dbt),'r--','linewidth',2)
title('T diff between X and ERA'); hl = legend('STEMP L2-ERA','STEMP UMBC-ERA','850 L2-ERA','850 UMBC-ERA','500 L2-ERA','500 UMBC-ERA','250 L2-ERA','250 UMBC-ERA','location','best','fontsize',10); grid
plotaxis2;

figure(2)
dq = 0 : 0.01 : 2;
plot(dq,histc(g175_L2.gprof.gas_13levs(1,l2_good)./p.gas_13levs(1,l2_good),dq)/inumL2,'k',dq,histc(latbins_results(8).wv_umbc(1,:)./latbins_results(8).wv_era(1,:),dq)/inumP,'k--',...
     dq,histc(g175_L2.gprof.gas_13levs(2,l2_good)./p.gas_13levs(2,l2_good),dq)/inumL2,'b',dq,histc(latbins_results(8).wv_umbc(2,:)./latbins_results(8).wv_era(2,:),dq)/inumP,'b--',...
     dq,histc(g175_L2.gprof.gas_13levs(3,l2_good)./p.gas_13levs(3,l2_good),dq)/inumL2,'r',dq,histc(latbins_results(8).wv_umbc(3,:)./latbins_results(8).wv_era(3,:),dq)/inumP,'r--','linewidth',2)
title('WV ratio between X and ERA'); hl = legend('850 L2/ERA','850 UMBC/ERA','500 L2/ERA','500 UMBC/ERA','250 L2/ERA','250 UMBC/ERA','location','best','fontsize',10); grid
plotaxis2;

