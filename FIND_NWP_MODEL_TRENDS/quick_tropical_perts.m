if ~exist('iOffSet')
  addpath /home/sergio/MATLABCODE
  addpath /home/sergio/MATLABCODE/COLORMAP
  addpath /asl/matlib/h4tools
  addpath /asl/matlib/aslutil
  iOffSet = 40;
  load spectral_rate_avgs_umbc_obs_era5_merra2_airsL3_climcapsL3.mat
end

  
figure(iOffSet+40); clf;
  plot(savetherates.lfavg,savetherates.latavg); plotaxis2; xlabel('landfrac'); ylabel('Latitude');
  plot(savetherates.latavg,savetherates.lfavg); plotaxis2; ylabel('landfrac'); xlabel('Latitude');

figure(iOffSet+40); clf;
  pcolor(savetherates.fchanx,savetherates.latavg,savetherates.obs); shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1);xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640]); title('Obs Rates');
  pcolor(savetherates.fchanx,savetherates.latavg,savetherates.obs-savetherates.era5); 
    shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1/10);xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640]); title('Obs-ERA5 Rates');

boo = input('Enter \n     (1) boo = find(abs(savetherates.latavg) < 30); \n     (2) boo = find (0 < savetherates.latavg & savetherates.latavg < 30) : ');
if length(boo) == 0
  boo = 1;
end
if boo == 1
  boo = find(abs(savetherates.latavg) < 30);
elseif boo == 2
  boo = find (0 < savetherates.latavg & savetherates.latavg < 30);
end
whos boo

%% mean and std dev of differences btween obs and ERA5
plot(savetherates.fchanx,nanmean(savetherates.obs(boo,:)-savetherates.era5(boo,:),1),savetherates.fchanx,nanstd(savetherates.obs(boo,:)-savetherates.era5(boo,:),[],1))
  plotaxis2; 
  axis([600 1600 -0.01 +0.01])

%% mean obs and mean ERA5, with 0.01 K dSKt/dt offset
plot(savetherates.fchanx,nanmean(savetherates.obs(boo,:)-0.01,1),savetherates.fchanx,nanmean(savetherates.era5(boo,:),1)-0.01,'r')
  plotaxis2; title('(b) tropical obs (r) tropical ERA5')
    axis([800-15 1200+15 -0.02 +0.02]); %% window region
    axis([800-15 1250+15 -0.02 +0.02]); %% window region
    axis([800-15 1250+15 -0.04 +0.02]); %% window region
    axis([1220 1240 -0.04 +0.02]); %% window region
    axis([800 830 -0.04 +0.02]); %% window region
    axis([800 830 -0.01 +0.01]); %% window region
    axis([950 970 -0.01 +0.01]); %% window region
    axis([890 910 -0.01 +0.01]); %% window region
    axis([800 1210 -0.01 +0.01]); %% window region

%% mean obs and mean ERA5, raw
plot(savetherates.fchanx,nanmean(savetherates.obs(boo,:),1),savetherates.fchanx,nanmean(savetherates.era5(boo,:),1),'r',savetherates.fchanx,nanmean(savetherates.umbc(boo,:),1),'g')
  plotaxis2; legend('tropical obs','tropical ERA5','umbc fits','location','best','fontsize',10);
  axis([1400-15 1600+15 -0.02 +0.02]); %% shows something similar to Fig 5 in paper_jgr.pdf
  axis([600 1600 -0.03 +0.03])
  axis([800 1210 -0.02 +0.02]); %% window region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('contrast 1.001 and 0.999 changes in WV only')
[hpert,ha,ppert,pa] = rtpread('/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/regr49_1013_400ppm_unitemiss.op.rtp');
[hpert,ppert] = replicate_rtp_headprof(hpert,ppert,1,3);
ppert.gas_1(:,1) =  ppert.gas_1(:,1) * 1.0;
ppert.gas_1(:,2) =  ppert.gas_1(:,2) * 1.002;
ppert.gas_1(:,3) =  ppert.gas_1(:,3) * 0.998;
rtpwrite('pertjunk.op.rtp',hpert,ha,ppert,pa);
sartaer = ['!//home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=pertjunk.op.rtp fout=pertjunk.rp.rtp'];
eval(sartaer);
[hpert,ha,ppert,pa] = rtpread('pertjunk.rp.rtp');

tpert = rad2bt(hpert.vchan,ppert.rcalc);
figure(iOffSet+41); clf; plot(hpert.vchan,tpert(:,2)-tpert(:,1),hpert.vchan,tpert(:,3)-tpert(:,1))
  title('WV only change')
  plotaxis2;
  legend('Increased water-water0','Reduced water-water0','location','best','fontsize',10);
  xlim([800 1200]);

%%%%%%%%%%%%%%%%%%%%%%%%%

[hpert,ha,ppert,pa] = rtpread('/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/regr49_1013_400ppm_unitemiss.op.rtp');
[hpert,ppert] = replicate_rtp_headprof(hpert,ppert,1,3);
ppert.gas_1(:,1) =  ppert.gas_1(:,1) * 1.0;
ppert.gas_1(:,2) =  ppert.gas_1(:,2) * 1.002; ppert.ptemp(:,2) =  ppert.ptemp(:,2) + 0.02; ppert.stemp(2) =  ppert.stemp(2) + 0.02;
ppert.gas_1(:,3) =  ppert.gas_1(:,3) * 0.998; ppert.ptemp(:,3) =  ppert.ptemp(:,3) + 0.02; ppert.stemp(3) =  ppert.stemp(3) + 0.02;
rtpwrite('pertjunk.op.rtp',hpert,ha,ppert,pa);
sartaer = ['!//home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=pertjunk.op.rtp fout=pertjunk.rp.rtp'];
eval(sartaer);
[hpert,ha,ppert,pa] = rtpread('pertjunk.rp.rtp');

tpert = rad2bt(hpert.vchan,ppert.rcalc);
figure(iOffSet+42); clf; plot(hpert.vchan,tpert(:,2)-tpert(:,1),hpert.vchan,tpert(:,3)-tpert(:,1))
  title('WV + T change')
  plotaxis2;
  legend('Increased water-water0','Reduced water-water0','location','best','fontsize',10);
  xlim([800 1200]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('contrast 1.001 and 1.002 changes in WV, based on Table 1 in JGR paper')
paper_dSKT_dt = 0.018;
paper_dAT_dt  = 0.028;
paper_dWV_dt  = 0.002;

disp('contrast 1.001 and 1.002 changes in WV, based on Fig 12 in JGR paper')
paper_dSKT_dt = 0.018;
paper_dAT_dt  = 0.020;
paper_dWV_dt  = 0.003;

[hpert,ha,ppert,pa] = rtpread('/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/regr49_1013_400ppm_unitemiss.op.rtp');
[hpert,ppert] = replicate_rtp_headprof(hpert,ppert,1,3);
ppert.gas_1(:,1) =  ppert.gas_1(:,1) * 1.0;
ppert.gas_1(:,2) =  ppert.gas_1(:,2) * (1 + 1*paper_dWV_dt);
ppert.gas_1(:,3) =  ppert.gas_1(:,3) * (1 + 2*paper_dWV_dt);
ppert.gas_2(:,2) =  ppert.gas_2(:,2) * (1 + 2.2/385);
ppert.gas_2(:,3) =  ppert.gas_2(:,3) * (1 + 2.2/385);
rtpwrite('pertjunk.op.rtp',hpert,ha,ppert,pa);
sartaer = ['!//home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=pertjunk.op.rtp fout=pertjunk.rp.rtp'];
eval(sartaer);
[hpert,ha,ppert,pa] = rtpread('pertjunk.rp.rtp');

tpert = rad2bt(hpert.vchan,ppert.rcalc);
figure(iOffSet+41); clf; plot(hpert.vchan,tpert(:,2)-tpert(:,1),hpert.vchan,tpert(:,3)-tpert(:,1))
  title('WV only change')
  plotaxis2;
  legend('1.002 Increased water-water0','1.004 Increased water-water0','location','best','fontsize',10);
  xlim([800 1200]);

%%%%%%%%%%%%%%%%%%%%%%%%%

[hpert,ha,ppert,pa] = rtpread('/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/regr49_1013_400ppm_unitemiss.op.rtp');
[hpert,ppert] = replicate_rtp_headprof(hpert,ppert,1,4);
ppert.gas_1(:,1) =  ppert.gas_1(:,1) * 1.0;
ppert.gas_1(:,2) =  ppert.gas_1(:,2) * (1 + 0.5*paper_dWV_dt); ppert.ptemp(:,2) =  ppert.ptemp(:,2) + paper_dAT_dt; ppert.stemp(2) =  ppert.stemp(2) + paper_dSKT_dt;
ppert.gas_1(:,3) =  ppert.gas_1(:,3) * (1 + 1.0*paper_dWV_dt); ppert.ptemp(:,3) =  ppert.ptemp(:,3) + paper_dAT_dt; ppert.stemp(3) =  ppert.stemp(3) + paper_dSKT_dt;
ppert.gas_1(:,4) =  ppert.gas_1(:,4) * (1 + 2.0*paper_dWV_dt); ppert.ptemp(:,4) =  ppert.ptemp(:,4) + paper_dAT_dt; ppert.stemp(4) =  ppert.stemp(4) + paper_dSKT_dt;
ppert.gas_2(:,2) =  ppert.gas_2(:,2) * (1 + 2.2/385);
ppert.gas_2(:,3) =  ppert.gas_2(:,3) * (1 + 2.2/385);
ppert.gas_2(:,4) =  ppert.gas_2(:,4) * (1 + 2.2/385);
rtpwrite('pertjunk.op.rtp',hpert,ha,ppert,pa);
sartaer = ['!//home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=pertjunk.op.rtp fout=pertjunk.rp.rtp'];
eval(sartaer);
[hpert,ha,ppert,pa] = rtpread('pertjunk.rp.rtp');

tpert = rad2bt(hpert.vchan,ppert.rcalc);
figure(iOffSet+42); clf; plot(hpert.vchan,tpert(:,2)-tpert(:,1),hpert.vchan,tpert(:,3)-tpert(:,1),hpert.vchan,tpert(:,4)-tpert(:,1),...
                              savetherates.fchanx,nanmean(savetherates.obs(boo,:),1),'k')
  title('WV + T change')
  plotaxis2;
  legend('0.5 Increased water-water0','1.0 Increased water-water0','2.0 Increased water-water0','Tropical Obs','location','best','fontsize',10);
  xlim([800 1200]);
  ylim([-1 +1]*0.02)

