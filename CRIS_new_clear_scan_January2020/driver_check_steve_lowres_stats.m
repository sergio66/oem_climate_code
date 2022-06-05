clear all

%% see clust_make_sergio_anomaly_nonucal_2012.m
ii = 20;
finX = ['/home/strow/Work/Cris/Stability/Data/Desc/statlat' num2str(ii) '.mat'];
aX = load(finX);
robsX = squeeze(aX.robs(:,5,:));
rclrX  = squeeze(aX.rclr(:,5,:));
rtimeX = aX.rtime_mean(:,5);
  rtimeX = rtimeX';
  robsX = robsX';
  rclrX = rclrX';

inst = 'airs';
inst = 'al1c';
inst = 'cris';  %% this is automatically CIRS NSR

%% load_fairs;  %% load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/airs_f.mat'); basically give f2378 and f2645
usethese0 = 1 : 1317;

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/cris_ichan_vchan_nsr_fsr.mat
usethese = find(nsr.ichan <= 1305);
fairsX = nsr.vchan(usethese);

[yyx,mmx,ddx] = tai2utcSergio(rtimeX);
days_since_2012x = change2days(yyx,mmx,ddx,2012);
ttimeX = yyx + (mmx-1)/12 + (ddx-1)/12/30;
plot(ttimeX(1:end-1),diff(ttimeX))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtime = [];
robs = [];
rclr = [];

for ii = 2012 : 2021
  ii
  fname = ['//asl/stats/cris/npp/clear/lowres/rtp_cris_lowres_rad_' num2str(ii) '_clear_desc_ocean.mat'];
  a = load(fname);
  rtime = [rtime; squeeze(a.rtime_mean(:,20,5))];
  robs = [robs; squeeze(a.robs(:,20,5,:))];
  rclr = [rclr; squeeze(a.rclr(:,20,5,:))];
end

robs = robs';
rclr = rclr';

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/cris_ichan_vchan_nsr_fsr.mat
fcris = nsr.vchan;

[yyx,mmx,ddx] = tai2utcSergio(rtime);
days_since_2012x = change2days(yyx,mmx,ddx,2012);
ttime = yyx + (mmx-1)/12 + (ddx-1)/12/30;
plot(ttime(1:end-1),diff(ttime))

i2355 = find(fcris >= 2355,1);
plot(ttime,rad2bt(2355,robs(i2355,:)))
plot(ttime,rad2bt(2355,robs(i2355,:))-rad2bt(2355,rclr(i2355,:)));  title('BT2355 : obs - calc')

i900 = find(fcris >= 900,1);
plot(ttime,rad2bt(900,robs(i900,:)))
plot(ttime,rad2bt(900,robs(i900,:))-rad2bt(900,rclr(i900,:))); title('BT900 : obs - calc')

i1419 = find(fcris >= 1419,1);
plot(ttime,rad2bt(1419,robs(i1419,:)))
plot(ttime,rad2bt(1419,robs(i1419,:))-rad2bt(1419,rclr(i1419,:))); title('BT1419 : obs-calc')

tobs = rad2bt(fcris,robs);
tclr = rad2bt(fcris,rclr);
plot(fcris,nanmean(tobs'-tclr'),'b',fcris,nanstd(tobs'-tclr'),'c')
 plotaxis2; ylim([-2 +2]); legend('bias','std','location','best')

plot(ttime,rad2bt(900,robs(i900,:)),ttimeX,rad2bt(900,robsX(i900,:)))
plot(ttime,rad2bt(900,robs(i900,:))-rad2bt(900,robsX(i900,:)))
plot(ttime,rad2bt(900,rclr(i900,:))-rad2bt(900,rclrX(i900,:)))

plot(ttime,rad2bt(900,robs(i900,:))-rad2bt(900,robsX(i900,:)),ttime,rad2bt(900,rclr(i900,:))-rad2bt(900,rclrX(i900,:))); 
  title('BT900 new file-old file (b) obs (r)sarta')
plot(ttime,rad2bt(900,robsX(i900,:))-rad2bt(900,rclrX(i900,:)),'b',ttime,rad2bt(900,robs(i900,:))-rad2bt(900,rclr(i900,:)),'r')
   title('BT900 bias (b) old (r) new')

plot(ttime,rad2bt(1419,robs(i1419,:))-rad2bt(1419,robsX(i1419,:)),ttime,rad2bt(1419,rclr(i1419,:))-rad2bt(1419,rclrX(i1419,:))); 
   title('BT1419 new file-old file (b) obs (r)sarta')
plot(ttime,rad2bt(1419,robsX(i1419,:))-rad2bt(1419,rclrX(i1419,:)),'b',ttime,rad2bt(1419,robs(i1419,:))-rad2bt(1419,rclr(i1419,:)),'r')
   title('BT1419 bias (b) old (r) new')

plot(ttime,rad2bt(2355,robs(i2355,:))-rad2bt(2355,robsX(i2355,:)),ttime,rad2bt(2355,rclr(i2355,:))-rad2bt(2355,rclrX(i2355,:))); 
   title('BT2355 new file-old file (b) obs (r)sarta')
plot(ttime,rad2bt(2355,robsX(i2355,:))-rad2bt(2355,rclrX(i2355,:)),'b',ttime,rad2bt(2355,robs(i2355,:))-rad2bt(2355,rclr(i2355,:)),'r')
   title('BT2355 bias (b) old (r) new')
