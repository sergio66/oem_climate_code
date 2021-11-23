%{
 /asl/rtp/rtp_airicrad_v6/clear/   has rtp files
 /asl/data/stats/airs/clear        has stats mat files
%}


addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/PLOTTER

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
iaTropics = find(abs(latbins) <= 30);

[h,ha,p,pa] = rtpread('/asl/rtp/rtp_airicrad_v6/clear/2016/era_airicrad_day061_clear.rtp');

iaX = find(p.solzen > 90 & p.landfrac == 0 & abs(p.rlat) <= 30); whos iaX
iaX = find(p.solzen > 90 & p.landfrac == 0 & p.rlat > latbins(20) & p.rlat <= latbins(21)); whos iaX

tobs = rad2bt(h.vchan,p.robs1(:,iaX));
figure(1); plot(h.vchan,mean(tobs'))
figure(2); plot(h.vchan,std(tobs'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = load('/asl/data/stats/airs/clear/rtp_airicrad_era_rad_kl_2016_clear_desc.mat');
figure(1); plot(h.vchan,mean(tobs'),'b.-',h.vchan,rad2bt(h.vchan,squeeze(stats.robs(61,20,:))),'r')
figure(2); plot(h.vchan,mean(tobs')'-rad2bt(h.vchan,squeeze(stats.robs(61,20,:))),'k')

wahP = rad2bt(h.vchan,squeeze(stats.robs(61,20,:)) + squeeze(stats.rclrbias_std(61,20,:))) - rad2bt(h.vchan,squeeze(stats.robs(61,20,:)));
wahM = rad2bt(h.vchan,squeeze(stats.robs(61,20,:)) - squeeze(stats.rclrbias_std(61,20,:))) - rad2bt(h.vchan,squeeze(stats.robs(61,20,:)));
figure(3); plot(h.vchan,std(tobs'),'b',h.vchan,wahP,'k',h.vchan,abs(wahM),'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now run driver_run_quick_stats.m')
