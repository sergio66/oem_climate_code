addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/PLOTTER

latbins = equal_area_spherical_bands(20);

iLatbin = 30;
    fname = ['/home/strow/Work/Airs/Stability/Data_old/Desc_fits/fit_robs_lat' num2str(iLatbin) '.mat'];
    fname = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_robs_lat' num2str(iLatbin) '.mat'];
anomaly30obs = load(fname);
    fname = ['/home/strow/Work/Airs/Stability/Data_old/Desc_fits/fit_rclr_lat' num2str(iLatbin) '.mat'];
    fname = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_rclr_lat' num2str(iLatbin) '.mat'];
anomaly30clr = load(fname);

if ~exist('p')
  [h,ha,p,pa] = rtpread('/asl/rtp/rtp_airicrad_v6/clear/2002/era_airicrad_day244_clear.rtp');
end

load('../../f2645.mat');
figure(1); plot(f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))); title('robs')
figure(2); plot(f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,anomaly30clr.all_b(:,1))); title('robs-rclr')

boo30 = find(p.rlat >= latbins(30) & p.rlat <= latbins(31) & p.solzen > 90 & p.landfrac == 0);

[hbefore,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
[hbefore,pbefore] = subset_rtp(hbefore,pbefore,[],[],30);
hbefore.ptype = 0;
hbefore.ngas = 2;
hbefore.gunit = h.gunit;
hbefore.glist = h.glist;
hbefore.nchan = h.nchan;
hbefore.ichan = h.ichan;
hbefore.vchan = h.vchan;
pbefore.nlevs = 60;
pbefore.gas_1 = nanmean(p.gas_1(:,boo30),2);
pbefore.gas_3 = nanmean(p.gas_3(:,boo30),2);
pbefore.plevs = nanmean(p.plevs(:,boo30),2);
pbefore.ptemp = nanmean(p.ptemp(:,boo30),2);
pbefore.stemp = nanmean(p.stemp(boo30));
pbefore.spres = nanmean(p.spres(boo30));
pbefore.salti = nanmean(p.salti(boo30));
pbefore.scanang = nanmean(abs(p.scanang(boo30)));
pbefore.satzen = nanmean(abs(p.satzen(boo30)));
pbefore.solzen = nanmean(abs(p.solzen(boo30)));
pbefore.robs1 = anomaly30obs.all_b(:,1);
pbefore.rcalc  = anomaly30clr.all_b(:,1);
pbefore.robs1X = nanmean(p.robs1(:,boo30),2);
pbefore.rcalcX = nanmean(p.rclr(:,boo30),2);
pbefore = rmfield(pbefore,'gas_2');
pbefore = rmfield(pbefore,'gas_4');
pbefore = rmfield(pbefore,'gas_5');
pbefore = rmfield(pbefore,'gas_6');
pbefore = rmfield(pbefore,'gas_9');
pbefore = rmfield(pbefore,'gas_12');

sarta   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

rtpwrite('run_bin30_day244_385ppm.ip.rtp',hbefore,ha,pbefore,pa);
klayerser = ['!' klayers ' fin=run_bin30_day244_385ppm.ip.rtp fout=run_bin30_day244_385ppm.op.rtp'];
sartaer   = ['!' sarta '   fin=run_bin30_day244_385ppm.op.rtp fout=run_bin30_day244_385ppm.rp.rtp'];
eval(klayerser)
eval(sartaer)
[hx,hax,px385,pax] = rtpread('run_bin30_day244_385ppm.rp.rtp');
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(hx,px385,1,2);

px370 = px385;
px370.gas_2 = px370.gas_2 * 370/385;
rtpwrite('run_bin30_day244_370ppm.op.rtp',hx,ha,px370,pa);
sartaer   = ['!' sarta '   fin=run_bin30_day244_370ppm.op.rtp fout=run_bin30_day244_370ppm.rp.rtp'];
eval(sartaer)
[hx,hax,px370,pax] = rtpread('run_bin30_day244_370ppm.rp.rtp');

figure(2); plot(f2645,rad2bt(f2645,anomaly30obs.all_b(:,1)),'b',f2645,rad2bt(f2645,anomaly30clr.all_b(:,1)),'c',...
                f2645,rad2bt(f2645,px385.rcalc),'g',f2645,rad2bt(f2645,px370.rcalc),'r'); 
                title('robs and rclr')

figure(3); plot(f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,anomaly30clr.all_b(:,1)),'b',...
                f2645,rad2bt(f2645,pbefore.robs1X)-rad2bt(f2645,pbefore.rcalcX),'k',... 
                f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,px385.rcalc),'g',...
                f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,px370.rcalc),'r')
                title('robs-rclr')
hl = legend('stats allbobs-allbcal','mean rtp(obs-cal)','mean prof-->sarta 385ppm','mean prof-->sarta 370ppm','location','best');
set(hl,'fontsize',10);

figure(4); plot(f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,pbefore.robs1X),'b',...
                f2645,rad2bt(f2645,anomaly30clr.all_b(:,1))-rad2bt(f2645,pbefore.rcalcX),'k')
hl = legend('stats allbobs-<obs>','stats allcal-<clr>','location','best');
set(hl,'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

