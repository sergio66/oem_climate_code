addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/PLOTTER

latbins = equal_area_spherical_bands(20);

load('../../f2645.mat');
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

p.rcalc = p.rclr;
[h30,p30] = subset_rtp(h,p,[],[],boo30);

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
%pbefore.robs1 = anomaly30obs.all_b(:,1);
%pbefore.rcalc  = anomaly30clr.all_b(:,1);
%pbefore.robs1X = nanmean(p.robs1(:,boo30),2);
%pbefore.rcalcX = nanmean(p.rclr(:,boo30),2);
pbefore.robs1 = nanmean(p.robs1(:,boo30),2);
pbefore.rcalc = nanmean(p.rclr(:,boo30),2);
pbefore = rmfield(pbefore,'gas_2');
pbefore = rmfield(pbefore,'gas_4');
pbefore = rmfield(pbefore,'gas_5');
pbefore = rmfield(pbefore,'gas_6');
pbefore = rmfield(pbefore,'gas_9');
pbefore = rmfield(pbefore,'gas_12');

sarta   = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

rtpwrite('run_bin30_day244_385ppm.ip.rtp',h30,ha,pbefore,pa);
klayerser = ['!' klayers ' fin=run_bin30_day244_385ppm.ip.rtp fout=run_bin30_day244_385ppm.op.rtp'];
sartaer   = ['!' sarta '   fin=run_bin30_day244_385ppm.op.rtp fout=run_bin30_day244_385ppm.rp.rtp'];
eval(klayerser)
eval(sartaer)
[hx,hax,px385,pax] = rtpread('run_bin30_day244_385ppm.rp.rtp');
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(hx,px385,1,2);

rtpwrite('run_bin30_day244_385ppm_latbin30_manyprofs.ip.rtp',h30,ha,p30,pa);
klayerser = ['!' klayers ' fin=run_bin30_day244_385ppm_latbin30_manyprofs.ip.rtp fout=run_bin30_day244_385ppm_latbin30_manyprofs.op.rtp'];
sartaer   = ['!' sarta '   fin=run_bin30_day244_385ppm_latbin30_manyprofs.op.rtp fout=run_bin30_day244_385ppm_latbin30_manyprofs.rp.rtp'];
eval(klayerser)
eval(sartaer)
[hx,hax,px385_many,pax] = rtpread('run_bin30_day244_385ppm_latbin30_manyprofs.rp.rtp');
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(hx,px385,1,2);

figure(1); 
  plot(hbefore.vchan,rad2bt(hbefore.vchan,pbefore.robs1)-rad2bt(hbefore.vchan,pbefore.rcalc),'b.-',...
       f2645,rad2bt(f2645,anomaly30obs.all_b(:,1))-rad2bt(f2645,anomaly30clr.all_b(:,1)),'k',...
       hx.vchan,nanmean(rad2bt(hx.vchan,px385_many.robs1)'-rad2bt(hx.vchan,px385_many.rcalc)'),'g',...
       hx.vchan,rad2bt(hx.vchan,px385.robs1)-rad2bt(hx.vchan,px385.rcalc),'r')
hl = legend('rtp <obs-cal>','statsfile (obs-cal)','<SARTA(prof)>','SARTA(<prof>)','location','best');
set(hl,'fontsize',10)

figure(2); 
  plot(hbefore.vchan,rad2bt(hbefore.vchan,pbefore.robs1)-rad2bt(f2645,anomaly30obs.all_b(:,1)),'b',...
       f2645,rad2bt(hbefore.vchan,pbefore.rcalc)-rad2bt(f2645,anomaly30clr.all_b(:,1)),'c',...
       hx.vchan,nanmean(rad2bt(hx.vchan,px385_many.robs1)')-rad2bt(hx.vchan,px385.robs1)','g',...
       hx.vchan,nanmean(rad2bt(hx.vchan,px385_many.rcalc)')-rad2bt(hx.vchan,px385.rcalc)','r')
hl = legend('rtp <obs>-stats obs','rtp <cal>- stats cal)','obs <SARTA(prof)>-SARTA(<prof>)','cal <SARTA(prof)>-SARTA(<prof>)','location','best');
set(hl,'fontsize',10)
