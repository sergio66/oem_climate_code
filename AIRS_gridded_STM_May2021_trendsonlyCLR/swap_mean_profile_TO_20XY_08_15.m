function [hout,pout] = swap_mean_profile_TO_20XY_08_15(hIN,pIN,YYYY,HHHH);

hx = hIN; px = pIN;

px.rtime = utc2taiSergio(YYYY,08,15,HHHH) * ones(size(px.rtime));

hx.ngas = 0;
hx = rmfield(hx,'glist');
hx = rmfield(hx,'gunit');
hx.ptype = 0;

px = rmfield(px,'nlevs');
px = rmfield(px,'plevs');
px = rmfield(px,'stemp');
px = rmfield(px,'ptemp');
px = rmfield(px,'gas_1');
px = rmfield(px,'gas_2');
px = rmfield(px,'gas_3');
px = rmfield(px,'gas_4');
px = rmfield(px,'gas_5');
px = rmfield(px,'gas_6');
px = rmfield(px,'gas_9');
px = rmfield(px,'gas_12');

addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

hax = struct;
pax = struct;
[~,hax,~,pax] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2003/10/01/interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice_timeoffset_0000_2003.10.01.240_cumsum_-1.rtp');

run_sarta.clear = -1;
run_sarta.cloud = -1;
run_sarta.cumsum = 9999;
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
run_sarta.sartacloud_code = code1;

[hx,hax,px,pax] = make_generic_interp_ERA_rtp(hx,hax,px,pax);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

rtpwrite('junk.ip.rtp',hx,hax,px,pax);
klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp >& ugh'];
eval(klayerser);

[hout,~,pout,~] = rtpread('junk.op.rtp');

eval(['!/bin/rm junk.ip.rtp junk.op.rtp ugh'])
