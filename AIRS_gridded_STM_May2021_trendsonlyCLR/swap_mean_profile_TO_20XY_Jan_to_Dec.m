function [hout,pavg] = swap_mean_profile_TO_20XY_08_15(hIN,pIN,YYYY,HHHH);

addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA/
addpath /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB
addpath /home/sergio/MATLABCODE/matlib/clouds
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta

run_sarta.clear = -1;
run_sarta.cloud = -1;
run_sarta.cumsum = 9999;
code1 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
run_sarta.sartacloud_code = code1;
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

rtime  = zeros(001,4608);
stemp  = zeros(001,4608);
ptemp  = zeros(101,4608);
gas_1  = zeros(101,4608);
gas_2  = zeros(101,4608);
gas_3  = zeros(101,4608);
gas_4  = zeros(101,4608);
gas_5  = zeros(101,4608);
gas_6  = zeros(101,4608);
gas_9  = zeros(101,4608);
gas_12 = zeros(101,4608);

hax = struct;
pax = struct;
[~,hax,~,pax] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2003/10/01/interp_analysis_cloudy_airs_l1c_era_sarta_baum_ice_timeoffset_0000_2003.10.01.240_cumsum_-1.rtp');

for mm = 1 : 12
  hx = hIN; px = pIN;

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

  px.rtime = utc2taiSergio(YYYY,mm,15,HHHH) * ones(size(px.rtime));
      
  [hx,hax,px,pax] = make_generic_interp_ERA_rtp(hx,hax,px,pax);
    
  rtpwrite('junk.ip.rtp',hx,hax,px,pax);
  klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp >& ugh'];
  eval(klayerser);
  
  [hout,~,pout,~] = rtpread('junk.op.rtp');
  str = ['pout' num2str(mm,'%02d') ' = pout;'];
  eval(str);

  rtime = rtime + pout.rtime;
  stemp = stemp + pout.stemp;
  ptemp = ptemp + pout.ptemp;
  gas_1 = gas_1 + pout.gas_1;
  gas_2 = gas_2 + pout.gas_2;
  gas_3 = gas_3 + pout.gas_3;
  gas_4 = gas_4 + pout.gas_4;
  gas_5 = gas_5 + pout.gas_5;
  gas_6 = gas_6 + pout.gas_6;
  gas_9 = gas_9 + pout.gas_9;
  gas_12= gas_12+ pout.gas_12;

end

pavg = pout;
pavg.rtime = rtime/12;
pavg.stemp = stemp/12;
pavg.ptemp = ptemp/12;
pavg.gas_1 = gas_1/12;
pavg.gas_2 = gas_2/12;
pavg.gas_3 = gas_3/12;
pavg.gas_4 = gas_4/12;
pavg.gas_5 = gas_5/12;
pavg.gas_6 = gas_6/12;
pavg.gas_9 = gas_9/12;
pavg.gas_12= gas_12/12;

eval(['!/bin/rm junk.ip.rtp junk.op.rtp ugh'])
