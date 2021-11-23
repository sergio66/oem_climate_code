addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
addpath /home/sergio/MATLABCODE
addpath /asl/matlib/rtptools

clear all
[h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002.rtp');
for ii = 1 : 64
  ind = (1:72)+(ii-1)*72;
  [hx,px] = subset_rtp_allcloudfields(h1_4608,p1_4608,[],[],ind);
  px = find_average_rtp(px);
  if ii == 1
    havg64 = hx;
    pavg64 = px;
  else
   [havg64,pavg64] = cat_rtp(havg64,pavg64,hx,px);
  end
end

rtpwrite('pavg64cld.op.rtp',havg64,[],pavg64,[]);
error('ksjgskgs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
[h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
for ii = 1 : 64
  ind = (1:72)+(ii-1)*72;
  [hx,px] = subset_rtp_allcloudfields(h1_4608,p1_4608,[],[],ind);
  px = find_average_rtp(px);
  if ii == 1
    havg64 = hx;
    pavg64 = px;
  else
   [havg64,pavg64] = cat_rtp(havg64,pavg64,hx,px);
  end
end

rtpwrite('pavg64.op.rtp',havg64,[],pavg64,[]);
