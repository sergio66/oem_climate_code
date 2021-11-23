addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op.rtp');

[hnew,pnew] = replicate_rtp_headprof(h,p,20,45);    %% replicate(h,p,iRTP,numtimes);

[hx,hax,px,pax] = rtpread('/asl/rtp/rtp_airicrad_v6/allfov/2020/004/allfov_ecmwf_airicrad_d2020004_061.rtp');
pnew.scanang = px.scanang(1:45);
pnew.satzen  = px.satzen(1:45);
pnew.zobs    = px.zobs(1:45);

rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin20_45angles.op.rtp',hnew,ha,pnew,pa)


