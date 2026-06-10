umbc100_spectral_olr = struct;    %% so it has no fields

choose100 = 1:100;
choose100 = 1:length(p.stemp);

if length(choose100 ) ~= length(p.stemp)
  [htest,ptest] = subset_rtp_allcloudfields(h,p,[],[],choose100);
else
  htest = h;
  ptest = p;
end

if htest.pfields == 1 | ~isfield(ptest,'rcalc')
  disp('making ptest.rcalc ....')
  pcopy = ptest;
  hcopy = htest;      
  klayers = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_80km/klayersV205_0_80km/BinV221/klayers_airs';
  sarta   = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_apr26_H2024';
  rtpwrite('olrjunk.op.rtp',htest,ha,ptest,pa);
  sartaer = ['!time ' sarta ' fin=olrjunk.op.rtp fout=olrjunk.rp.rtp'];
  eval(sartaer)
  [htest,~,ptest,~] = rtpread('olrjunk.rp.rtp');
  rmer = ['!/bin/rm olrjunk.op.rtp olrjunk.rp.rtp'];
  eval(rmer)
end
if length(choose100 ) ~= length(p.stemp)
  umbc100_spectral_olr = compute_feedbacks_generic_ecRad(htest,ptest,results(choose100,:),results(choose100,6)',deltaT(:,choose100),fracWV(:,choose100),fracO3(:,choose100),umbc100_spectral_olr,-1,rlat65,rlon73,-1,'UMBC');
else
  umbc100_spectral_olr = compute_feedbacks_generic_ecRad(htest,ptest,results,results(:,6)',deltaT,fracWV,fracO3,umbc100_spectral_olr,-1,rlat65,rlon73,-1,'UMBC');
end
