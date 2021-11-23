%{
sarta_clear = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';

[hx,hax,px,pax] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/regr49_1013_400ppm.op.rtp');

[h,ha,p,pa] = rtpread('pbeforeavg.op.rtp');
%h.ngas = 6; h.gunit = h.gunit(1:6); h.glist = h.glist(1:6);
p.gas_9 = px.gas_9(:,49) * ones(1,40);
p.gas_12 = px.gas_12(:,49) * ones(1,40);
p.plat = p.rlat; p.plon = p.rlon;
[h,p] = subset_rtp(h,p,[],[],1:39); rtpwrite('pbeforeavg1_39.op.rtp',h,ha,p,pa);
/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pbeforeavg1_39.op.rtp fout=pbeforeavg1_39.rp.rtp

[h,ha,p,pa] = rtpread('pafteravg.op.rtp');
%h.ngas = 6; h.gunit = h.gunit(1:6); h.glist = h.glist(1:6);
p.gas_9 = px.gas_9(:,49) * ones(1,40);
p.gas_12 = px.gas_12(:,49) * ones(1,40);
p.plat = p.rlat; p.plon = p.rlon;
[h,p] = subset_rtp(h,p,[],[],1:39); rtpwrite('pafteravg1_39.op.rtp',h,ha,p,pa);
/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pafteravg1_39.op.rtp fout=pafteravg1_39.rp.rtp
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');

%{
pert = pbefore; pert.stemp = pert.stemp + 1;
rtpwrite('pbeforeavg1_39_pertstempP.op.rtp',h,ha,pert,pa);
!/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pbeforeavg1_39_pertstempP.op.rtp fout=pbeforeavg1_39_pertstempP.rp.rtp
%}
[h,ha,pbeforepertP,pa] = rtpread('pbeforeavg1_39_pertstempP.rp.rtp');

%{
pert = pbefore; pert.stemp = pert.stemp - 1;
rtpwrite('pbeforeavg1_39_pertstempM.op.rtp',h,ha,pert,pa);
!/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pbeforeavg1_39_pertstempM.op.rtp fout=pbeforeavg1_39_pertstempM.rp.rtp
%}
[h,ha,pbeforepertM,pa] = rtpread('pbeforeavg1_39_pertstempM.rp.rtp');

%{
pert = pbefore; pert.ptemp = pert.ptemp + 1;
rtpwrite('pbeforeavg1_39_pertptempP.op.rtp',h,ha,pert,pa);
!/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pbeforeavg1_39_pertptempP.op.rtp fout=pbeforeavg1_39_pertptempP.rp.rtp
%}
[h,ha,pbeforepertPTz,pa] = rtpread('pbeforeavg1_39_pertptempP.rp.rtp');

%{
pert = pbefore; pert.ptemp = pert.ptemp - 1;
rtpwrite('pbeforeavg1_39_pertptempM.op.rtp',h,ha,pert,pa);
!/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003 fin=pbeforeavg1_39_pertptempM.op.rtp fout=pbeforeavg1_39_pertptempM.rp.rtp
%}
[h,ha,pbeforepertMTz,pa] = rtpread('pbeforeavg1_39_pertptempM.rp.rtp');


[h,ha,pafter,pa] = rtpread('pafteravg1_39.rp.rtp');
f2378 = h.vchan;

tbefore        = rad2bt(h.vchan,pbefore.rcalc);
tbeforepertP   = rad2bt(h.vchan,pbeforepertP.rcalc);
tbeforepertM   = rad2bt(h.vchan,pbeforepertM.rcalc);
tbeforepertPTz = rad2bt(h.vchan,pbeforepertPTz.rcalc);
tbeforepertMTz = rad2bt(h.vchan,pbeforepertMTz.rcalc);
tafter         = rad2bt(h.vchan,pafter.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /asl/matlib/h4tools
[h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
[h,ha,pafter,pa]  = rtpread('pafteravg1_39.rp.rtp');

iiBin = 30;
iiBin = 20;

w0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);    w1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);
a0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']); a1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']);

figure(1); plot(f2378,tbeforepertP(:,iiBin)-tbefore(:,iiBin),'b',a0.f,a0.jacST,'r'); title('Stemp Jacobian')
  hl = legend('SARTA finite diff','kCARTA analytic','location','best');
figure(2); plot(f2378,tbeforepertP(:,iiBin)-2*tbefore(:,iiBin)+tbeforepertM(:,iiBin),a0.f,a1.jacST-a0.jacST); title('Stemp Hessian')
  hl = legend('SARTA finite diff','kCARTA finite diff','location','best');

figure(4); plot(f2378,tbeforepertPTz(:,iiBin)-tbefore(:,iiBin),'b',a0.f,sum(a0.jacT'),'r'); title('Ptemp Jacobian')
  hl = legend('SARTA finite diff','kCARTA analytic','location','best');
figure(5); plot(f2378,tbeforepertPTz(:,iiBin)-2*tbefore(:,iiBin)+tbeforepertMTz(:,iiBin),a0.f,sum(a1.jacT'-a0.jacT')); title('Ptemp Hessian')
  hl = legend('SARTA finite diff','kCARTA finite diff','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
https://en.wikipedia.org/wiki/Broyden%27s_method

J(n) = J(n-1) +  dy(n) - J(n-1) dx(n)
                 -------------------- dx(n)'
                       |dx(n)|^2

units dy = BT
      J = d(BT)/dx
thus  BT - Bt/dx dx  dx      = dBT/dx = jacobian
      -------------
          dx^2
%}

dy = tbeforepertP(:,iiBin)-tbefore(:,iiBin); dx = 1.0;
dJ = (dy - a0.jacST(1:2378)*dx)/(dx*dx) * dx;
figure(3); plot(f2378,tbeforepertP(:,iiBin)-2*tbefore(:,iiBin)+tbeforepertM(:,iiBin),'b',a0.f,a1.jacST-a0.jacST,'r',f2378,dJ,'k'); title('Stemp Hessian + estimate')
  hl = legend('SARTA finite diff','kCARTA finite diff','Broyden','location','best');

dy = tbeforepertPTz(:,iiBin)-tbefore(:,iiBin); dx = ones(97,1);
dJ = (dy - a0.jacT(1:2378,:)*dx)/(dx'*dx) * dx';
figure(6); plot(f2378,tbeforepertPTz(:,iiBin)-2*tbefore(:,iiBin)+tbeforepertMTz(:,iiBin),'b',a0.f,sum(a1.jacT'-a0.jacT'),'r',f2378,sum(dJ'),'k'); title('Ptemp Hessian + estimate')
  hl = legend('SARTA finite diff','kCARTA finite diff','Broyden','location','best');
