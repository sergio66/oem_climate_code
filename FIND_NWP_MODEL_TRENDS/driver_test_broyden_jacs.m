addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

%% see eqn 17, https://amt.copernicus.org/preprints/amt-2016-181/amt-2016-181.pdf

if ~exist('all')
  load ERA5_atm_data_2002_09_to_2021_08_desc.mat
end

ll = 32;
ind = (1:72) + (ll-1)*72;

plot(all.rlon(ind),all.rlat(ind),'x')

psub.plevs = all.nwp_plevs(:,:,ind);
psub.stemp = all.stemp(:,ind);
psub.ptemp = all.nwp_ptemp(:,:,ind);
psub.gas_1 = all.nwp_gas_1(:,:,ind);
psub.gas_3 = all.nwp_gas_3(:,:,ind);

ii = 1 : 228;
ii1 = 229;
pjunk.plevs = squeeze(psub.plevs(ii,:,1))';
pjunk.stemp = squeeze(psub.stemp(ii,1))';
pjunk.ptemp = squeeze(psub.ptemp(ii,:,1))';
pjunk.gas_1 = squeeze(psub.gas_1(ii,:,1))';
pjunk.gas_3 = squeeze(psub.gas_3(ii,:,1))';
  pjunk.plevs(:,ii1) = nanmean(pjunk.plevs,2);
  pjunk.ptemp(:,ii1) = nanmean(pjunk.ptemp,2);
  pjunk.gas_1(:,ii1) = nanmean(pjunk.gas_1,2);
  pjunk.gas_3(:,ii1) = nanmean(pjunk.gas_3,2);
  pjunk.stemp(ii1) = nanmean(pjunk.stemp);
pjunk.nlevs = ones(size(pjunk.stemp)) * 37;
pjunk.spres = ones(size(pjunk.stemp)) * 1013;
pjunk.nemis = ones(size(pjunk.stemp)) * 2;
pjunk.efreq(1,:) = ones(size(pjunk.stemp)) * 500;
pjunk.efreq(2,:) = ones(size(pjunk.stemp)) * 3500;
pjunk.emis(1,:) = ones(size(pjunk.stemp)) * 0.98;
pjunk.emis(2,:) = ones(size(pjunk.stemp)) * 0.98;
pjunk.rho = (1-pjunk.emis)/pi;
pjunk.scanang = ones(size(pjunk.stemp)) * 0;
pjunk.satzen  = ones(size(pjunk.stemp)) * 0;
pjunk.rlon    = ones(size(pjunk.stemp)) * -177.5;
pjunk.rlat    = ones(size(pjunk.stemp)) * 0;
pjunk.plon = pjunk.rlon;
pjunk.plat = pjunk.rlat;

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';

klayerser = ['!' klayers ' fin=TEST_BROYDEN_JACS/broyden.ip.rtp fout=TEST_BROYDEN_JACS/broyden.op.rtp'];
sartaer = ['!' sarta '     fin=TEST_BROYDEN_JACS/broyden.op.rtp fout=TEST_BROYDEN_JACS/broyden.rp.rtp'];

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat
hjunk = h;
hjunk.ngas = 2;
hjunk.glist = [ 1  3]';
hjunk.gunit = [21 21]';
hjunk.ptype = 0;
hjunk.pfields = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtpwrite('TEST_BROYDEN_JACS/broyden.ip.rtp',hjunk,[],pjunk,[]);
eval(klayerser);
eval(sartaer);

[h,ha,p,pa] = rtpread('TEST_BROYDEN_JACS/broyden.rp.rtp');
prof0 = p;
tcal0 = rad2bt(h.vchan,p.rcalc);
plot(h.vchan,nanmean(tcal0(:,1:228),2),h.vchan,tcal0(:,ii1))
plot(h.vchan,nanmean(tcal0(:,1:228),2)-tcal0(:,ii1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = pjunk;
dT = 1;
px.stemp = pjunk.stemp + dT;

rtpwrite('TEST_BROYDEN_JACS/broyden.ip.rtp',hjunk,[],px,[]);
eval(klayerser);
eval(sartaer);

[h,ha,p,pa] = rtpread('TEST_BROYDEN_JACS/broyden.rp.rtp');
tcalST = rad2bt(h.vchan,p.rcalc);
plot(h.vchan,nanmean(tcalST(:,1:228)-tcal0(:,1:228),2)); title('<ST> jac')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = pjunk;
dT = 1;
px.ptemp = pjunk.ptemp + dT;

rtpwrite('TEST_BROYDEN_JACS/broyden.ip.rtp',hjunk,[],px,[]);
eval(klayerser);
eval(sartaer);

[h,ha,p,pa] = rtpread('TEST_BROYDEN_JACS/broyden.rp.rtp');
tcalTz = rad2bt(h.vchan,p.rcalc);
plot(h.vchan,nanmean(tcalTz(:,1:228)-tcal0(:,1:228),2)); title('<Tz> jac')
jacTz = tcalTz - tcal0; jacTz0 = jacTz(:,ii1); jacTz = jacTz(:,ii);
plot(h.vchan,nanmean(jacTz')-jacTz0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
px = pjunk;
px.gas_1 = pjunk.gas_1 * 1.1;

rtpwrite('TEST_BROYDEN_JACS/broyden.ip.rtp',hjunk,[],px,[]);
eval(klayerser);
eval(sartaer);

[h,ha,p,pa] = rtpread('TEST_BROYDEN_JACS/broyden.rp.rtp');
tcalWVz = rad2bt(h.vchan,p.rcalc);
plot(h.vchan,nanmean(tcalWVz(:,1:228)-tcal0(:,1:228),2)); title('<WVz> jac')
jacWV = tcalWVz - tcal0; jacWV0 = jacWV(:,ii1); jacWV = jacWV(:,ii);
plot(h.vchan,nanmean(jacWV')-jacWV0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = pjunk;
px.gas_3 = pjunk.gas_3 * 1.1;

rtpwrite('TEST_BROYDEN_JACS/broyden.ip.rtp',hjunk,[],px,[]);
eval(klayerser);
eval(sartaer);

[h,ha,p,pa] = rtpread('TEST_BROYDEN_JACS/broyden.rp.rtp');
tcalO3z = rad2bt(h.vchan,p.rcalc);
plot(h.vchan,nanmean(tcalO3z(:,1:228)-tcal0(:,1:228),2)); title('<O3z> jac')
jacO3 = tcalO3z - tcal0; jacO30 = jacO3(:,ii1); jacO3 = jacO3(:,ii);
plot(h.vchan,nanmean(jacO3')-jacO30')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jacST = (tcalST - tcal0)/dT;
  jacST0 = jacST(:,ii1);
  jacST = jacST(:,ii);

jacALL = 
plot(h.vchan,(nanmean(jacTz')-jacTz0')./jacTz0' * 100,...
     h.vchan,(nanmean(jacO3')-jacO30')./jacO30' * 100,...
     h.vchan,(nanmean(jacWV')-jacWV0')./jacWV0' * 100,...
     h.vchan,(nanmean(jacST')-jacST0')./jacST0' * 100)
hl = legend('Tz','O3','WV','stemp','location','best','fontsize',8); ylabel('Percent Error')
axis([640 1640 -10 +10])

plot(h.vchan,(nanmean(jacTz')-jacTz0')./jacTz0' * 100,...
     h.vchan,(nanmean(jacWV')-jacWV0')./jacWV0' * 100,...
     h.vchan,(nanmean(jacST')-jacST0')./jacST0' * 100)
hl = legend('Tz','WV','stemp','location','best','fontsize',8); ylabel('Percent Error')
axis([640 1640 -10 +10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see eqn 17, https://amt.copernicus.org/preprints/amt-2016-181/amt-2016-181.pdfd
x = ones(2645,1) * (pjunk.stemp(ii)-pjunk.stemp(ii1));

Fxnew = tcal0(:,ii);
Fx0   = tcal0(:,ii1) * ones(1,length(ii));
djac = (Fxnew-Fx0)-jacST.*dT;
djac = djac.*dT; 
djac = djac ./ (dT.*dT);
jacSTnew = jacST + djac;

plot(h.vchan,nanmean(jacST'),h.vchan,jacST0)
plot(h.vchan,nanmean(jacST')-jacST0')
plot(h.vchan,nanmean(jacST'),h.vchan,nanmean(jacSTnew'),h.vchan,jacST0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do case 1 by hand v1
jj = 10;
Fxnew = tcal0(:,jj);
Fx0   = tcal0(:,ii1);
dxjj = dx(jj);
djac = (Fxnew-Fx0)-jacST(:,jj)*dxjj;
djac = djac*dxjj; 
djac = djac / (dxjj*dxjj);
jacSTnewjj = jacST(:,jj) + djac;

%% do case 1 by hand v2
jj = 10;
Fxnew = tcal0(:,jj);
Fx0   = tcal0(:,ii1);
dxjj = [prof0.ptemp(1:97,jj)-prof0.ptemp(1:97,ii1)' dx(jj)];
djac = (Fxnew-Fx0)-jacST(:,jj)*dxjj;
djac = djac*dxjj; 
djac = djac / (dxjj*dxjj);
jacSTnewjj = jacST(:,jj) + djac;

plot(h.vchan,jacST(:,jj),h.vchan,jacSTnewjj,h.vchan,jacST0)
