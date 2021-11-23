addpath /home/sergio/MATLABCODE
sarta = ['/home/chepplew/gitLib/sarta/bin//airs_l1c_2834_may19'];
sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';
sarta = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';

sartaer = ['!' sarta ' fin=junk.op.rtp fout=junk.rp.rtp'];

hnew = h;
pnew = p;
%hnew.vchan = 1 : 2834;

rtpwrite('junk.op.rtp',hnew,ha,pnew,pa);
eval(sartaer);
[hjunk,ha,pjunk,pa] = rtpread('junk.rp.rtp');
tjunk0 = rad2bt(hjunk.vchan,pjunk.rcalc);

for iRTP = 1 : length(p.stemp)
  fprintf(1,'Tjac for iRTP = %3i \n',iRTP)
  [hnew,pnew] = replicate_rtp_headprof(h,p,iRTP,101);
  for jj = 1 : 101
    pnew.ptemp(jj,jj) = pnew.ptemp(jj,jj) + 1;
  end
  rtpwrite('junk.op.rtp',hnew,ha,pnew,pa);
  eval(sartaer);
  [hjunk,ha,pjunk,pa] = rtpread('junk.rp.rtp');
  tjunk = rad2bt(hjunk.vchan,pjunk.rcalc);
  tjacSARTA(iRTP,:,:) = tjunk - tjunk0(:,iRTP)*ones(1,101);
end

figure(1)
plot(hjunk.vchan,sum(squeeze(tjacSARTA(1,:,:))'),'b',hjunk.vchan,sum(squeeze(tjacSARTA(22,:,:))'),'g',hjunk.vchan,sum(squeeze(tjacSARTA(45,:,:))'),'r')
  hl = legend('45 deg','22 deg','0 deg','location','best'); title('SARTA col Temp')

i900 = find(hjunk.vchan >= 900,1);
iInd = (1:97)+2*97;
tjunkS = squeeze(tjacSARTA(:,i900,:)); tjunkK = squeeze(jac(:,i900,iInd));
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-',p.scanang,0.445*1./cos(p.scanang*pi/180),'k'); 
  title('col T at 900 cm-1'); hl = legend('SARTA','KCARTA','0.445*1/cos(scanang)','location','best'); xlabel('scanang')

i725 = find(hjunk.vchan >= 725,1);
iInd = (1:97)+2*97;
tjunkS = squeeze(tjacSARTA(:,i725,:)); tjunkK = squeeze(jac(:,i725,iInd));
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-',p.scanang,0.445*1./cos(p.scanang*pi/180),'k'); 
  title('col T at 725 cm-1'); hl = legend('SARTA','KCARTA','0.445*1/cos(scanang)','location','best'); xlabel('scanang')
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-'); 
  title('col T at 725 cm-1'); hl = legend('SARTA','KCARTA','location','best'); xlabel('scanang')

i700 = find(hjunk.vchan >= 700,1);
iInd = (1:97)+2*97;
tjunkS = squeeze(tjacSARTA(:,i700,:)); tjunkK = squeeze(jac(:,i700,iInd));
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-',p.scanang,0.445*1./cos(p.scanang*pi/180),'k'); 
  title('col T at 700 cm-1'); hl = legend('SARTA','KCARTA','0.445*1/cos(scanang)','location','best'); xlabel('scanang')
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-'); 
  title('col T at 700 cm-1'); hl = legend('SARTA','KCARTA','location','best'); xlabel('scanang')

i665 = find(hjunk.vchan >= 665,1);
iInd = (1:97)+2*97;
tjunkS = squeeze(tjacSARTA(:,i665,:)); tjunkK = squeeze(jac(:,i665,iInd));
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-',p.scanang,0.445*1./cos(p.scanang*pi/180),'k'); 
  title('col T at 665 cm-1'); hl = legend('SARTA','KCARTA','0.445*1/cos(scanang)','location','best'); xlabel('scanang')
plot(p.scanang,sum(tjunkS,2),'bo-',p.scanang,sum(tjunkK,2),'rx-'); 
  title('col T at 665 cm-1'); hl = legend('SARTA','KCARTA','location','best'); xlabel('scanang')


