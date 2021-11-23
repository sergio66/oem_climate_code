addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE

sarta = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';

junk = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs/individual_prof_convolved_kcarta_crisHI_crisMED_20.mat');

[h,ha,pafter,pa] = rtpread('pafteravg1_39.rp.rtp');
[appmvLAY1,appmvAVG1,appmvMAX1,apavgLAY1,atavgLAY1,appmv5001,appmv751] = layers2ppmv(h,pafter,1:length(pafter.stemp),1);
mmw_after = mmwater_rtp(h,pafter);

[h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
[ppmvLAY1,ppmvAVG1,ppmvMAX1,pavgLAY1,tavgLAY1,ppmv5001,ppmv751] = layers2ppmv(h,pbefore,1:length(pbefore.stemp),1);
mmw_before = mmwater_rtp(h,pbefore);

h = rmfield(h,'vchan');
h.vchan = junk.fKc;
h.nchan = 2834;
h.ichan = (1:2834)';
pbefore.robs1 = zeros(2834,length(pbefore.stemp));
pbefore.rcalc = zeros(2834,length(pbefore.stemp));

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,pbefore,1:length(pbefore.stemp),2);

dT = 0:0.05:1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hnew,pnew0] = replicate_rtp_headprof(h,pbefore,20,length(dT));
rtpwrite('pert0.op.rtp',h,ha,pnew0,pa);
sartaer = ['!' sarta ' fin=pert0.op.rtp fout=pert0.rp.rtp'];
eval(sartaer);
[hx,hax,px0,pax] = rtpread('pert0.rp.rtp');
bt0 = rad2bt(h.vchan,px0.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnewT] = replicate_rtp_headprof(h,pbefore,20,length(dT));
for ii = 1 : length(dT)
  pnewT.ptemp(:,ii) = pnewT.ptemp(:,ii) + dT(ii);
end

delT = ones(2834,1) * dT;

rtpwrite('pertT.op.rtp',h,ha,pnewT,pa);
sartaer = ['!' sarta ' fin=pertT.op.rtp fout=pertT.rp.rtp'];
eval(sartaer);
[hx,hax,pT,pax] = rtpread('pertT.rp.rtp');
btT = rad2bt(h.vchan,pT.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnewST] = replicate_rtp_headprof(h,pbefore,20,length(dT));
for ii = 1 : length(dT)
  pnewST.stemp(ii) = pnewST.stemp(ii) + dT(ii);
end

delT = ones(2834,1) * dT;

rtpwrite('pertST.op.rtp',h,ha,pnewST,pa);
sartaer = ['!' sarta ' fin=pertST.op.rtp fout=pertST.rp.rtp'];
eval(sartaer);
[hx,hax,pST,pax] = rtpread('pertST.rp.rtp');
btST = rad2bt(h.vchan,pST.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnew2] = replicate_rtp_headprof(h,pbefore,20,length(dT));

for ii = 1 : length(dT)
  pnew2.ptemp(:,ii) = pnew2.ptemp(:,ii) + dT(ii);
  pnew2.gas_2(:,ii) = pnew2.gas_2(:,ii) * 1.1;;
end

del2 = ones(2834,1) * ones(size(dT)) * log10(1.1);

rtpwrite('pertT2.op.rtp',h,ha,pnew2,pa);
sartaer = ['!' sarta ' fin=pertT2.op.rtp fout=pertT2.rp.rtp'];
eval(sartaer);
[hx,hax,pT2,pax] = rtpread('pertT2.rp.rtp');
btT2 = rad2bt(h.vchan,pT2.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnewTm] = replicate_rtp_headprof(h,pbefore,20,length(dT));
for ii = 1 : length(dT)
  pnewTm.ptemp(:,ii) = pnewTm.ptemp(:,ii) - dT(ii);
end

delT = ones(2834,1) * dT;

rtpwrite('pertTm.op.rtp',h,ha,pnewTm,pa);
sartaer = ['!' sarta ' fin=pertTm.op.rtp fout=pertTm.rp.rtp'];
eval(sartaer);
[hx,hax,pm,pax] = rtpread('pertTm.rp.rtp');
btm = rad2bt(h.vchan,pm.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnewW] = replicate_rtp_headprof(h,pbefore,20,length(dT));
for ii = 1 : length(dT)
  pnewW.gas_1(:,ii) = pnewW.gas_1(:,ii) * (1+dT(ii)/5);
end

delWV = ones(2834,1) * log(1+dT/5);

rtpwrite('pertW.op.rtp',h,ha,pnewW,pa);
sartaer = ['!' sarta ' fin=pertW.op.rtp fout=pertW.rp.rtp'];
eval(sartaer);
[hx,hax,pW,pax] = rtpread('pertW.rp.rtp');
btW = rad2bt(h.vchan,pW.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
g101 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g101_jac.mat');
g102 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g102_jac.mat');

Tjac = (btT-bt0)./delT;
figure(1); plot(h.vchan,(btT-bt0)); title('T(new)-T(0) for different dT')
figure(2); plot(h.vchan,(btT-bt0)./delT); title('T jac for different dT')
%figure(3); plot(h.vchan,(btT2-btT)); title('CO2 jac for different dT')

%Tjac = (btT-bt0)./delT;
%figure(1); plot(h.vchan,Tjac); title('T jac for different dT')
%figure(2); plot(h.vchan,Tjac-Tjac(:,2)); title('dT jac for different dT')

co2jac = (btT2-btT);
figure(3); plot(h.vchan,co2jac); title('CO2 jac for different dT')
figure(4); plot(h.vchan,co2jac-co2jac(:,1)); title('dCO2 jac for different dT')

figure(5); 
Tmjac = (bt0-btm)./delT;
T2jac = (btT-2*bt0+btm)./delT./delT;
figure(5); plot(h.vchan,T2jac); title('d2BT/dT2 jac for different dT')

figure(6); plot(h.vchan,btW-bt0); title('WV jac for different dW')
figure(6); plot(h.vchan,(btW-bt0)./delWV,'b',g1.fout,sum(g1.jout'+g101.jout'+g102.jout'),'g'); title('WV jac for different dW')
%figure(6); plot(h.vchan,btW(:,21)-bt0(:,21)); title('BT change for +20% dW')

figure(7); plot(h.vchan,(btST-bt0)./delT); title('ST jac for different dT')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); plot(h.vchan,  ((btW-bt0)./delWV) ./ ((btST-bt0)./delT) ); title('dST/dWV ')
axis([700 1300 -50 0])

mooWV = ((btW-bt0)./delWV) ./ ((btST-bt0)./delT);
figure(9); plot(h.vchan,mooWV-mooWV(:,2)*ones(1,21))
figure(9); plot(h.vchan,mooWV(:,21)-mooWV(:,2))
title('Sensitivity of jac ST to WV')
axis([700 1300 -2 +2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10); plot(h.vchan,  ((btT-bt0)./delT) ./ ((btST-bt0)./delT) ); title('dST/dT ')
axis([700 1300 0 +10])

mooT = ((btT-bt0)./delT) ./ ((btST-bt0)./delT);
figure(11); plot(h.vchan,mooT-mooT(:,2)*ones(1,21))
figure(11); plot(h.vchan,mooT(:,21)-mooT(:,2))
title('Sensitivity of jac ST to T')
axis([700 1300 -0.05 +0.05])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
