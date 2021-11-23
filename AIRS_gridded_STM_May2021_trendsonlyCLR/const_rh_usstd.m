addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/

% /home/sergio/MATLABCODE/QUICKTASKS_TELECON/ConstRH/const_rh.m

sarta = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_pge_v6_tunmlt';
fop = 'junk.op.rtp';
frp = 'junk.rp.rtp';

[h,ha,p0,pa] = rtpread('/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.so2.latlon.const_emiss.rtp');
iProf = input('Enter which of reg profiles to subset (1-49) : ');
if length(iProf) == 0
  iProf = 49;
end
[h,p0] = subset_rtp(h,p0,[],[],iProf);

playsN = p0.plevs(1:100)-p0.plevs(2:101);
playsD = log(p0.plevs(1:100)./p0.plevs(2:101));
plays = playsN ./ playsD;

[RH0,RH1km,colwater0] = layeramt2RH(h,p0);
figure(1); semilogy(RH0,plays); set(gca,'ydir','reverse')
xlabel('RH0'); title('Default Humidity')
rtpwrite(fop,h,ha,p0,pa);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
eval(sartaer);
[hx,hax,px0,pax] = rtpread(frp);
figure(2); plot(hx.vchan,rad2bt(hx.vchan,px0.rcalc));
xlabel('wavenumber'); ylabel('BT(K)')
disp('hit <ret> 0'); pause

%% suppose T(z) increases everywhere by 0.1 K
%%%% >>>>>>>>>>> you can put in your favorite dT(z) vector here
pnewT = p0;
dTz = 0.1 * ones(101,1);

plot(p0.ptemp,1:101)   %% can see trop is layers 58-101,strat is above that layers 10-58
dTz = zeros(101,1);
dTz(10:50) = -0.01;  %% strat cooling
dTz(60:p0.nlevs) = +0.01;  %% strat cooling

%%%% >>>>>>>>>>> you can put in your favorite dT(z) vector here
semilogy(dTz,p0.plevs,'o-'); title('this is the dTz profile');
set(gca,'ydir','reverse')
disp('hit <ret>'); pause
pnewT.ptemp = pnewT.ptemp + dTz;
[RHT,RH1km,colwaterT] = layeramt2RH(h,pnewT);
figure(1)
semilogy(RHT,plays); set(gca,'ydir','reverse')
semilogy(RHT-RH0,plays); set(gca,'ydir','reverse')
xlabel('RHT - RH0')
title('Change T by 0.1 K')
rtpwrite(fop,h,ha,pnewT,pa);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
eval(sartaer);
[hx,hax,pxT,pax] = rtpread(frp);
figure(2); plot(hx.vchan,rad2bt(hx.vchan,pxT.rcalc)-rad2bt(hx.vchan,px0.rcalc));
title('Increase T by dT(z)')
xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
disp('hit <ret> T'); pause

%% now see what happens if we CHANGE ie INCREASE wv by the amount (RHT - RH0)./RH0
pnewW = p0;
pnewW.gas_1(1:100) = pnewW.gas_1(1:100) .* (1+(RHT - RH0)./RH0);
[RHW,RH1km,colwaterW] = layeramt2RH(h,pnewW);
figure(1)
semilogy(RHT,plays,'b',RHW,plays,'r'); set(gca,'ydir','reverse')
semilogy(RHT-RH0,plays,'bo-',RHW-RH0,plays,'r'); set(gca,'ydir','reverse')
xlabel('Change in RH')
title('Change T by +0.1K or W by (1+(RHT - RH0)./RH0)')
hl=legend('increase T by dT(z)','INcrease W by (1+(RHT - RH0)./RH0)','location','northeast','fontsize',8);
rtpwrite(fop,h,ha,pnewW,pa);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
eval(sartaer);
[hx,hax,pxW,pax] = rtpread(frp);
figure(2); plot(hx.vchan,rad2bt(hx.vchan,pxW.rcalc)-rad2bt(hx.vchan,px0.rcalc));
title('Change W by (1+(RHT - RH0)./RH0)')
xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
disp('hit <ret> W'); pause

%% finally, what happens if we increase T by 0.1 K and DECREASE WV the above amount?
pnewF = p0;
pnewF.ptemp = pnewF.ptemp + dTz;
pnewF.gas_1(1:100) = pnewF.gas_1(1:100) .* (1 - (RHT - RH0)./RH0);
[RHF,RH1km,colwaterF] = layeramt2RH(h,pnewF);
figure(1)
semilogy(RH0,plays,'b',RHF,plays,'r'); set(gca,'ydir','reverse')
semilogy(RHT-RH0,plays,'b',RH0-RHW,plays,'k',RHF-RH0,plays,'rx-'); set(gca,'ydir','reverse')
hl=legend('increase T by dT(z)','DEcrease W by (1-(RHT - RH0)./RH0)','increase T by dT(z) and DEcrease W by (1-(RHT - RH0)./RH0)','location','northeast','fontsize',8);
xlabel('Change in RH')
title('Change T by +0.1K and W by (1-(RHT - RH0)./RH0)')
grid
rtpwrite(fop,h,ha,pnewF,pa);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
eval(sartaer);
[hx,hax,pxF,pax] = rtpread(frp);
figure(2);
  plot(hx.vchan,rad2bt(hx.vchan,pxT.rcalc)-rad2bt(hx.vchan,px0.rcalc),'b',...
       hx.vchan,rad2bt(hx.vchan,px0.rcalc)-rad2bt(hx.vchan,pxW.rcalc),'k',...  
       hx.vchan,rad2bt(hx.vchan,pxF.rcalc)-rad2bt(hx.vchan,px0.rcalc),'rx-');
title('Change T by dT(z) and W by (1-(RHT - RH0)./RH0)')
xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
hl=legend('increase T by dT(z)','DEcrease W by (1-(RHT - RH0)./RH0)','increase T by dT(z) and DEcrease W by (1-(RHT - RH0)./RH0)','location','northeast','fontsize',8);
grid
disp('hit <ret> F'); pause


