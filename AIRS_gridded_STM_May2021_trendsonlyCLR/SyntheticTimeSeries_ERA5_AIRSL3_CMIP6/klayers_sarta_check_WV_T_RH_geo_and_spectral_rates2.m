rtpwrite(fip,h72,[],p72,[]);
klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ugh'];
sartaer   = ['!' sarta '   fin=' fop ' fout=' frp];

%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing klayers ...')
eval(klayerser);
[h72I,ha72I,p72I,pa72I] = rtpread(fop);
ppmvLAY_2 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),2);
ppmvLAY_4 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),4);
ppmvLAY_6 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),6);

i500 = find(p72I.plevs(:,1) >= 500,1);
p72I.gas_4 = p72I.gas_4 .* (ones(101,1)*(n2oppm_t./ppmvLAY_4(i500,:)));
p72I.gas_6 = p72I.gas_6 .* (ones(101,1)*(ch4ppm_t./ppmvLAY_6(i500,:)));

ppmvLAY_2 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),2);
ppmvLAY_4 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),4);
ppmvLAY_6 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),6);

disp('doing sarta ...')
rtpwrite(fop,h72I,ha72I,p72I,pa72I);
eval(sartaer);
%%%%%%%%%%%%%%%%%%%%%%%%%

% numtimesteps = 144
% [h72,~,p72,~] = rtpread(fip);
[h72x,hax,p72x,pax] = rtpread(frp);
p72x.rh = layeramt2RH(h72x,p72x);
p72x.mmw = mmwater_rtp(h72x,p72x);
if isfield(p72,'verybad')
  p72x.verybad = p72.verybad;
  p72x.lonbin  = p72.lonbin;
end
plot(p72x.rlat,p72x.mmw,'.')
plot(p72x.rlat,p72x.stemp,'.')
plot(p72x.rlat,p72x.ptemp(80,:),'.')
plot(p72x.rlat,p72x.gas_1(80,:),'.')

ppmvLAY_1 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),1);
ppmvLAY_2 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),2);
ppmvLAY_3 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),3);
ppmvLAY_4 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),4);
ppmvLAY_5 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),5);
ppmvLAY_6 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),6);

tcalc = reshape(rad2bt(h72x.vchan,p72x.rcalc),2645,72,numtimesteps);;
tcalcavg = squeeze(nanmean(tcalc,2));

plot(squeeze(nanmean(tcalc(1520,:,:),3)))
plot(1:numtimesteps,tcalcavg(1520,:),1:numtimesteps,squeeze(tcalc(1520,:,:)))
plot(1:numtimesteps,squeeze(tcalc(1520,:,:)),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
plot(1:numtimesteps,nanmean(squeeze(tcalc(1520,:,:))),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%
days = (1:numtimesteps)*30/365;

%%   polyfit(days,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
%%   %Math_tsfit_lin_robust(days*365,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)
%%   
%%   stempjunk = reshape(p72x.stemp,72,numtimesteps);
%%   polyfit(days,nanmean(stempjunk,1),1); ans(1)
%%   %Math_tsfit_lin_robust(days*365,nanmean(stempjunk),4); ans(2)

%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(yy,mm,dd,2002);
%%   disp('dude I just computed dayOFtime')
%% 
%%   polyfit(dayOFtime/365.25,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
%%   %Math_tsfit_lin_robust(dayOFtime,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)
%%   
%%   stempjunk = reshape(p72x.stemp,72,numtimesteps);
%%   polyfit(dayOFtime/365.25,nanmean(stempjunk,1),1); ans(1)
%%   %Math_tsfit_lin_robust(dayOFtime,nanmean(stempjunk),4); ans(2)
  
%%%%%%%%%%%%%%%%%%%%%%%%%

rtpwrite(frp,h72x,hax,p72x,pax);
