function compute_deltaRH = const_rh_4608(deltaST,deltaT,fracWV,rates,mask,daAllSartaSims);

%% compute_deltaRH.orig        = RH0;
%% compute_deltaRH.T_Only      = RHT;
%% compute_deltaRH.Wminus_Only = RHW;
%% compute_deltaRH.TandWminus  = RHF;
%% %compute_deltaRH.umbc        = RHAIRS;
%% compute_deltaRH.final        = RHAIRS;
%%  
%% compute_deltaRH.BT_orig        = rad2bt(hx.vchan,px0.rcalc);        %% always computed
%% compute_deltaRH.BT_T_Only      = rad2bt(hx.vchan,pxT.rcalc);        %% only computed if doAllSartaSims > 0
%% compute_deltaRH.BT_Wminus_Only = rad2bt(hx.vchan,pxWMinus.rcalc);   %% only computed if doAllSartaSims > 0
%% compute_deltaRH.BT_TandWminus  = rad2bt(hx.vchan,pxF.rcalc);        %% only computed if doAllSartaSims > 0
%% %compute_deltaRH.BT_AIRS        = rad2bt(hx.vchan,pxAIRS.rcalc);     %% always computed
%% compute_deltaRH.BT_final       = rad2bt(hx.vchan,pxAIRS.rcalc);     %% always computed
%% 
%% plus compute_deltaRH.surfacelayer and compute_deltaRH.final_surfacetrend

%% daAllSartaSims = +1 = do all sarta calcs, at each test, pretty slow
if nargin == 4
  mask = 1 : length(deltaST);
  daAllSartaSims = -1;
elseif nargin == 5
  daAllSartaSims = -1;
end

[mm,nn] = size(deltaT);
if mm == 100
  deltaT(101,:) = 0;
end
[mm,nn] = size(fracWV);
if mm == 100
  fracWV(101,:) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see derivation in my notebook 42
%% dRH = 0 = dRH/dT dT + dRH/dq dq ==> dq = -(dRH/dT) dT / dRH/dq ==> dq/q = [-(dRH/dT) dT] / [dRH/dq q]
%% now -(dRH/dT) dT ~ -dRH ~ -(RHT-RH0)
%%       dRH/dq q   ~ RH0 since RH ~ q LRT/psat(T)
%% so dq/q = -(RHT-RH0)/RH0
%% ==> qf = q0(1+dq/q0) = q0(1-(RHT-RH0)/RH0)
%% so look, if T goes up, RH goes down rHT < RH0 and so -(RHT-RH0) > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/

%% call it using : const_rh_4608(results(:,6)',deltaT,fracWV,rates,mask); axis([640 1640 -0.1 +0.05])
%% copied from /home/sergio/MATLABCODE/QUICKTASKS_TELECON/ConstRH/const_rh.m

sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
fop = 'junk.op.rtp';
frp = 'junk.rp.rtp';

[h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp'); %% line 23 of see_clust_put_together_jacs_clr.m
[~,kcarta.subjac.ppmv2] = layers2ppmv(h1_4608,p1_4608,1:4608,2);
[~,kcarta.subjac.ppmv4] = layers2ppmv(h1_4608,p1_4608,1:4608,4);
[~,kcarta.subjac.ppmv6] = layers2ppmv(h1_4608,p1_4608,1:4608,6);

h = h1_4608;
p0 = p1_4608;

playsN = p0.plevs(1:100,:)-p0.plevs(2:101,:);
playsD = log(p0.plevs(1:100,:)./p0.plevs(2:101,:));
plays = playsN ./ playsD;
p0.plays = plays;

figure(1); clf
figure(2); clf
figure(3); clf

[RH0,RH1km,colwater0] = layeramt2RH(h,p0);
figure(1); semilogy(nanmean(RH0(:,mask)'),nanmean(plays(:,mask)')); set(gca,'ydir','reverse')
xlabel('RH0'); title('Default Humidity')
rtpwrite(fop,h,[],p0,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
eval(sartaer);
[hx,hax,px0,pax] = rtpread(frp);
figure(2); plot(hx.vchan,nanmean(rad2bt(hx.vchan,px0.rcalc(:,mask))'));
xlabel('wavenumber'); ylabel('BT(K)')
disp('hit <ret> 0'); pause(0.1); %pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% suppose T(z) increases everywhere by 0.1 K
%%%% >>>>>>>>>>> you can put in your favorite dT(z) vector here
pnewT = p0;
dTz = deltaT; 
%%%% >>>>>>>>>>> you can put in your favorite dT(z) vector here
figure(3); semilogy(nanmean(dTz(1:100,mask)'),nanmean(p0.plays(1:100,mask)'),'o-'); title('dTz/dt profile');
set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000])
disp('hit <ret>'); pause(0.1); %pause

pnewT.gas_2 = pnewT.gas_2 * (1 + 2.2/400);
pnewT.gas_4 = pnewT.gas_4 * (1 + 0.8/300);
pnewT.gas_6 = pnewT.gas_6 * (1 + 5.0/1700);
pnewT.stemp = pnewT.stemp + deltaST;
pnewT.ptemp = pnewT.ptemp + dTz;
[RHT,RH1km,colwaterT] = layeramt2RH(h,pnewT);
figure(1)
semilogy(nanmean(RHT(1:100,mask)'),nanmean(plays(1:100,mask)')); set(gca,'ydir','reverse')
semilogy(nanmean(RHT(1:100,mask)'-RH0(1:100,mask)'),nanmean(plays(1:100,mask)'),'o-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000])
xlabel('RHT-RH0')
title('Change T to retrieved deltaT')
if daAllSartaSims > 0
  rtpwrite(fop,h,[],pnewT,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
  eval(sartaer);
  [hx,hax,pxT,pax] = rtpread(frp);
  figure(2); plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxT.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'));
  title('T -> T+dT(z)')
  xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
else
  pxT.rcalc = nan(2645,length(p0.stemp));
end
disp('hit <ret> T'); pause(0.1); %pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% now see what happens if we CHANGE ie INCREASE wv by the amount (RHT-RH0)./RH0 NOT what we want because the explicitly is a negative sign in derivation
%% if T goes up, RHT < RH0 and (RHT-RH0)./RH0 < 0 so this code DEcreases water amount
pnewWPlus = p0;
pnewWPlus.gas_2 = pnewWPlus.gas_2 * (1 + 2.2/400);
pnewWPlus.gas_4 = pnewWPlus.gas_4 * (1 + 0.8/300);
pnewWPlus.gas_6 = pnewWPlus.gas_6 * (1 + 5.0/1700);
pnewWPlus.stemp = pnewWPlus.stemp + deltaST;
pnewWPlus.gas_1(1:100,:) = pnewWPlus.gas_1(1:100,:) .* (1+(RHT-RH0)./RH0);
[RHW,RH1km,colwaterW] = layeramt2RH(h,pnewWPlus);  %% note we are overwriting from the above
figure(1)
semilogy(nanmean(RHT(1:100,mask)'),nanmean(plays(1:100,mask)'),'b',nanmean(RHW(1:100,mask)'),nanmean(plays(1:100,mask)'),'r'); set(gca,'ydir','reverse')
semilogy(nanmean(RHT(1:100,mask)'-RH0(1:100,mask)'),nanmean(plays(1:100,mask)'),'bo-',nanmean(RHW(1:100,mask)'-RH0(1:100,mask)'),nanmean(plays(1:100,mask)'),'r'); set(gca,'ydir','reverse')
xlabel('Change in RH')
title('W->W1=W*(1+(RHT-RH0)./RH0) = W(1-dW)')
hl=legend('T->T1=T+dT(z)','W->W1=W*(1+(RHT-RH0)./RH0)','location','northeast','fontsize',8);
if daAllSartaSims > 0
  rtpwrite(fop,h,[],pnewWPlus,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
  eval(sartaer);
  [hx,hax,pxWPlus,pax] = rtpread(frp);
  figure(2); plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxWPlus.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'));
  title('W->W1=W*(1+(RHT-RH0)./RH0) = W(1-dW)')
  xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
else
  pxWPlus.rcalc = nan(2645,length(p0.stemp));
end
disp('hit <ret> Wplus'); pause(0.1); %pause
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now see what happens if we CHANGE ie DECREASE wv by the amount (RHT-RH0)./RH0, what we want
%% if T goes up, RHT < RH0 and -(RHT-RH0)./RH0 > 0 so this code INcreases water amount
pnewWMinus = p0;
pnewWMinus.gas_2 = pnewWMinus.gas_2 * (1 + 2.2/400);
pnewWMinus.gas_4 = pnewWMinus.gas_4 * (1 + 0.8/300);
pnewWMinus.gas_6 = pnewWMinus.gas_6 * (1 + 5.0/1700);
pnewWMinus.stemp = pnewWMinus.stemp + deltaST;
pnewWMinus.gas_1(1:100,:) = pnewWMinus.gas_1(1:100,:) .* (1-(RHT-RH0)./RH0);
[RHW,RH1km,colwaterW] = layeramt2RH(h,pnewWMinus);  %% note we overwrite this down below
figure(1)
semilogy(nanmean(RHT(1:100,mask)'),nanmean(plays(1:100,mask)'),'b',nanmean(RHW(1:100,mask)'),nanmean(plays(1:100,mask)'),'r'); set(gca,'ydir','reverse')
semilogy(nanmean(RHT(1:100,mask)'-RH0(1:100,mask)'),nanmean(plays(1:100,mask)'),'bo-',nanmean(RHW(1:100,mask)'-RH0(1:100,mask)'),nanmean(plays(1:100,mask)'),'r'); set(gca,'ydir','reverse')
xlabel('Change in RH')
title('W->W1=W*(1-(RHT-RH0)./RH0) = W(1+dW)')
hl=legend('T->T1=T+dT(z)','W->W1=W*(1-(RHT-RH0)./RH0)','location','northeast','fontsize',8);
if daAllSartaSims > 0
  rtpwrite(fop,h,[],pnewWMinus,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
  eval(sartaer);
  [hx,hax,pxWMinus,pax] = rtpread(frp);
  figure(2); plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxWMinus.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'));
  title('W->W1=W*(1-(RHT-RH0)./RH0) = W(1+dW)')
  xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
else
  pxWMinus.rcalc = nan(2645,length(p0.stemp));
end
disp('hit <ret> Wminus'); pause(0.1); %pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, what happens if we increase T by 0.1 K and DE/IN CREASE WV the above amount?
pnewF = p0;
pnewF.gas_2 = pnewF.gas_2 * (1 + 2.2/400);
pnewF.gas_4 = pnewF.gas_4 * (1 + 0.8/300);
pnewF.gas_6 = pnewF.gas_6 * (1 + 5.0/1700);
pnewF.stemp = pnewF.stemp + deltaST;
pnewF.ptemp = pnewF.ptemp + dTz;
pnewF.gas_1(1:100,:) = pnewF.gas_1(1:100,:) .* (1-(RHT-RH0)./RH0);    %% this is correct
%%%%%pnewF.gas_1(1:100,:) = pnewF.gas_1(1:100,:) .* (1+(RHT-RH0)./RH0);    %% this is wrong

[RHF,RH1km,colwaterF] = layeramt2RH(h,pnewF);
figure(1)
semilogy(nanmean(RH0(:,mask)'),nanmean(plays(:,mask)'),'b',nanmean(RHF(:,mask)'),nanmean(plays(:,mask)'),'r'); set(gca,'ydir','reverse')
semilogy(nanmean(RHT(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'b',nanmean(RHW(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'k',...
         nanmean(RHF(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'rx-'); set(gca,'ydir','reverse')
hl=legend('T->T1=T+dT','W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)','T->T1,W->W1','location','northeast','fontsize',8);
xlabel('Change in RH')
title('T->T1=T+dT and W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)')
grid

if daAllSartaSims > 0
  rtpwrite(fop,h,[],pnewF,[]);
  sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
  eval(sartaer);
  [hx,hax,pxF,pax] = rtpread(frp);
else
  pxF.rcalc = nan(2645,length(p0.stemp));
end

figure(2);
  plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxT.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'b',...
       hx.vchan,nanmean(rad2bt(hx.vchan,px0.rcalc(:,mask))'-rad2bt(hx.vchan,pxWMinus.rcalc(:,mask))'),'k',...  
       hx.vchan,nanmean(rad2bt(hx.vchan,pxF.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'rx-');

  plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxT.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'b',...
       hx.vchan,nanmean(rad2bt(hx.vchan,pxWMinus.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'k',...  
       hx.vchan,nanmean(rad2bt(hx.vchan,pxF.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'rx-');
title('T->T1=T+dT and W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)')
xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
hl=legend('T->T1=T+dT','W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)','T->T1,W->W1','location','northeast','fontsize',8);
grid
disp('hit <ret> F'); pause(0.1); %pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% grand finale! put in deltaT and fracWV change to see if we reproduce AIRS trends
pnewAIRS = p0;
pnewAIRS.gas_2 = pnewAIRS.gas_2 * (1 + 2.2/400);
pnewAIRS.gas_4 = pnewAIRS.gas_4 * (1 + 0.8/300);
pnewAIRS.gas_6 = pnewAIRS.gas_6 * (1 + 5.0/1700);
pnewAIRS.stemp = pnewAIRS.stemp + deltaST;
pnewAIRS.ptemp = pnewAIRS.ptemp + dTz;
pnewAIRS.gas_1 = pnewAIRS.gas_1 .* (1 + fracWV);
[RHAIRS,RH1km,colwaterAIRS] = layeramt2RH(h,pnewAIRS);
figure(1)
semilogy(nanmean(RH0(:,mask)'),nanmean(plays(:,mask)'),'b',nanmean(RHAIRS(:,mask)'),nanmean(plays(:,mask)'),'r'); set(gca,'ydir','reverse')
semilogy(nanmean(RHT(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'b',nanmean(RHW(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'k',...
         nanmean(RHF(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'rx-',...
         nanmean(RHAIRS(:,mask)'-RH0(:,mask)'),nanmean(plays(:,mask)'),'cx-'); 
          set(gca,'ydir','reverse')
hl=legend('T->T1=T+dT','W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)','T->T1,W->W1','input \delta T,\delta WV','location','northeast','fontsize',8);
xlabel('Change in RH')
xlabel('dRH/dt')
grid; 
xlim([-0.05 +0.05]); ylim([1 1000])
xlim([-0.10 +0.10]); ylim([1 1000])

rtpwrite(fop,h,[],pnewAIRS,[]);
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];
eval(sartaer);
[hx,hax,pxAIRS,pax] = rtpread(frp);
figure(2);
  plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxT.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'b',...
       hx.vchan,nanmean(rad2bt(hx.vchan,px0.rcalc(:,mask))'-rad2bt(hx.vchan,pxWMinus.rcalc(:,mask))'),'k',...  
       hx.vchan,nanmean(rad2bt(hx.vchan,pxF.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'rx-');

  plot(hx.vchan,nanmean(rad2bt(hx.vchan,pxT.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'b',...
       hx.vchan,nanmean(rad2bt(hx.vchan,pxWMinus.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'k',...  
       hx.vchan,nanmean(rad2bt(hx.vchan,pxF.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'rx-',...
       hx.vchan,nanmean(rad2bt(hx.vchan,pxAIRS.rcalc(:,mask))'-rad2bt(hx.vchan,px0.rcalc(:,mask))'),'cx-',...
       hx.vchan,nanmean(rates'),'g');
title('dT(z),dW=(1-(RHT-RH0)./RH0) and AIRS')
xlabel('wavenumber'); ylabel('\delta BT(K) = New-Orig')
hl=legend('T->T1=T+dT','W->W1=W*(1-(RHT-RH0)./RH0)=W(1+dW)','T->T1,W->W1','input \delta T,\delta WV','AIRS rates','location','northeast','fontsize',8);
grid; axis([640 1640 -0.1 +0.05]);
disp('hit <ret> FA'); pause(0.1); %pause
%% keyboard_nowindow

compute_deltaRH.orig        = RH0;
compute_deltaRH.T_Only      = RHT;
compute_deltaRH.Wminus_Only = RHW;
compute_deltaRH.TandWminus  = RHF;
%compute_deltaRH.umbc        = RHAIRS;
compute_deltaRH.final        = RHAIRS;
 
compute_deltaRH.BT_orig        = rad2bt(hx.vchan,px0.rcalc);        %% always computed
compute_deltaRH.BT_T_Only      = rad2bt(hx.vchan,pxT.rcalc);        %% only computed if doAllSartaSims > 0
compute_deltaRH.BT_Wminus_Only = rad2bt(hx.vchan,pxWMinus.rcalc);   %% only computed if doAllSartaSims > 0
compute_deltaRH.BT_TandWminus  = rad2bt(hx.vchan,pxF.rcalc);        %% only computed if doAllSartaSims > 0
%compute_deltaRH.BT_AIRS        = rad2bt(hx.vchan,pxAIRS.rcalc);     %% always computed
compute_deltaRH.BT_final       = rad2bt(hx.vchan,pxAIRS.rcalc);     %% always computed

compute_deltaRH.surfacelayer   = px0.nlevs-1;
compute_deltaRH.surfacelayer   = px0.nlevs-2;
for ii = 1 : 4608
  %compute_deltaRH.umbc_surfacetrend(ii) = compute_deltaRH.umbc(compute_deltaRH.surfacelayer(ii),ii)-compute_deltaRH.orig(compute_deltaRH.surfacelayer(ii),ii);
  compute_deltaRH.final_surfacetrend(ii) = compute_deltaRH.final(compute_deltaRH.surfacelayer(ii),ii)-compute_deltaRH.orig(compute_deltaRH.surfacelayer(ii),ii);
end
