addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE

clear all

sarta = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';

iVers = 1
if iVers == -1
  [h,ha,p,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/Desc/16DayAvgNoS_withanom/latbin30_16day_avg.rp.rtp');
  [h,p] = subset_rtp(h,p,[],[],1);
else
  [h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.rp.rtp');
  [h,p] = subset_rtp(h,pbefore,[],[],20);
end

[hnew,pnew0] = replicate_rtp_headprof(h,p,1,4);  

pnew0 = pnew0; %% unchanged WV,T
pnew0.gas_1(:,2) = pnew0.gas_1(:,2) * 1.1;  %% increasing WV

pnew0.stemp(3) = pnew0.stemp(3) + 1; %% unchanged WV, changed T
pnew0.stemp(4) = pnew0.stemp(4) + 1.0;  %% increasing T
pnew0.gas_1(:,4) = pnew0.gas_1(:,4) * 1.1;  %% increasing WV

rtpwrite('pert_testKC_wv.op.rtp',h,ha,pnew0,pa);
sartaer = ['!' sarta ' fin=pert_testKC_wv.op.rtp fout=pert_testKC_wv.rp.rtp'];
eval(sartaer);
[hx,hax,px,pax] = rtpread('pert_testKC_wv.rp.rtp');
tx = rad2bt(hx.vchan,px.rcalc);

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,px,1:length(px.stemp),2);

diffW1_2 = sum(pnew0.gas_1(1:97,2)-pnew0.gas_1(1:97,1));
diffW3_4 = sum(pnew0.gas_1(1:97,4)-pnew0.gas_1(1:97,3));
wvjac_t0 = (tx(:,2)-tx(:,1))/(diffW1_2);
wvjac_t1 = (tx(:,4)-tx(:,3))/(diffW3_4);

figure(1); 
subplot(211); plot(hx.vchan,wvjac_t1-wvjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,sum(pnew0.gas_1(1:97,1))*(wvjac_t1-wvjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

error('jlkshgslkhjsg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = load('Kcarta_rads/individual_prof_convolved_kcarta_crisHI_crisMED_1.mat');
a2 = load('Kcarta_rads/individual_prof_convolved_kcarta_crisHI_crisMED_2.mat');
a3 = load('Kcarta_rads/individual_prof_convolved_kcarta_crisHI_crisMED_3.mat');
a4 = load('Kcarta_rads/individual_prof_convolved_kcarta_crisHI_crisMED_4.mat');

pkc.rcalc(:,1) = a1.rKc(1:2378);
pkc.rcalc(:,2) = a2.rKc(1:2378);
pkc.rcalc(:,3) = a3.rKc(1:2378);
pkc.rcalc(:,4) = a4.rKc(1:2378);
fkc = a1.fKc(1:2378);
tkc = rad2bt(fkc,pkc.rcalc);

wvjackc_t0 = (tkc(:,2)-tkc(:,1))/(ppmvAVG(2)-ppmvAVG(1));
wvjackc_t1 = (tkc(:,4)-tkc(:,3))/(ppmvAVG(4)-ppmvAVG(3));

figure(2); 
subplot(211); plot(fkc,wvjackc_t1-wvjackc_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(fkc,ppmvAVG(1)*(wvjackc_t1-wvjackc_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

