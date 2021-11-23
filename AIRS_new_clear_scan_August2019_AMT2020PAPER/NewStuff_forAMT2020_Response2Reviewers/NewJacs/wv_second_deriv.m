[hnew,pnew0] = replicate_rtp_headprof(hall2,pall2,ii,4);  

%%%%%%%%%%%%%%%%%%%%%%%%%
pnew0 = pnew0; %% unchanged WV,ST
pnew0.gas_1(:,2) = pnew0.gas_1(:,2) * 1.1;  %% increasing WV

pnew0.stemp(3) = pnew0.stemp(3) + 1; %% unchanged WV, changed ST
pnew0.stemp(4) = pnew0.stemp(4) + 1.0;  %% increasing ST
pnew0.gas_1(:,4) = pnew0.gas_1(:,4) * 1.1;  %% increasing WV

rtpwrite('pert_testKC_wv.op.rtp',h,ha,pnew0,pa);
sartaer = ['!' sarta ' fin=pert_testKC_wv.op.rtp fout=pert_testKC_wv.rp.rtp'];
eval(sartaer);
[hx,hax,px,pax] = rtpread('pert_testKC_wv.rp.rtp');
tx = rad2bt(hx.vchan,px.rcalc);

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,px,1:length(px.stemp),1);

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = sum(pnew0.gas_1(1:97,2)-pnew0.gas_1(1:97,1));
diff3_4 = sum(pnew0.gas_1(1:97,4)-pnew0.gas_1(1:97,3));
wvjac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
wvjac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(1); 
subplot(211); plot(hx.vchan,wvjac_t1-wvjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,sum(pnew0.gas_1(1:97,1))*(wvjac_t1-wvjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

savex.water_stemp_true(:,ii) = wvjac_t1-wvjac_t0;
savex.water_stemp_log(:,ii) = sum(pnew0.gas_1(1:97,1))*(wvjac_t1-wvjac_t0);

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = sum(ppmvLAY(1:97,2)-ppmvLAY(1:97,1));
diff3_4 = sum(ppmvLAY(1:97,4)-ppmvLAY(1:97,3));
wvjac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
wvjac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(4); 
subplot(211); plot(hx.vchan,wvjac_t1-wvjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT) ppmv')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,sum(ppmvLAY(1:97,1))*(wvjac_t1-wvjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT) ppmv')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = ppmvAVG(2)-ppmvAVG(1);
diff3_4 = ppmvAVG(4)-ppmvAVG(3);
wvjac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
wvjac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(5); 
subplot(211); plot(hx.vchan,wvjac_t1-wvjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT) <ppmv>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,ppmvAVG(1)*(wvjac_t1-wvjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT) <ppmv>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%
mmw = mmwater_rtp(hnew,pnew0);

diff1_2 = mmw(2)-mmw(1);
diff3_4 = mmw(4)-mmw(3);
wvjac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
wvjac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(6); 
subplot(211); plot(hx.vchan,wvjac_t1-wvjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dWV dT) <colW>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,mmw(1)*(wvjac_t1-wvjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT) <colW>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

savex.water_stemp_true_colW(:,ii) = wvjac_t1-wvjac_t0;
savex.water_stemp_log_colW(:,ii) = mmw(1)*(wvjac_t1-wvjac_t0);
