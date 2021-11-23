[hnew,pnew0] = replicate_rtp_headprof(hall2,pall2,ii,4);  

%%%%%%%%%%%%%%%%%%%%%%%%%
pnew0 = pnew0; %% unchanged WV,ST
pnew0.gas_2(:,2) = pnew0.gas_2(:,2) * 1.1;  %% increasing CO2

pnew0.ptemp(:,3) = pnew0.ptemp(:,3) + 1; %% unchanged CO2, changed ST
pnew0.ptemp(:,4) = pnew0.ptemp(:,4) + 1.0;  %% increasing ST
pnew0.gas_2(:,4) = pnew0.gas_2(:,4) * 1.1;  %% increasing CO2

rtpwrite('pert_testKC_wv.op.rtp',h,ha,pnew0,pa);
sartaer = ['!' sarta ' fin=pert_testKC_wv.op.rtp fout=pert_testKC_wv.rp.rtp'];
eval(sartaer);
[hx,hax,px,pax] = rtpread('pert_testKC_wv.rp.rtp');
tx = rad2bt(hx.vchan,px.rcalc);

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,px,1:length(px.stemp),2);

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = sum(pnew0.gas_2(1:97,2)-pnew0.gas_2(1:97,1));
diff3_4 = sum(pnew0.gas_2(1:97,4)-pnew0.gas_2(1:97,3));
co2jac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
co2jac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(3); 
subplot(211); plot(hx.vchan,co2jac_t1-co2jac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dCO2 dT) molec/cm2')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,sum(pnew0.gas_2(1:97,1))*(co2jac_t1-co2jac_t0))  %% dT = 1 so no need to worry about that
  title('CO2 (d^2BT)/(dCO2 dT) molec/cm2')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

%savex.co2_ta_true(:,ii) = co2jac_t1-co2jac_t0;
%savex.co2_ta_log(:,ii) = sum(pnew0.gas_2(1:97,1))*(co2jac_t1-co2jac_t0);

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = sum(ppmvLAY(1:97,2)-ppmvLAY(1:97,1));
diff3_4 = sum(ppmvLAY(1:97,4)-ppmvLAY(1:97,3));
co2jac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
co2jac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(4); 
subplot(211); plot(hx.vchan,co2jac_t1-co2jac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dCO2 dT) ppmv')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,sum(ppmvLAY(1:97,1))*(co2jac_t1-co2jac_t0))  %% dT = 1 so no need to worry about that
  title('CO2 (d^2BT)/(dCO2 dT) ppmv')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

%savex.co2_ta_true(:,ii) = co2jac_t1-co2jac_t0;
%savex.co2_ta_log(:,ii) = sum(pnew0.gas_2(1:97,1))*(co2jac_t1-co2jac_t0);
%savex.ptemp(:,ii) = pnew0.ptemp(:,1); 
%savex.gas_1(:,ii) = pnew0.gas_1(:,1); 
%savex.gas_2(:,ii) = pnew0.gas_2(:,1); 
%savex.stemp(ii)   = pnew0.stemp(1); 

%%%%%%%%%%%%%%%%%%%%%%%%%
diff1_2 = ppmvAVG(2)-ppmvAVG(1);
diff3_4 = ppmvAVG(4)-ppmvAVG(3);
co2jac_t0 = (tx(:,2)-tx(:,1))/(diff1_2);
co2jac_t1 = (tx(:,4)-tx(:,3))/(diff3_4);

figure(5); 
subplot(211); plot(hx.vchan,co2jac_t1-co2jac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dCO2 dT) <ppmv>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,ppmvAVG(1)*(co2jac_t1-co2jac_t0))  %% dT = 1 so no need to worry about that
  title('CO2 (d^2BT)/(dCO2 dT) <ppmv>')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

savex.co2_ta_true(:,ii) = co2jac_t1-co2jac_t0;
savex.co2_ta_log(:,ii) = sum(pnew0.gas_2(1:97,1))*(co2jac_t1-co2jac_t0);
savex.ptemp(:,ii) = pnew0.ptemp(:,1); 
savex.gas_1(:,ii) = pnew0.gas_1(:,1); 
savex.gas_2(:,ii) = pnew0.gas_2(:,1); 
savex.stemp(ii)   = pnew0.stemp(1); 

