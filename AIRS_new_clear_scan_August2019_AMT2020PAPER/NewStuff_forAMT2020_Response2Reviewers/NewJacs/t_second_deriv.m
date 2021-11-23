[hnew,pnew0] = replicate_rtp_headprof(hall2,pall2,ii,4);  

%%%%%%%%%%%%%%%%%%%%%%%%%
pnew0 = pnew0; %% unchanged Ta,ST
pnew0.ptemp(:,2) = pnew0.ptemp(:,2) + 1;  %% increasing Ta

pnew0.stemp(3) = pnew0.stemp(3) + 1; %% unchanged Ta, changed ST
pnew0.stemp(4) = pnew0.stemp(4) + 1;  %% increasing ST
pnew0.ptemp(:,4) = pnew0.ptemp(:,4) + 1;  %% increasing Ta

rtpwrite('pert_testKC_wv.op.rtp',h,ha,pnew0,pa);
sartaer = ['!' sarta ' fin=pert_testKC_wv.op.rtp fout=pert_testKC_wv.rp.rtp'];
eval(sartaer);
[hx,hax,px,pax] = rtpread('pert_testKC_wv.rp.rtp');
tx = rad2bt(hx.vchan,px.rcalc);

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,px,1:length(px.stemp),2);

tjac_t0 = (tx(:,2)-tx(:,1))/1;
tjac_t1 = (tx(:,4)-tx(:,3))/1;

figure(2); 
subplot(211); plot(hx.vchan,tjac_t1-tjac_t0)  %% dT = 1 so no need to worry about that
  title('true (d^2BT)/(dTa dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

subplot(212); plot(hx.vchan,(tjac_t1-tjac_t0))  %% dT = 1 so no need to worry about that
  title('WV (d^2BT)/(dWV dT)')
  ax = axis; ax(1) = 640; ax(2) = 1640; axis(ax);

savex.ta_stemp_true(:,ii) = tjac_t1-tjac_t0;
