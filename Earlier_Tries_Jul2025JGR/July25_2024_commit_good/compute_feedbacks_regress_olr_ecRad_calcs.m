function xout = compute_feedbacks_regress_olr_ecRad_calcs(x0,indSST,iLambda_UseGlobalSST_regress,caModelStr)

xout = x0;

if isfield(xout,'feedback')
  xout = rmfield(xout,'feedback');
end
if isfield(xout,'feedback_ecRad')
  xout = rmfield(xout,'feedback_ecRad');
end
if isfield(xout,'ecRadresults')
  xout = rmfield(xout,'ecRadresults');
end

%indSST    = results(:,6)';
ix = 0;

clear nmeanval nstdval

%boo0 = -(umbc_spectral_olr.allperts_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr);
%scatter_coast(p.rlon,p.rlat,100,boo0./indSST); colormap jet
%junk0U = polyfit(indSST,boo0,1);
%[nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo0,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
%xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk0U(1)) ' d(SST) + ' num2str(junk0U(2))])

dx = -1 : 0.025 : +1;
ind0 = 1 : 72;

% OHOH before Sept 28, 2023 I had this backwards so may need to recompute this part all over again and save
do_XX_YY_from_X_Y

coslat  = cos(YY*pi/180);

allind = 1 : length(YY);
tropics = find(abs(YY) < 30);
midlats = find(abs(YY) >= 30 & abs(YY) < 60);
poles   = find(abs(YY) >= 60);

xout.feedback_ecRad.global_coslat_wgt_skt  = sum(indSST .* coslat)/sum(coslat);

globalSST = nanmean(indSST);
globalSST = xout.feedback_ecRad.global_coslat_wgt_skt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Query replace       _ecRad_ -> .
%%    umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit   --> umbc_spectral_olr.feedback_ecRad.planck.polyfit
%%    umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit --> umbc_spectral_olr.feedback_ecRad.planck.robustfit
%% Query replace  _ecRad =     --->   .

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;

% THIS IS UNIFORM PERTURBATION RESPONSE
ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'PLANCK';
xout.feedback_ecRad.savenums(ix,:) = -(xout.planck_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.planck.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.planck.individual = xout.feedback_ecRad.planck.individual./indSST;
else
  xout.feedback_ecRad.planck.individual = xout.feedback_ecRad.planck.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.planck.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.planck.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.planck.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.planck.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.planck.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.planck.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.planck.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.planck.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.planck.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.planck.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.planck.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.planck.globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

% The lapse-rate feedback  is minus the OLR response to the
% difference between the actual temperature response and the uniform
% Planck response, still holding qv fixed.
% SO HERE I COMPUTE ACTUAL RESPONSE TO TA(z) and SKT, AND LATER SUBTRACT OUT THE UNIFORM PLANCK RESPONSE
ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'LAPSE';
xout.feedback_ecRad.savenums(ix,:) = -(xout.lapse_ecRad.clr-xout.planck_ecRad.clr);
xout.feedback_ecRad.lapse.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.lapse.individual = xout.feedback_ecRad.lapse.individual./indSST;
else
  xout.feedback_ecRad.lapse.individual = xout.feedback_ecRad.lapse.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.lapse.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.lapse.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.lapse.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.lapse.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.lapse.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.lapse.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.lapse.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.lapse.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.lapse.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
  wah = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
  plot(indSST,xout.feedback_ecRad.savenums(ix,:),'b.',indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics),'ro'); title(['lapse rate : tropics : ' num2str(wah)]); plotaxis2; pause(0.1);  
xout.feedback_ecRad.lapse.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.lapse.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  wah = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.feedback_ecRad.lapse.globalSST_weighted_latbin(jj) = wah;
  plot(indSST,xout.feedback_ecRad.savenums(ix,:),'b.',indSST(ind),xout.feedback_ecRad.savenums(ix,ind),'ro'); title(['lapse rate : latbin' num2str(jj,'%02d') ' ' num2str(wah)]); plotaxis2; pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'O3';
xout.feedback_ecRad.savenums(ix,:) = -(xout.o3_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.o3.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.o3.individual = xout.feedback_ecRad.o3.individual./indSST;
else
  xout.feedback_ecRad.o3.individual = xout.feedback_ecRad.o3.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.o3.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.o3.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.o3.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.o3.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.o3.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.o3.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.o3.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.o3.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.o3.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.o3.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.o3.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.o3.globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'WV';
xout.feedback_ecRad.savenums(ix,:) = -(xout.wv_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.wv.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.wv.individual = xout.feedback_ecRad.wv.individual./indSST;
else
  xout.feedback_ecRad.wv.individual = xout.feedback_ecRad.wv.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.wv.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.wv.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.wv.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.wv.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.wv.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.wv.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.wv.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.wv.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.wv.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.wv.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.wv.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.wv.globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'SKT';
xout.feedback_ecRad.savenums(ix,:) = -(xout.skt_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.skt.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.skt.individual = xout.feedback_ecRad.skt.individual./indSST;
else
  xout.feedback_ecRad.skt.individual = xout.feedback_ecRad.skt.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.skt.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.skt.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.skt.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.skt.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.skt.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.skt.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.skt.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.skt.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.skt.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.skt.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.skt.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.skt.globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'PTEMP_CO2';
xout.feedback_ecRad.savenums(ix,:) = -(xout.ptemp_co2_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.ptemp_co2.individual = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.ptemp_co2.individual = xout.feedback_ecRad.ptemp_co2.individual./indSST;
else
  xout.feedback_ecRad.ptemp_co2.individual = xout.feedback_ecRad.ptemp_co2.individual/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.ptemp_co2.polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.ptemp_co2.polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.ptemp_co2.robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.ptemp_co2.robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.ptemp_co2.robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.ptemp_co2.robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.ptemp_co2.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.ptemp_co2.globalSST_weighted_all     = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics = sum(xout.feedback_ecRad.savenums(ix,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.ptemp_co2.globalSST_weighted_midlats = sum(xout.feedback_ecRad.savenums(ix,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
xout.feedback_ecRad.ptemp_co2.globalSST_weighted_poles   = sum(xout.feedback_ecRad.savenums(ix,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.ptemp_co2.globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(xout,'perts9999')
%   'yay'
% else
%   'nay'
% end
% 
% keyboard_nowindow
 
if isfield(xout,'perts9999')
  iy = 0;

  iy = iy + 1;
  xout.perts9999.savestr9999{iy}    = 'atm_skt_ghg';
  xout.perts9999.savenums9999(iy,:) = -(xout.perts9999.atm_skt_ghg_ecRad.clr - xout.olr0_ecRad.clr);
  xout.perts9999.atm_skt_ghg.individual = xout.perts9999.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.perts9999.atm_skt_ghg.individual = xout.perts9999.atm_skt_ghg.individual./indSST;
  else
    xout.perts9999.atm_skt_ghg.individual = xout.perts9999.atm_skt_ghg.individual/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.perts9999.savenums9999(iy,:),1);
  plot(indSST,xout.perts9999.savenums9999(iy,:),'.')
  xout.perts9999.atm_skt_ghg.polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.perts9999.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.perts9999.savenums9999(iy,ind),1);
    xout.perts9999.atm_skt_ghg.polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.perts9999.savenums9999(iy,:));
    xout.perts9999.atm_skt_ghg.robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.perts9999.savenums9999(iy,tropics));
    xout.perts9999.atm_skt_ghg.robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.perts9999.savenums9999(iy,midlats));
    xout.perts9999.atm_skt_ghg.robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.perts9999.savenums9999(iy,poles));
    xout.perts9999.atm_skt_ghg.robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.perts9999.savenums9999(iy,ind));
    xout.perts9999.atm_skt_ghg.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.perts9999.atm_skt_ghg.globalSST_weighted_all     = sum(xout.perts9999.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt_ghg.globalSST_weighted_tropics = sum(xout.perts9999.savenums9999(iy,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt_ghg.globalSST_weighted_midlats = sum(xout.perts9999.savenums9999(iy,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt_ghg.globalSST_weighted_poles   = sum(xout.perts9999.savenums9999(iy,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.perts9999.atm_skt_ghg.globalSST_weighted_latbin(jj) = sum(xout.perts9999.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  end

  %%%%%%%%%%

  iy = iy + 1;
  xout.perts9999.savestr9999{iy}    = 'atm_skt';
  xout.perts9999.savenums9999(iy,:) = -(xout.perts9999.atm_skt_ecRad.clr - xout.olr0_ecRad.clr);
  xout.perts9999.atm_skt.individual = xout.perts9999.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.perts9999.atm_skt.individual = xout.perts9999.atm_skt.individual./indSST;
  else
    xout.perts9999.atm_skt.individual = xout.perts9999.atm_skt.individual/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.perts9999.savenums9999(iy,:),1);
  plot(indSST,xout.perts9999.savenums9999(iy,:),'.')
  xout.perts9999.atm_skt.polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.perts9999.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.perts9999.savenums9999(iy,ind),1);
    xout.perts9999.atm_skt.polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.perts9999.savenums9999(iy,:));
    xout.perts9999.atm_skt.robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.perts9999.savenums9999(iy,tropics));
    xout.perts9999.atm_skt.robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.perts9999.savenums9999(iy,midlats));
    xout.perts9999.atm_skt.robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.perts9999.savenums9999(iy,poles));
    xout.perts9999.atm_skt.robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.perts9999.savenums9999(iy,ind));
    xout.perts9999.atm_skt.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.perts9999.atm_skt.globalSST_weighted_all     = sum(xout.perts9999.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt.globalSST_weighted_tropics = sum(xout.perts9999.savenums9999(iy,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt.globalSST_weighted_midlats = sum(xout.perts9999.savenums9999(iy,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_skt.globalSST_weighted_poles   = sum(xout.perts9999.savenums9999(iy,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.perts9999.atm_skt.globalSST_weighted_latbin(jj) = sum(xout.perts9999.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  end

  %%%%%%%%%%

  iy = iy + 1;
  xout.perts9999.savestr9999{iy}    = 'atm_only';
  xout.perts9999.savenums9999(iy,:) = -(xout.perts9999.atm_only_ecRad.clr - xout.olr0_ecRad.clr);
  xout.perts9999.atm_only.individual = xout.perts9999.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.perts9999.atm_only.individual = xout.perts9999.atm_only.individual./indSST;
  else
    xout.perts9999.atm_only.individual = xout.perts9999.atm_only.individual/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.perts9999.savenums9999(iy,:),1);
  plot(indSST,xout.perts9999.savenums9999(iy,:),'.')
  xout.perts9999.atm_only.polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.perts9999.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.perts9999.savenums9999(iy,ind),1);
    xout.perts9999.atm_only.polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.perts9999.savenums9999(iy,:));
    xout.perts9999.atm_only.robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.perts9999.savenums9999(iy,tropics));
    xout.perts9999.atm_only.robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.perts9999.savenums9999(iy,midlats));
    xout.perts9999.atm_only.robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.perts9999.savenums9999(iy,poles));
    xout.perts9999.atm_only.robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.perts9999.savenums9999(iy,ind));
    xout.perts9999.atm_only.robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.perts9999.atm_only.globalSST_weighted_all     = sum(xout.perts9999.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_only.globalSST_weighted_tropics = sum(xout.perts9999.savenums9999(iy,tropics).*coslat(tropics))/sum(coslat(tropics))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_only.globalSST_weighted_midlats = sum(xout.perts9999.savenums9999(iy,midlats).*coslat(midlats))/sum(coslat(midlats))/xout.feedback_ecRad.global_coslat_wgt_skt;
  xout.perts9999.atm_only.globalSST_weighted_poles   = sum(xout.perts9999.savenums9999(iy,poles).*coslat(poles))/sum(coslat(poles))/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.perts9999.atm_only.globalSST_weighted_latbin(jj) = sum(xout.perts9999.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
factor = [0.2 0.2 1 1 0.2 1]'; factor = factor * ones(1,length(dx));
plot(dx,nmeanval,'o-','linewidth',2);         
  plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10);             ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
plot(dx,nmeanval.*factor,'o-','linewidth',2); 
  axis([-1 +1 -1 +1]*0.25); 
  plotaxis2; hl= legend('0.2*planck','0.2*lapse','o3','wv','0.2*skt','t/co2','location','south','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');

disp(' ')
%% used to be called "regressed_feedbacks"
polyfit_feedbacks = [xout.feedback_ecRad.planck.polyfit xout.feedback_ecRad.lapse.polyfit ...
                       xout.feedback_ecRad.o3.polyfit     xout.feedback_ecRad.wv.polyfit ...
                       xout.feedback_ecRad.skt.polyfit    xout.feedback_ecRad.ptemp_co2.polyfit];
fprintf(1,'polyfit feedbacks ALL   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',polyfit_feedbacks);
fprintf(1,'Jevanjee 2021 GRL : sum(polyfit_feedbacks([1 2 3 4])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(polyfit_feedbacks([3 4 5 6])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([3 4 5 6])))

disp('<<<<       >>>>')
disp(' ')
globalSST_weighted_feedbacks_all = [xout.feedback_ecRad.planck.globalSST_weighted_all xout.feedback_ecRad.lapse.globalSST_weighted_all ...
                                   xout.feedback_ecRad.o3.globalSST_weighted_all     xout.feedback_ecRad.wv.globalSST_weighted_all ...
                                   xout.feedback_ecRad.skt.globalSST_weighted_all    xout.feedback_ecRad.ptemp_co2.globalSST_weighted_all];
fprintf(1,'globalSST_weighted feedbacks ALL   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',globalSST_weighted_feedbacks_all);
fprintf(1,'Jevanjee 2021 GRL : sum(globalSST_weighted_feedbacks_all([1 2 3 4])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_all([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(globalSST_weighted_feedbacks_all([3 4 5 6])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_all([3 4 5 6])))

disp(' ')
globalSST_weighted_feedbacks_tropics = [xout.feedback_ecRad.planck.globalSST_weighted_tropics xout.feedback_ecRad.lapse.globalSST_weighted_tropics ...
                                   xout.feedback_ecRad.o3.globalSST_weighted_tropics     xout.feedback_ecRad.wv.globalSST_weighted_tropics ...
                                   xout.feedback_ecRad.skt.globalSST_weighted_tropics    xout.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics];
fprintf(1,'globalSST_weighted feedbacks TROPICS   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',globalSST_weighted_feedbacks_tropics);
fprintf(1,'Jevanjee 2021 GRL : sum(globalSST_weighted_feedbacks_tropics([1 2 3 4])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_tropics([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(globalSST_weighted_feedbacks_tropics([3 4 5 6])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_tropics([3 4 5 6])))

disp(' ')
globalSST_weighted_feedbacks_midlats = [xout.feedback_ecRad.planck.globalSST_weighted_midlats xout.feedback_ecRad.lapse.globalSST_weighted_midlats ...
                                   xout.feedback_ecRad.o3.globalSST_weighted_midlats     xout.feedback_ecRad.wv.globalSST_weighted_midlats ...
                                   xout.feedback_ecRad.skt.globalSST_weighted_midlats    xout.feedback_ecRad.ptemp_co2.globalSST_weighted_midlats];
fprintf(1,'globalSST_weighted feedbacks MIDLATS   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',globalSST_weighted_feedbacks_midlats);
fprintf(1,'Jevanjee 2021 GRL : sum(globalSST_weighted_feedbacks_midlats([1 2 3 4])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_midlats([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(globalSST_weighted_feedbacks_midlats([3 4 5 6])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_midlats([3 4 5 6])))

disp(' ')
globalSST_weighted_feedbacks_poles = [xout.feedback_ecRad.planck.globalSST_weighted_poles xout.feedback_ecRad.lapse.globalSST_weighted_poles ...
                                   xout.feedback_ecRad.o3.globalSST_weighted_poles     xout.feedback_ecRad.wv.globalSST_weighted_poles ...
                                   xout.feedback_ecRad.skt.globalSST_weighted_poles    xout.feedback_ecRad.ptemp_co2.globalSST_weighted_poles];
fprintf(1,'globalSST_weighted feedbacks POLES   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',globalSST_weighted_feedbacks_poles);
fprintf(1,'Jevanjee 2021 GRL : sum(globalSST_weighted_feedbacks_poles([1 2 3 4])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_poles([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(globalSST_weighted_feedbacks_poles([3 4 5 6])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks_poles([3 4 5 6])))

disp('<<<<       >>>>')
disp(' ')
ixix = [1];  %% values
iyiy = [2];  %% unc
robustfit_feedbacks_all = [xout.feedback_ecRad.planck.robustfit_all(ixix) xout.feedback_ecRad.lapse.robustfit_all(ixix) ...
                           xout.feedback_ecRad.o3.robustfit_all(ixix)     xout.feedback_ecRad.wv.robustfit_all(ixix) ...
                           xout.feedback_ecRad.skt.robustfit_all(ixix)    xout.feedback_ecRad.ptemp_co2.robustfit_all(ixix)];
fprintf(1,'robustfit feedbacks ALL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_all);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_all([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_all([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([3 4 5 6])))

disp(' ')
robustfit_feedbacks_tropics = [xout.feedback_ecRad.planck.robustfit_tropics(ixix) xout.feedback_ecRad.lapse.robustfit_tropics(ixix) ...
                               xout.feedback_ecRad.o3.robustfit_tropics(ixix)     xout.feedback_ecRad.wv.robustfit_tropics(ixix) ...
                               xout.feedback_ecRad.skt.robustfit_tropics(ixix)    xout.feedback_ecRad.ptemp_co2.robustfit_tropics(ixix)];
fprintf(1,'robustfit feedbacks TRP : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_tropics);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_tropics([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_tropics([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([3 4 5 6])))

disp(' ')
robustfit_feedbacks_midlats = [xout.feedback_ecRad.planck.robustfit_midlats(ixix) xout.feedback_ecRad.lapse.robustfit_midlats(ixix) ...
                               xout.feedback_ecRad.o3.robustfit_midlats(ixix)     xout.feedback_ecRad.wv.robustfit_midlats(ixix) ...
                               xout.feedback_ecRad.skt.robustfit_midlats(ixix)    xout.feedback_ecRad.ptemp_co2.robustfit_midlats(ixix)];
fprintf(1,'robustfit feedbacks MID : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_midlats);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_midlats([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_midlats([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([3 4 5 6])))

disp(' ')
robustfit_feedbacks_poles = [xout.feedback_ecRad.planck.robustfit_poles(ixix) xout.feedback_ecRad.lapse.robustfit_poles(ixix) ...
                             xout.feedback_ecRad.o3.robustfit_poles(ixix)     xout.feedback_ecRad.wv.robustfit_poles(ixix) ...
                             xout.feedback_ecRad.skt.robustfit_poles(ixix)    xout.feedback_ecRad.ptemp_co2.robustfit_poles(ixix)];
fprintf(1,'robustfit feedbacks POL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_poles);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_poles([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_poles([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([3 4 5 6])))

xout.ecRadresults.savestr                              = xout.feedback_ecRad.savestr;
xout.ecRadresults.polyfit_feedbacks                    = polyfit_feedbacks;
xout.ecRadresults.globalSST_weighted_feedbacks_all     = globalSST_weighted_feedbacks_all;
xout.ecRadresults.globalSST_weighted_feedbacks_tropics = globalSST_weighted_feedbacks_tropics;
xout.ecRadresults.globalSST_weighted_feedbacks_midlats = globalSST_weighted_feedbacks_midlats;
xout.ecRadresults.globalSST_weighted_feedbacks_poles   = globalSST_weighted_feedbacks_poles;
xout.ecRadresults.robustfit_feedbacks_all              = robustfit_feedbacks_all;
xout.ecRadresults.robustfit_feedbacks_tropics          = robustfit_feedbacks_tropics;
xout.ecRadresults.robustfit_feedbacks_midlats          = robustfit_feedbacks_midlats;
xout.ecRadresults.robustfit_feedbacks_poles            = robustfit_feedbacks_poles;
xout.ecRadresults.dbins                                = dx;
xout.ecRadresults.nmeanval                             = nmeanval;
xout.ecRadresults.nstdval                              = nstdval;

%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 3
  caModelStr = 'XYZ';
end

figure(3); clf; 
plot(indSST,xout.feedback_ecRad.savenums,'.'); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); 
ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
title(caModelStr);

figure(4); clf
plot(rlat,xout.feedback_ecRad.planck.polyfit_latbin,...
     rlat,xout.feedback_ecRad.lapse.polyfit_latbin,...
     rlat,xout.feedback_ecRad.o3.polyfit_latbin,...
     rlat,xout.feedback_ecRad.wv.polyfit_latbin,...
     rlat,xout.feedback_ecRad.skt.polyfit_latbin,...
     rlat,xout.feedback_ecRad.ptemp_co2.polyfit_latbin,...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' POLYFIT'])
axis([-90 +90 -5 +5])

figure(5); clf
plot(rlat,xout.feedback_ecRad.planck.globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.lapse.globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.o3.globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.skt.globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.ptemp_co2.globalSST_weighted_latbin,...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' GLOBALSST WEIGHTED'])
axis([-90 +90 -5 +5])

figure(6); clf
plot(rlat,xout.feedback_ecRad.planck.robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.lapse.robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.o3.robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.wv.robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.skt.robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.ptemp_co2.robustfit_latbin(:,1),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' ROBUSTFIT'])
axis([-90 +90 -5 +5])

figure(7); clf
plot(rlat,xout.feedback_ecRad.planck.robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.lapse.robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.o3.robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.wv.robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.skt.robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.ptemp_co2.robustfit_latbin(:,2),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' ROBUSTFIT UNC'])
xlim([-90 +90])

figure(8); clf
factor = [0.2 0.2 1 1 0.2 1]'; factor = factor * ones(1,length(dx));
plot(dx,nmeanval(1:4,:),'o-','linewidth',2);         
  plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10);             ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
plot(dx,nmeanval(1:4,:).*factor(1:4,:),'o-','linewidth',2); 
  plotaxis2; hl= legend('0.2*Planck','0.2*Lapse','Ozone','WV','location','south','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
title(caModelStr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
