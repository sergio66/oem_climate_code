function xout = compute_feedbacks_regress_olr_ecRad_calcs(x0,indSST,iLambda_UseGlobalSST_regress,caModelStr)

xout = x0;

if isfield(xout,'feedback')
  xout = rmfield(xout,'feedback');
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

globalSST = nanmean(indSST);

load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
  YY = Y(:)';

allind = 1 : length(YY);
tropics = find(abs(YY) < 30);
midlats = find(abs(YY) >= 30 & abs(YY) < 60);
poles   = find(abs(YY) >= 60);
coslat  = cos(YY*pi/180);

xout.feedback_ecRad.global_coslat_wgt_skt  = sum(indSST .* coslat)/sum(coslat);

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;

% THIS IS UNIFORM PERTURBATION RESPONSE
ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'PLANCK';
xout.feedback_ecRad.savenums(ix,:) = -(xout.planck_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.planck_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.planck_ecRad = xout.feedback_ecRad.planck_ecRad./indSST;
else
  xout.feedback_ecRad.planck_ecRad = xout.feedback_ecRad.planck_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.planck_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.planck_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.planck_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.planck_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.planck_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.planck_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.planck_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.planck_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

% The lapse-rate feedback  is minus the OLR response to the
% difference between the actual temperature response and the uniform
% Planck response, still holding qv fixed.
% SO HERE I COMPUTE ACTUAL RESPONSE TO TA(z) and SKT, AND LATER SUBTRACT OUT THE UNIFORM PLANCK RESPONSE
ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'LAPSE';
xout.feedback_ecRad.savenums(ix,:) = -(xout.lapse_ecRad.clr-xout.planck_ecRad.clr);
xout.feedback_ecRad.lapse_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.lapse_ecRad = xout.feedback_ecRad.lapse_ecRad./indSST;
else
  xout.feedback_ecRad.lapse_ecRad = xout.feedback_ecRad.lapse_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.lapse_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.lapse_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.lapse_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.lapse_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.lapse_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.lapse_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.lapse_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.lapse_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'O3';
xout.feedback_ecRad.savenums(ix,:) = -(xout.o3_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.o3_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.o3_ecRad = xout.feedback_ecRad.o3_ecRad./indSST;
else
  xout.feedback_ecRad.o3_ecRad = xout.feedback_ecRad.o3_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.o3_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.o3_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.o3_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.o3_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.o3_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.o3_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.o3_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.o3_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'WV';
xout.feedback_ecRad.savenums(ix,:) = -(xout.wv_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.wv_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.wv_ecRad = xout.feedback_ecRad.wv_ecRad./indSST;
else
  xout.feedback_ecRad.wv_ecRad = xout.feedback_ecRad.wv_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.wv_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.wv_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.wv_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.wv_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.wv_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.wv_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.wv_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.wv_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'SKT';
xout.feedback_ecRad.savenums(ix,:) = -(xout.skt_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.skt_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.skt_ecRad = xout.feedback_ecRad.skt_ecRad./indSST;
else
  xout.feedback_ecRad.skt_ecRad = xout.feedback_ecRad.skt_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.skt_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.skt_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.skt_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.skt_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.skt_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.skt_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.skt_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.skt_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.skt_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

ix = ix + 1;
xout.feedback_ecRad.savestr{ix}    = 'PTEMP_CO2';
xout.feedback_ecRad.savenums(ix,:) = -(xout.ptemp_co2_ecRad.clr - xout.olr0_ecRad.clr);
xout.feedback_ecRad.ptemp_co2_ecRad = xout.feedback_ecRad.savenums(ix,:);
if iLambda_UseGlobalSST_regress == -1
  xout.feedback_ecRad.ptemp_co2_ecRad = xout.feedback_ecRad.ptemp_co2_ecRad./indSST;
else
  xout.feedback_ecRad.ptemp_co2_ecRad = xout.feedback_ecRad.ptemp_co2_ecRad/globalSST;
end
%%%%%
junk = polyfit(indSST,xout.feedback_ecRad.savenums(ix,:),1);
plot(indSST,xout.feedback_ecRad.savenums(ix,:),'.')
xout.feedback_ecRad.ptemp_co2_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind),1);
  xout.feedback_ecRad.ptemp_co2_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums(ix,:));
  xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums(ix,tropics));
  xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums(ix,midlats));
  xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums(ix,poles));
  xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums(ix,ind));
  xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end
%%%%
xout.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums(ix,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  xout.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums(ix,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(xout,'perts9999')
  iy = 0;

  iy = iy + 1;
  xout.feedback_ecRad.savestr9999{iy}    = 'atm_skt_ghg';
  xout.feedback_ecRad.savenums9999(iy,:) = -(xout.perts9999.atm_skt_ghg_ecRad.clr - xout.olr0_ecRad.clr);
  xout.feedback9999_ecRad.atm_skt_ghg_ecRad = xout.feedback_ecRad.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad = xout.feedback9999_ecRad.atm_skt_ghg_ecRad./indSST;
  else
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad = xout.feedback9999_ecRad.atm_skt_ghg_ecRad/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.feedback_ecRad.savenums9999(iy,:),1);
  plot(indSST,xout.feedback_ecRad.savenums9999(iy,:),'.')
  xout.feedback9999_ecRad.atm_skt_ghg_ecRad_polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind),1);
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums9999(iy,:));
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums9999(iy,tropics));
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums9999(iy,midlats));
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums9999(iy,poles));
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind));
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.feedback9999_ecRad.atm_skt_ghg_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.feedback9999_ecRad.atm_skt_ghg_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  end

  %%%%%%%%%%

  iy = iy + 1;
  xout.feedback_ecRad.savestr9999{iy}    = 'atm_skt';
  xout.feedback_ecRad.savenums9999(iy,:) = -(xout.perts9999.atm_skt_ecRad.clr - xout.olr0_ecRad.clr);
  xout.feedback9999_ecRad.atm_skt_ecRad = xout.feedback_ecRad.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.feedback9999_ecRad.atm_skt_ecRad = xout.feedback9999_ecRad.atm_skt_ecRad./indSST;
  else
    xout.feedback9999_ecRad.atm_skt_ecRad = xout.feedback9999_ecRad.atm_skt_ecRad/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.feedback_ecRad.savenums9999(iy,:),1);
  plot(indSST,xout.feedback_ecRad.savenums9999(iy,:),'.')
  xout.feedback9999_ecRad.atm_skt_ecRad_polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind),1);
    xout.feedback9999_ecRad.atm_skt_ecRad_polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums9999(iy,:));
    xout.feedback9999_ecRad.atm_skt_ecRad_robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums9999(iy,tropics));
    xout.feedback9999_ecRad.atm_skt_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums9999(iy,midlats));
    xout.feedback9999_ecRad.atm_skt_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums9999(iy,poles));
    xout.feedback9999_ecRad.atm_skt_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind));
    xout.feedback9999_ecRad.atm_skt_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.feedback9999_ecRad.atm_skt_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.feedback9999_ecRad.atm_skt_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
  end

  %%%%%%%%%%

  iy = iy + 1;
  xout.feedback_ecRad.savestr9999{iy}    = 'atm_only';
  xout.feedback_ecRad.savenums9999(iy,:) = (xout.perts9999.atm_only_ecRad.clr - xout.olr0_ecRad.clr);
  xout.feedback9999_ecRad.atm_only_ecRad = xout.feedback_ecRad.savenums9999(iy,:);
  if iLambda_UseGlobalSST_regress == -1
    xout.feedback9999_ecRad.atm_only_ecRad = xout.feedback9999_ecRad.atm_only_ecRad./indSST;
  else
    xout.feedback9999_ecRad.atm_only_ecRad = xout.feedback9999_ecRad.atm_only_ecRad/globalSST;
  end
  %%%%%
  junk = polyfit(indSST,xout.feedback_ecRad.savenums9999(iy,:),1);
  plot(indSST,xout.feedback_ecRad.savenums9999(iy,:),'.')
  xout.feedback9999_ecRad.atm_only_ecRad_polyfit = junk(1);
  [nn,nx,ny,nmeanval(iy,:),nstdval(iy,:)] = myhist2d(indSST,xout.feedback_ecRad.savenums9999(iy,:),dx,dx);
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = polyfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind),1);
    xout.feedback9999_ecRad.atm_only_ecRad_polyfit_latbin(jj) = junk(1);
  end
  %%%%%
  [junk junkB] = robustfit(indSST,xout.feedback_ecRad.savenums9999(iy,:));
    xout.feedback9999_ecRad.atm_only_ecRad_robustfit_all = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(tropics),xout.feedback_ecRad.savenums9999(iy,tropics));
    xout.feedback9999_ecRad.atm_only_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(midlats),xout.feedback_ecRad.savenums9999(iy,midlats));
    xout.feedback9999_ecRad.atm_only_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
  [junk junkB] = robustfit(indSST(poles),xout.feedback_ecRad.savenums9999(iy,poles));
    xout.feedback9999_ecRad.atm_only_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    junk = robustfit(indSST(ind),xout.feedback_ecRad.savenums9999(iy,ind));
    xout.feedback9999_ecRad.atm_only_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
  end
  %%%%
  xout.feedback9999_ecRad.atm_only_ecRad_globalSST_weighted = sum(xout.feedback_ecRad.savenums9999(iy,:).*coslat)/sum(coslat)/xout.feedback_ecRad.global_coslat_wgt_skt;
  for jj = 1 : 64
    ind = (jj-1)*72 + ind0;
    xout.feedback9999_ecRad.atm_only_ecRad_globalSST_weighted_latbin(jj) = sum(xout.feedback_ecRad.savenums9999(iy,ind).*coslat(ind))/sum(coslat(ind))/xout.feedback_ecRad.global_coslat_wgt_skt;
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
polyfit_feedbacks = [xout.feedback_ecRad.planck_ecRad_polyfit xout.feedback_ecRad.lapse_ecRad_polyfit ...
                       xout.feedback_ecRad.o3_ecRad_polyfit     xout.feedback_ecRad.wv_ecRad_polyfit ...
                       xout.feedback_ecRad.skt_ecRad_polyfit    xout.feedback_ecRad.ptemp_co2_ecRad_polyfit];
fprintf(1,'polyfit feedbacks ALL   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',polyfit_feedbacks);
fprintf(1,'Jevanjee 2021 GRL : sum(polyfit_feedbacks([1 2 3 4])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(polyfit_feedbacks([3 4 5 6])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([3 4 5 6])))

disp(' ')
globalSST_weighted_feedbacks = [xout.feedback_ecRad.planck_ecRad_globalSST_weighted xout.feedback_ecRad.lapse_ecRad_globalSST_weighted ...
                       xout.feedback_ecRad.o3_ecRad_globalSST_weighted     xout.feedback_ecRad.wv_ecRad_globalSST_weighted ...
                       xout.feedback_ecRad.skt_ecRad_globalSST_weighted    xout.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted];
fprintf(1,'globalSST_weighted feedbacks ALL   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',globalSST_weighted_feedbacks);
fprintf(1,'Jevanjee 2021 GRL : sum(globalSST_weighted_feedbacks([1 2 3 4])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(globalSST_weighted_feedbacks([3 4 5 6])) = %5.2f W/m2/K \n',sum(globalSST_weighted_feedbacks([3 4 5 6])))

disp(' ')
ixix = [1];  %% values
iyiy = [2];  %% unc
robustfit_feedbacks_all = [xout.feedback_ecRad.planck_ecRad_robustfit_all(ixix) xout.feedback_ecRad.lapse_ecRad_robustfit_all(ixix) ...
                           xout.feedback_ecRad.o3_ecRad_robustfit_all(ixix)     xout.feedback_ecRad.wv_ecRad_robustfit_all(ixix) ...
                           xout.feedback_ecRad.skt_ecRad_robustfit_all(ixix)    xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_all(ixix)];
fprintf(1,'robustfit feedbacks ALL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_all);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_all([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_all([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([3 4 5 6])))

disp(' ')
robustfit_feedbacks_tropics = [xout.feedback_ecRad.planck_ecRad_robustfit_tropics(ixix) xout.feedback_ecRad.lapse_ecRad_robustfit_tropics(ixix) ...
                               xout.feedback_ecRad.o3_ecRad_robustfit_tropics(ixix)     xout.feedback_ecRad.wv_ecRad_robustfit_tropics(ixix) ...
                               xout.feedback_ecRad.skt_ecRad_robustfit_tropics(ixix)    xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(ixix)];
fprintf(1,'robustfit feedbacks TRP : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_tropics);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_tropics([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_tropics([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([3 4 5 6])))

disp(' ')
robustfit_feedbacks_midlats = [xout.feedback_ecRad.planck_ecRad_robustfit_midlats(ixix) xout.feedback_ecRad.lapse_ecRad_robustfit_midlats(ixix) ...
                               xout.feedback_ecRad.o3_ecRad_robustfit_midlats(ixix)     xout.feedback_ecRad.wv_ecRad_robustfit_midlats(ixix) ...
                               xout.feedback_ecRad.skt_ecRad_robustfit_midlats(ixix)    xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_midlats(ixix)];
fprintf(1,'robustfit feedbacks MID : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_midlats);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_midlats([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_midlats([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([3 4 5 6])))

disp(' ')
robustfit_feedbacks_poles = [xout.feedback_ecRad.planck_ecRad_robustfit_poles(ixix) xout.feedback_ecRad.lapse_ecRad_robustfit_poles(ixix) ...
                             xout.feedback_ecRad.o3_ecRad_robustfit_poles(ixix)     xout.feedback_ecRad.wv_ecRad_robustfit_poles(ixix) ...
                             xout.feedback_ecRad.skt_ecRad_robustfit_poles(ixix)    xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_poles(ixix)];
fprintf(1,'robustfit feedbacks POL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_poles);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_poles([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_poles([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([3 4 5 6])))

xout.ecRadresults.savestr                       = xout.feedback_ecRad.savestr;
xout.ecRadresults.polyfit_feedbacks             = polyfit_feedbacks;
xout.ecRadresults.globalSST_weighted_feedbacks  = globalSST_weighted_feedbacks;
xout.ecRadresults.robustfit_feedbacks_all       = robustfit_feedbacks_all;
xout.ecRadresults.robustfit_feedbacks_tropics   = robustfit_feedbacks_tropics;
xout.ecRadresults.robustfit_feedbacks_midlats   = robustfit_feedbacks_midlats;
xout.ecRadresults.robustfit_feedbacks_poles     = robustfit_feedbacks_poles;
xout.ecRadresults.dbins                         = dx;
xout.ecRadresults.nmeanval                      = nmeanval;
xout.ecRadresults.nstdval                       = nstdval;

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
plot(rlat,xout.feedback_ecRad.planck_ecRad_polyfit_latbin,...
     rlat,xout.feedback_ecRad.lapse_ecRad_polyfit_latbin,...
     rlat,xout.feedback_ecRad.o3_ecRad_polyfit_latbin,...
     rlat,xout.feedback_ecRad.wv_ecRad_polyfit_latbin,...
     rlat,xout.feedback_ecRad.skt_ecRad_polyfit_latbin,...
     rlat,xout.feedback_ecRad.ptemp_co2_ecRad_polyfit_latbin,...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' POLYFIT'])

figure(5); clf
plot(rlat,xout.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.skt_ecRad_globalSST_weighted_latbin,...
     rlat,xout.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted_latbin,...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' GLOBALSST WEIGHTED'])

figure(6); clf
plot(rlat,xout.feedback_ecRad.planck_ecRad_robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.lapse_ecRad_robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.o3_ecRad_robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.wv_ecRad_robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.skt_ecRad_robustfit_latbin(:,1),...
     rlat,xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_latbin(:,1),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' ROBUSTFIT'])

figure(7); clf
plot(rlat,xout.feedback_ecRad.planck_ecRad_robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.lapse_ecRad_robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.o3_ecRad_robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.wv_ecRad_robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.skt_ecRad_robustfit_latbin(:,2),...
     rlat,xout.feedback_ecRad.ptemp_co2_ecRad_robustfit_latbin(:,2),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title([caModelStr ' ROBUSTFIT UNC'])

figure(8); clf
factor = [0.2 0.2 1 1 0.2 1]'; factor = factor * ones(1,length(dx));
plot(dx,nmeanval(1:4,:),'o-','linewidth',2);         
  plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10);             ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
plot(dx,nmeanval(1:4,:).*factor(1:4,:),'o-','linewidth',2); 
  plotaxis2; hl= legend('0.2*Planck','0.2*Lapse','Ozone','WV','location','south','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
title(caModelStr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
