function xout = compute_feedbacks_regress_olr_ecRad_calcs(x0,indSST,iLambda_UseGlobalSST_regress)

xout = x0;

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
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;

% SO THIS IS UNIFORM PERTURBATION RESPONSE
xout.feedback.planck_ecRad = xout.planck_ecRad.clr-xout.olr0_ecRad.clr; 
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.planck_ecRad = -xout.feedback.planck_ecRad./indSST;
else
  xout.feedback.planck_ecRad = -xout.feedback.planck_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.planck_ecRad.clr-xout.olr0_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.planck_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.planck_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.planck_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.planck_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.planck_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.planck_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.planck_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

% The lapse-rate feedback λlapse is minus the OLR response to the
% difference between the actual temperature response and the uniform
% Planck response, still holding qv fixed.
% SO HERE I COMPUTE ACTUAL RESPONSE TO TA(z) and SKT, AND LATER SUBTRACT OUT THE UNIFORM PLANCK RESPONSE
xout.feedback.lapse_ecRad = xout.lapse_ecRad.clr-xout.planck_ecRad.clr;
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.lapse_ecRad = -xout.feedback.lapse_ecRad./indSST;
else
  xout.feedback.lapse_ecRad = -xout.feedback.lapse_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.lapse_ecRad.clr-xout.planck_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.lapse_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.lapse_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.lapse_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.lapse_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.lapse_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.lapse_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.lapse_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

xout.feedback.o3_ecRad = xout.o3_ecRad.clr-xout.olr0_ecRad.clr; 
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.o3_ecRad = -xout.feedback.o3_ecRad./indSST;
else
  xout.feedback.o3_ecRad = -xout.feedback.o3_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.o3_ecRad.clr-xout.olr0_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.o3_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.o3_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.o3_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.o3_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.o3_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.o3_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.o3_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

xout.feedback.wv_ecRad = xout.wv_ecRad.clr-xout.olr0_ecRad.clr;
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.wv_ecRad = -xout.feedback.wv_ecRad./indSST;
else
  xout.feedback.wv_ecRad = -xout.feedback.wv_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.wv_ecRad.clr-xout.olr0_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.wv_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.wv_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.wv_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.wv_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.wv_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.wv_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.wv_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

xout.feedback.skt_ecRad = xout.skt_ecRad.clr-xout.olr0_ecRad.clr;
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.skt_ecRad = -xout.feedback.skt_ecRad./indSST;
else
  xout.feedback.skt_ecRad = -xout.feedback.skt_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.skt_ecRad.clr-xout.olr0_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.skt_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.skt_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.skt_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.skt_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.skt_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.skt_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.skt_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

xout.feedback.ptemp_co2_ecRad = xout.ptemp_co2_ecRad.clr-xout.olr0_ecRad.clr; 
if iLambda_UseGlobalSST_regress == -1
  xout.feedback.ptemp_co2_ecRad = -xout.feedback.ptemp_co2_ecRad./indSST;
else
  xout.feedback.ptemp_co2_ecRad = -xout.feedback.ptemp_co2_ecRad/globalSST;
end
ix = ix + 1;
savenums(ix,:) = (xout.ptemp_co2_ecRad.clr-xout.olr0_ecRad.clr);
%%%%%
junk = polyfit(indSST,-savenums(ix,:),1);
plot(indSST,-savenums(ix,:),'.')
xout.feedback.ptemp_co2_ecRad_polyfit = junk(1);
[nn,nx,ny,nmeanval(ix,:),nstdval(ix,:)] = myhist2d(indSST,-savenums(ix,:),dx,dx);
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = polyfit(indSST(ind),-savenums(ix,ind),1);
  xout.feedback.ptemp_co2_ecRad_polyfit_latbin(jj) = junk(1);
end
%%%%%
[junk junkB] = robustfit(indSST,-savenums(ix,:));
  xout.feedback.ptemp_co2_ecRad_robustfit_all = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(tropics),-savenums(ix,tropics));
  xout.feedback.ptemp_co2_ecRad_robustfit_tropics = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(midlats),-savenums(ix,midlats));
  xout.feedback.ptemp_co2_ecRad_robustfit_midlats = [junk(2) junkB.se(2)];
[junk junkB] = robustfit(indSST(poles),-savenums(ix,poles));
  xout.feedback.ptemp_co2_ecRad_robustfit_poles = [junk(2) junkB.se(2)];
for jj = 1 : 64
  ind = (jj-1)*72 + ind0;
  junk = robustfit(indSST(ind),-savenums(ix,ind));
  xout.feedback.ptemp_co2_ecRad_robustfit_latbin(jj,:) = [junk(2) junkB.se(2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf
factor = [0.2 0.2 1 1 0.2 1]'; factor = factor * ones(1,length(dx));
plot(dx,nmeanval,'o-','linewidth',2);         
  plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10);             ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
plot(dx,nmeanval.*factor,'o-','linewidth',2); 
  axis([-1 +1 -1 +1]*0.25); 
  plotaxis2; hl= legend('0.2*planck','0.2*lapse','o3','wv','0.2*skt','t/co2','location','south','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');

disp(' ')
%% used to be regressed_feedbacks
polyfit_feedbacks = [xout.feedback.planck_ecRad_polyfit xout.feedback.lapse_ecRad_polyfit ...
                       xout.feedback.o3_ecRad_polyfit     xout.feedback.wv_ecRad_polyfit ...
                       xout.feedback.skt_ecRad_polyfit    xout.feedback.ptemp_co2_ecRad_polyfit];
fprintf(1,'polyfit feedbacks ALL   : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',polyfit_feedbacks);
fprintf(1,'Jevanjee 2021 GRL : sum(polyfit_feedbacks([1 2 3 4])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(polyfit_feedbacks([3 4 5 6])) = %5.2f W/m2/K \n',sum(polyfit_feedbacks([3 4 5 6])))

disp(' ')
ixix = [1];  %% values
iyiy = [2];  %% unc
robustfit_feedbacks_all = [xout.feedback.planck_ecRad_robustfit_all(ixix) xout.feedback.lapse_ecRad_robustfit_all(ixix) ...
                           xout.feedback.o3_ecRad_robustfit_all(ixix)     xout.feedback.wv_ecRad_robustfit_all(ixix) ...
                           xout.feedback.skt_ecRad_robustfit_all(ixix)    xout.feedback.ptemp_co2_ecRad_robustfit_all(ixix)];
fprintf(1,'robustfit feedbacks ALL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_all);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_all([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_all([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_all([3 4 5 6])))

disp(' ')
robustfit_feedbacks_tropics = [xout.feedback.planck_ecRad_robustfit_tropics(ixix) xout.feedback.lapse_ecRad_robustfit_tropics(ixix) ...
                               xout.feedback.o3_ecRad_robustfit_tropics(ixix)     xout.feedback.wv_ecRad_robustfit_tropics(ixix) ...
                               xout.feedback.skt_ecRad_robustfit_tropics(ixix)    xout.feedback.ptemp_co2_ecRad_robustfit_tropics(ixix)];
fprintf(1,'robustfit feedbacks TRP : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_tropics);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_tropics([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_tropics([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_tropics([3 4 5 6])))

disp(' ')
robustfit_feedbacks_midlats = [xout.feedback.planck_ecRad_robustfit_midlats(ixix) xout.feedback.lapse_ecRad_robustfit_midlats(ixix) ...
                               xout.feedback.o3_ecRad_robustfit_midlats(ixix)     xout.feedback.wv_ecRad_robustfit_midlats(ixix) ...
                               xout.feedback.skt_ecRad_robustfit_midlats(ixix)    xout.feedback.ptemp_co2_ecRad_robustfit_midlats(ixix)];
fprintf(1,'robustfit feedbacks MID : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_midlats);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_midlats([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_midlats([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_midlats([3 4 5 6])))

disp(' ')
robustfit_feedbacks_poles = [xout.feedback.planck_ecRad_robustfit_poles(ixix) xout.feedback.lapse_ecRad_robustfit_poles(ixix) ...
                             xout.feedback.o3_ecRad_robustfit_poles(ixix)     xout.feedback.wv_ecRad_robustfit_poles(ixix) ...
                             xout.feedback.skt_ecRad_robustfit_poles(ixix)    xout.feedback.ptemp_co2_ecRad_robustfit_poles(ixix)];
fprintf(1,'robustfit feedbacks POL : planck lapse o3 wv skt tz/co2 = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  W/m2/K \n',robustfit_feedbacks_poles);
fprintf(1,'Jevanjee 2021 GRL : sum(robustfit_feedbacks_poles([1 2 3 4])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([1 2 3 4])))
fprintf(1,'Sergio Test       : sum(robustfit_feedbacks_poles([3 4 5 6])) = %5.2f W/m2/K \n',sum(robustfit_feedbacks_poles([3 4 5 6])))

xout.ecRadresults.polyfit_feedbacks           = polyfit_feedbacks;
xout.ecRadresults.robustfit_feedbacks_all     = robustfit_feedbacks_all;
xout.ecRadresults.robustfit_feedbacks_tropics = robustfit_feedbacks_tropics;
xout.ecRadresults.robustfit_feedbacks_midlats = robustfit_feedbacks_midlats;
xout.ecRadresults.robustfit_feedbacks_poles   = robustfit_feedbacks_poles;
xout.ecRadresults.dbins                       = dx;
xout.ecRadresults.nmeanval                    = nmeanval;
xout.ecRadresults.nstdval                     = nstdval;

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf; 
plot(indSST,-savenums,'.'); plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');

figure(4); clf
plot(1:64,xout.feedback.planck_ecRad_polyfit_latbin,...
     1:64,xout.feedback.lapse_ecRad_polyfit_latbin,...
     1:64,xout.feedback.o3_ecRad_polyfit_latbin,...
     1:64,xout.feedback.wv_ecRad_polyfit_latbin,...
     1:64,xout.feedback.skt_ecRad_polyfit_latbin,...
     1:64,xout.feedback.ptemp_co2_ecRad_polyfit_latbin,...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title('POLYFIT')

figure(5); clf
plot(1:64,xout.feedback.planck_ecRad_robustfit_latbin(:,1),...
     1:64,xout.feedback.lapse_ecRad_robustfit_latbin(:,1),...
     1:64,xout.feedback.o3_ecRad_robustfit_latbin(:,1),...
     1:64,xout.feedback.wv_ecRad_robustfit_latbin(:,1),...
     1:64,xout.feedback.skt_ecRad_robustfit_latbin(:,1),...
     1:64,xout.feedback.ptemp_co2_ecRad_robustfit_latbin(:,1),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title('ROBUSTFIT')

figure(6); clf
plot(1:64,xout.feedback.planck_ecRad_robustfit_latbin(:,2),...
     1:64,xout.feedback.lapse_ecRad_robustfit_latbin(:,2),...
     1:64,xout.feedback.o3_ecRad_robustfit_latbin(:,2),...
     1:64,xout.feedback.wv_ecRad_robustfit_latbin(:,2),...
     1:64,xout.feedback.skt_ecRad_robustfit_latbin(:,2),...
     1:64,xout.feedback.ptemp_co2_ecRad_robustfit_latbin(:,2),...
     'linewidth',2); 
plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10); ylabel('\lambda W/m2/K'); xlabel('Latbin'); title('ROBUSTFIT UNC')

figure(7); clf
factor = [0.2 0.2 1 1 0.2 1]'; factor = factor * ones(1,length(dx));
plot(dx,nmeanval(1:4,:),'o-','linewidth',2);         
  plotaxis2; hl= legend('planck','lapse','o3','wv','skt','t/co2','location','best','fontsize',10);             ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');
plot(dx,nmeanval(1:4,:).*factor(1:4,:),'o-','linewidth',2); 
  plotaxis2; hl= legend('0.2*Planck','0.2*Lapse','Ozone','WV','location','south','fontsize',10); ylabel('\delta OLR W/m2'); xlabel('\delta SKT K');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
