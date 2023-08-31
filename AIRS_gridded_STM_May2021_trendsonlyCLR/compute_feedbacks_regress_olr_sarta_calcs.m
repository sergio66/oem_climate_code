function xout = compute_feedbacks_regress_sarta_calcs(x0,indSST,iLambda_UseGlobalSST_regress,h)

%% change radiance mW --> W and then multiply by pi for flux
%% THIS IS FOR 2645 AIRS only!!!!

ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

xout = x0;

if isfield(xout,'feedback')
  xout = rmfield(xout,'feedback');
end
if isfield(xout,'feedback_sarta')
  xout = rmfield(xout,'feedback_sarta');
end
%if isfield(xout,'ecRadresults')
%  xout = rmfield(xout,'ecRadresults');
%end

% junk = pi/1000*sum(xout.planck - xout.olr0,1);
% junk = -junk./indSST;
% xout.feedback_sarta.planck = junk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
  YY = Y(:)';
coslat  = cos(YY*pi/180);

xout.feedback_sarta.global_coslat_wgt_skt  = sum(indSST .* coslat)/sum(coslat);

globalSST = nanmean(indSST);
globalSST = xout.feedback_ecRad.global_coslat_wgt_skt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk1 = pi/1000*trapz(h.vchan(ix1),xout.planck(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.planck(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback_sarta.planck = junk;
good = abs(junk) < 10;
xout.feedback_sarta.planck_nanmean_global = nanmean(junk(good));
bunk = polyfit(indSST(good),-junk12(good),1);
xout.feedback_sarta.planck_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.o3(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.o3(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback_sarta.o3 = junk;
xout.feedback_sarta.o3_nanmean_global = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback_sarta.o3_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.ptemp_co2(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.ptemp_co2(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback_sarta.ptemp_co2 = junk;
xout.feedback_sarta.ptemp_co2_nanmean_global = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback_sarta.ptemp_co2_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.lapse(ix1,:) - xout.planck(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.lapse(ix2,:) - xout.planck(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback_sarta.lapse = junk;
xout.feedback_sarta.lapse_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback_sarta.lapse_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.wv(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.wv(ix2,:) - xout.olr0(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback_sarta.wv = junk;
xout.feedback_sarta.wv_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback_sarta.wv_polyfit_global = bunk(1);

xout.feedback_sarta.wv = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),xout.skt(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.skt(ix2,:) - xout.olr0(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback_sarta.skt = junk;
xout.feedback_sarta.skt_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback_sarta.skt_polyfit_global = bunk(1);
