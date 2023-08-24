function xout = compute_feedbacks_umbc_ecRad_regress_calcs_olr(x0,indSST,iLambda_UseGlobalSST_regress,h)

%% change radiance mW --> W and then multiply by pi for flux
%% THIS IS FOR 2645 AIRS only!!!!

ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

xout = x0;

%junk = pi/1000*sum(xout.planck - xout.olr0,1);
%junk = -junk./indSST;
%xout.feedback.planck = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),xout.planck(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.planck(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback.planck = junk;
good = abs(junk) < 10;
xout.feedback.planck_nanmean_global = nanmean(junk(good));
bunk = polyfit(indSST(good),-junk12(good),1);
xout.feedback.planck_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.o3(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.o3(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback.o3 = junk;
xout.feedback.o3_nanmean_global = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback.o3_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.ptemp_co2(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.ptemp_co2(ix2,:) - xout.olr0(ix2,:));
junk12 = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk12./indSST;
else
  junk = -junk12/globalSST;
end
xout.feedback.ptemp_co2 = junk;
xout.feedback.ptemp_co2_nanmean_global = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback.ptemp_co2_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.lapse(ix1,:) - xout.planck(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.lapse(ix2,:) - xout.planck(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback.lapse = junk;
xout.feedback.lapse_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback.lapse_polyfit_global = bunk(1);

junk1 = pi/1000*trapz(h.vchan(ix1),xout.wv(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.wv(ix2,:) - xout.olr0(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback.wv = junk;
xout.feedback.wv_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback.wv_polyfit_global = bunk(1);

xout.feedback.wv = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),xout.skt(ix1,:) - xout.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),xout.skt(ix2,:) - xout.olr0(ix2,:));
junk = junk1 + junk2;
if iLambda_UseGlobalSST_regress == -1
  junk = -junk./indSST;
else
  junk = -junk/globalSST;
end
xout.feedback.skt = junk;
xout.feedback.skt_nanmean_gobal = nanmean(junk);
bunk = polyfit(indSST,-junk12,1);
xout.feedback.skt_polyfit_global = bunk(1);
