addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

for ii = 1 : 7
  figure(ii); clf
end

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin20_45angles.op.rtp');

if ~exist('rad')
  disp('loading in kcarta-->airs jacs and rads')
  fdir = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs_varyingscanang_oneprof';
  for ii = 1 : 45
    fin = [fdir '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
    a = load(fin);
    rad(ii,:) = a.rKc;
    fin = [fdir '/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_jac.mat'];
    a = load(fin);
    jac(ii,:,:) = a.rKc;
  end
  fKc = a.fKc;
  clear a
  rad = rad';
end
  
ix = find(fKc >= 1231.1,1);
figure(1); plot(p.scanang,rad2bt(fKc(ix),rad(ix,:)),'o-'); title('limb darkening, BT1231'); xlabel('scanang')

ix1 = find(fKc >= 790.1,1);
ix2 = find(fKc >= 791.5,1);
figure(2); plot(p.scanang,rad2bt(fKc(ix1),rad(ix1,:)),'bx-',p.scanang,rad2bt(fKc(ix2),rad(ix2,:)),'rx-'); 
title('limb darkening, BT (b)790 (r)792'); xlabel('scanang')
figure(3); plot(p.scanang,rad2bt(fKc(ix1),rad(ix1,:))-rad2bt(fKc(ix2),rad(ix2,:)),'o-'); 
title('limb darkening, BT 790-792'); xlabel('scanang')

ix = find(fKc >= 667,1);
figure(1); plot(p.scanang,rad2bt(fKc(ix),rad(ix,:)),'o-'); title('limb darkening, BT 667'); xlabel('scanang')

meanrad = nanmean(rad,2);
meanjac = squeeze(nanmean(jac,1));
ix = 01; WV_45ang = squeeze(jac(ix,:,(01:97)+0*97)); O3_45ang = squeeze(jac(ix,:,(01:97)+1*97)); T_45ang = squeeze(jac(ix,:,(01:97)+2*97));
         Wgt_45ang = squeeze(jac(ix,:,(01:97)+3*97)); stemp_45ang = squeeze(jac(ix,:,389)); 
ix = 22; WV_22ang = squeeze(jac(ix,:,(01:97)+0*97)); O3_22ang = squeeze(jac(ix,:,(01:97)+1*97)); T_22ang = squeeze(jac(ix,:,(01:97)+2*97));
         Wgt_22ang = squeeze(jac(ix,:,(01:97)+3*97)); stemp_22ang = squeeze(jac(ix,:,389)); 
ix = 45; WV_00ang = squeeze(jac(ix,:,(01:97)+0*97)); O3_00ang = squeeze(jac(ix,:,(01:97)+1*97)); T_00ang = squeeze(jac(ix,:,(01:97)+2*97));
         Wgt_00ang = squeeze(jac(ix,:,(01:97)+3*97)); stemp_00ang = squeeze(jac(ix,:,389)); 
WV_avg = meanjac(:,(01:97)+0*97); O3_avg = meanjac(:,(01:97)+1*97); T_avg = meanjac(:,(01:97)+2*97);
  Wgt_avg = meanjac(:,(01:97)+3*97); stemp_avg = meanjac(:,389); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotjacs 
run_sarta_jacs
compare_dofs
compare_dofs_certainchans

model_rads_jacs

