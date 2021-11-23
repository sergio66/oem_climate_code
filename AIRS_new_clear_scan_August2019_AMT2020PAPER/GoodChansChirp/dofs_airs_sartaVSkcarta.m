%{

1.10 pm
For AIRS L1c, see
/home/motteler/shome/airs_decon/test/int_L1c_NEdN.mat.  This is
representative AIRS NeN (aka NEdN) averaged over a full day.  mNeN is
the raw mean, and nedn the mNeN values interpolated for the synthetic
channels.  The latter is what we've mostly been using, for example as
input for the AIRS parent CHIRP NEdN calc.

1:18 pm
For CrIS, see the same directory, and various files N20_HR_NEdN.mat,
cris_HR_NEdN.mat, cris_NEdN.mat, NPP_HR_NEdN.mat, cris_LR_NEdN.mat.
For CrIS, the NEdN is going to depend on both the resolution and the
apodization.

motteler  1:25 PM
For CHIRP, maybe the easiest thing is to see the little script
quick_look.m in /home/motteler/shome/chirp_test.  This plots AIRS and
CrIS parent NEdN.

sergio(opens in new tab)  6:03 PM
Howard, I notice you gave globally averaged NeDN .. so I presume I
should run sarta for a tropical profile, 22 deg scan angle (since that
is average profile over the earth, 50% is the tropics, and 22deg is
the mean scanagle), compute BT .. then add NeDN, compute new BT .. and
that would be mean NeDT?  I do have those dBT/drBT codes .. but since
these are all averages, they might be overkill compared to what I
sketched above? (edited)

motteler(opens in new tab)  6:52 PM
NEdN is mostly a property of the instrument, not the particular scene
you are looking at.  NEdT is a function of the scene, if you have dr
and r it is something like dt = B(v, r + dr) - B(v, r).  So don't
worry about NEdN for particular regions.  Averaging over a whole day
(as I was doing) will give you a good first estimate.  NEdN might vary
slightly with orbital phase, but that's a different concern.  And we
do calculate NEdN from the ICT, and so in effect at a particular
nominal temperature.  But don't worry about that for now, just use
NEdN as we have it, and calculate NEdT for whatever scenes you want.
New

motteler(opens in new tab)  7:12 PM
/home/motteler/shome/airs_decon/test/nedt_test1.m has some examples of
calculating NEdT for AIRS, AIRS to CrIS, and true Cris.

sergio(opens in new tab)  10:05 PM
Hi Howard, but of a messy plot but I've been able to get what I need
from looking at your code and ,mat files, thanks!  I did notice
AIRS2CHIRP and CRIS2CHIRP had different noise, but I guess that is
solely because AIRS and CrIS FSR have different noise anyway 

motteler(opens in new tab)  11:14 PM
Right, the translations start out with the noise functions of the
parents, modified by the translation, for AIRS, and by interpolation
and apodization for CrIS.  Those steps reduce noise, to varying
degrees.  AIRS L1c noise is complicated.  It is only really defined
for the channels with a direct L1b parent, in which case it is just
the L1b noise.  For the synthetic channels we are simply interpolating
NEdN in frequency.  That's only a rough approximation, since noise can
vary so much from channel to channel, but at least it preserves
properties of modules, to some extent.

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /home/sergio/MATLABCODE/AveragingKernel_DOF
addpath /home/sergio/MATLABCODE/NeDT_NeDN_AIRS_IASI_CRIS

for ii = 1 : 40
  rad20 = load(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat']);
  kcrad(:,ii) = rad20.rKc;
end

rad20 = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs/individual_prof_convolved_kcarta_crisHI_crisMED_20.mat');
jac20 = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_20_jac.mat');

sarta20 = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_CRIS_IASI_allres_WV_T_O3_stemp_jacs/sarta_finitediff_jacs_prof20.mat'); %% airs2834 jacs

d40dbt1 = rad2bt(rad20.fKc,kcrad) - rad2bt(rad20.fKc,sarta20.rads1_40);
d40dbt2 = rad2bt(rad20.fKc,kcrad) - rad2bt(rad20.fKc,sarta20.rads2_40);

plot(rad20.fKc,nanmean(d40dbt1'),'b',rad20.fKc,nanstd(d40dbt1'),'c--',...
     rad20.fKc,nanmean(d40dbt2'),'r',rad20.fKc,nanstd(d40dbt2'),'m--')

plot(rad20.fKc,rms(d40dbt1'),'b',rad20.fKc,rms(d40dbt2'),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airs = load('/home/motteler/shome/airs_decon/test/int_L1c_NEdN.mat');
crishr1 = load('/home/motteler/shome/airs_decon/test/N20_HR_NEdN.mat');
crishr2 = load('/home/motteler/shome/airs_decon/test/cris_HR_NEdN.mat');
cris3 = load('/home/motteler/shome/airs_decon/test/cris_NEdN.mat');
cris4 = load('/home/motteler/shome/airs_decon/test/NPP_HR_NEdN.mat');
crislo1 = load('/home/motteler/shome/airs_decon/test/cris_LR_NEdN.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AIRS noise
addpath /home/sergio/MATLABCODE
plot(airs.v1c,airs.nedn,'b.',airs.v1c,airs.mNeN,'r',instr_chans,instr_chans('airs',2),'k')

load ../../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/sarta_chans_for_l1c.mat
plot(rad20.fKc(ichan) - airs.v1c) %% should be tiny

airsNedT = dr2dt(rad20.fKc(ichan),rad20.rKc(ichan),airs.nedn);

nedt250orig = interp1(instr_chans,instr_chans('airs',2),rad20.fKc(ichan),[],'extrap');
airs_noiseTtrue_from_250 = nedt_T0_T1(rad20.fKc(ichan),nedt250orig,250*ones(2645,1),real(rad2bt(rad20.fKc(ichan),rad20.rKc(ichan))));

figure(1); clf
plot(airs.v1c,airsNedT,'r',instr_chans,instr_chans('airs',2),'k',airs.v1c,airs_noiseTtrue_from_250,'b')
  hl = legend('Howard <NeDN> -> NeDT for tropical profile','AIRS 2002 NeDT for 250 K','AIRS 2002 NeDT for 250 K -> NeDT for tropical profile','location','best','fontsize',8);
  axis([640 1640 0 1])
 title('AIRS NeDT');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHIRP noise
load chirpNeDN.mat

figure(2); clf
[Ychirp,iXchirp,iYchirp] = intersect(chirp.wnum,rad20.med_fcris);
chirpNedTA = dr2dt(chirp.wnum,rad20.med_rcris_all(iYchirp),chirp.nednA);
chirpNedTC = dr2dt(chirp.wnum,rad20.med_rcris_all(iYchirp),chirp.nednC);
plot(chirp.wnum,chirpNedTA,chirp.wnum,chirpNedTC); 
  hl = legend('AIRS2CHIRP','CRIS2CHIRP','location','best','fontsize',8);
  axis([640 1640 0 1])
 title('CHIRP NeDT');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CRIS FSR noise
figure(3); clf
[YcrisFSR,iXcrisFSR,iYcrisFSR] = intersect(crishr2.freq,rad20.hi_fcris);
crisFSR_NeDT = dr2dt(crishr2.freq,rad20.hi_rcris_all(iYcrisFSR),mean(crishr2.nedn,2));
plot(crishr2.freq,crisFSR_NeDT);
  title('CRIS FSR NeDT');
  axis([640 1640 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
plot(airs.v1c,airsNedT,'b',crishr2.freq,crisFSR_NeDT,'g',chirp.wnum,chirpNedTA,'r',chirp.wnum,chirpNedTC,'m')
hl = legend('AIRS L1C','CRIS FSR','AIRS2CHIRP','CRIS2CHIRP','location','best','fontsize',10);
  axis([640 1640 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); 
dWV = 0.0792; dT = 1; dO3 = 0.175; %% use 1K 20% 50% for errors in T,WV,O3
dWV = 0.1479; dT = 2; dO3 = 0.175; %% use 2K 40% 50% for errors in T,WV,O3

iChanAirs = 1 : length(airs.v1c);

fprintf(1,'\n AIRS L1C dofs KCARTA : \n')
Se = diag(airsNedT.^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,jac20.rKc(ichan,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,jac20.rKc(ichan,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,jac20.rKc(ichan,iInd),'O3');

ind = (1:97)+2*97; sarta2008(:,ind) = fliplr(sarta20.Tz1jac);
ind = (1:97)+0*97; sarta2008(:,ind) = fliplr(sarta20.WV1jac);
ind = (1:97)+1*97; sarta2008(:,ind) = fliplr(sarta20.O31jac);
fprintf(1,'\n AIRS L1C dofs SARTA 2008 : \n')
Se = diag((airsNedT).^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'O3');

ind = (1:97)+2*97; sarta2018(:,ind) = fliplr(sarta20.Tz2jac);
ind = (1:97)+0*97; sarta2018(:,ind) = fliplr(sarta20.WV2jac);
ind = (1:97)+1*97; sarta2018(:,ind) = fliplr(sarta20.O32jac);
fprintf(1,'\n AIRS L1C dofs SARTA 2018 : \n')
Se = diag((airsNedT).^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'O3');

ind = (1:97)+2*97; sarta2008(:,ind) = fliplr(sarta20.Tz1jac);
ind = (1:97)+0*97; sarta2008(:,ind) = fliplr(sarta20.WV1jac);
ind = (1:97)+1*97; sarta2008(:,ind) = fliplr(sarta20.O31jac);
fprintf(1,'\n AIRS L1C dofs SARTA 2008 + dBT1: \n')
dBT1 = abs(rad2bt(rad20.fKc,rad20.rKc)-rad2bt(rad20.fKc,sarta20.rads1));
dBT1 = abs(rad2bt(rad20.fKc,rad20.rKc)-sarta20.rads1); %% oops already in BT
dBT1 = rms(d40dbt1')';
Se = diag((airsNedT+dBT1(ichan)).^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2008(ichan,iInd),'O3');

ind = (1:97)+2*97; sarta2018(:,ind) = fliplr(sarta20.Tz2jac);
ind = (1:97)+0*97; sarta2018(:,ind) = fliplr(sarta20.WV2jac);
ind = (1:97)+1*97; sarta2018(:,ind) = fliplr(sarta20.O32jac);
fprintf(1,'\n AIRS L1C dofs SARTA 2018 + dBT2: \n')
dBT2 = abs(rad2bt(rad20.fKc,rad20.rKc)-rad2bt(rad20.fKc,sarta20.rads2));
dBT2 = abs(rad2bt(rad20.fKc,rad20.rKc)-sarta20.rads2); %% oops already in BT
dBT2 = rms(d40dbt2')';
Se = diag((airsNedT+dBT2(ichan)).^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(airs.v1c,Se,Sa,sarta2018(ichan,iInd),'O3');

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n CRIS FSR dofs : \n')
Se = diag(crisFSR_NeDT.^2);
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_crisFSR,dofs_T_crisFSR_AVG] = generic_ak_dof(crishr2.freq,Se,Sa,jac20.hi_rcris_all(iYcrisFSR,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_crisFSR,dofs_WV_crisFSR_AVG] = generic_ak_dof(crishr2.freq,Se,Sa,jac20.hi_rcris_all(iYcrisFSR,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_crisFSR,dofs_O3_crisFSR_AVG] = generic_ak_dof(crishr2.freq,Se,Sa,jac20.hi_rcris_all(iYcrisFSR,iInd),'O3');

fprintf(1,'\n CHIRP A dofs : \n')
bad = find(isinf(chirpNedTA) | isnan(chirpNedTA));
good = find(isfinite(chirpNedTA));
if length(bad) > 0
  fprintf(1,'ohoh found %6i bad numbers for CHIRP NeDT! interpolating \n',length(bad));  
  chirpNedTA(bad) = interp1(chirp.wnum(good),chirpNedTA(good),chirp.wnum(bad),[],'extrap');
end
Se = diag(chirpNedTA.^2);
bad = find(isnan(jac20.med_rcris_all) | isinf(jac20.med_rcris_all));
if length(bad) > 0
  fprintf(1,'ohoh found %6i bad jacs for CHIRP! setting to 0 \n',length(bad));
  jac20.med_rcris_all(bad) = 0;
end
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_chirpA,dofs_T_chirpA_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_chirpA,dofs_WV_chirpA_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_chirpA,dofs_O3_chirpA_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'O3');

fprintf(1,'\n CHIRP C dofs : \n')
bad = find(isinf(chirpNedTC) | isnan(chirpNedTC));
good = find(isfinite(chirpNedTC));
if length(bad) > 0
  fprintf(1,'ohoh found %6i bad numbers for CHIRP NeDT! interpolating \n',length(bad));  
  chirpNedTC(bad) = interp1(chirp.wnum(good),chirpNedTC(good),chirp.wnum(bad),[],'extrap');
end
Se = diag(chirpNedTC.^2);
bad = find(isnan(jac20.med_rcris_all) | isinf(jac20.med_rcris_all));
if length(bad) > 0
  fprintf(1,'ohoh found %6i bad jacs for CHIRP! setting to 0 \n',length(bad));
  jac20.med_rcris_all(bad) = 0;
end
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_T_chirpC,dofs_T_chirpC_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'T');
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_chirpC,dofs_WV_chirpC_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_O3_chirpC,dofs_O3_chirpC_AVG] = generic_ak_dof(chirp.wnum,Se,Sa,jac20.med_rcris_all(iYchirp,iInd),'O3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); plot(dofs_WV_airsAVG.cdofs,1:97,'b',dofs_WV_crisFSR_AVG.cdofs,1:97,'g',dofs_WV_chirpA_AVG.cdofs,1:97,'r',dofs_WV_chirpC_AVG.cdofs,1:97,'m','linewidth',2);
  hl = legend('AIRS L1C','CRIS FSR','AIRS 2 CHIRP','CRIS 2 CHIRP','location','best','fontsize',10);
  title('WV DOF')
figure(2); plot(dofs_O3_airsAVG.cdofs,1:97,'b',dofs_O3_crisFSR_AVG.cdofs,1:97,'g',dofs_O3_chirpA_AVG.cdofs,1:97,'r',dofs_O3_chirpC_AVG.cdofs,1:97,'m','linewidth',2);
  hl = legend('AIRS L1C','CRIS FSR','AIRS 2 CHIRP','CRIS 2 CHIRP','location','best','fontsize',10);
  title('O3 DOF')
figure(3); plot(dofs_T_airsAVG.cdofs,1:97,'b',dofs_T_crisFSR_AVG.cdofs,1:97,'g',dofs_T_chirpA_AVG.cdofs,1:97,'r',dofs_T_chirpC_AVG.cdofs,1:97,'m','linewidth',2);
  hl = legend('AIRS L1C','CRIS FSR','AIRS 2 CHIRP','CRIS 2 CHIRP','location','best','fontsize',10);
  title('T DOF')

