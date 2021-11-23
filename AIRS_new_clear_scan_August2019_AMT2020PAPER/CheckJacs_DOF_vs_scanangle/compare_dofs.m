addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /home/sergio/MATLABCODE/AveragingKernel_DOF
addpath /home/sergio/MATLABCODE/NeDT_NeDN_AIRS_IASI_CRIS

dWV = 0.1479; dT = 2; dO3 = 0.175; %% use 2K 40% 50% for errors in T,WV,O3
dWV = 0.0792; dT = 2; dO3 = 0.175; %% use 1K 20% 50% for errors in T,WV,O3

g = dogoodchan;
[finstr,nedT] = getBTnoise_instr(fKc(1:2378),meanrad(1:2378),'airs');
Se = diag(nedT(g)); Se = Se.*Se;
iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
  [ak_WV_airs,dofs_WV_airsAVG] = generic_ak_dof(finstr(g),Se,Sa,meanjac(g,iInd),'WV');
iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
  [ak_O3_airs,dofs_O3_airsAVG] = generic_ak_dof(finstr(g),Se,Sa,meanjac(g,iInd),'O3');
iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
  [ak_T_airs,dofs_T_airsAVG] = generic_ak_dof(finstr(g),Se,Sa,meanjac(g,iInd),'T');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 45
  [finstr,nedT] = getBTnoise_instr(fKc(1:2378),rad(1:2378,ii),'airs');
  Se = diag(nedT(g)); Se = Se.*Se;
  iInd = (1:97)+0*97; Sa = dWV * ones(1,97); Sa = diag(Sa.^2);   %% exp10(0.0792)-1 = 0.2 = 20% error
    [ak_WV_airs,dofs_WV_airs(ii)] = generic_ak_dof(finstr(g),Se,Sa,squeeze(jac(ii,g,iInd)),'WV');
  iInd = (1:97)+1*97; Sa = dO3 * ones(1,97); Sa = diag(Sa.^2);      %% 1 K error
    [ak_O3_airs,dofs_O3_airs(ii)] = generic_ak_dof(finstr(g),Se,Sa,squeeze(jac(ii,g,iInd)),'O3');
  iInd = (1:97)+2*97; Sa = dT * ones(1,97); Sa = diag(Sa.^2);    %% exp10(0.175)-1 = 0.2 = 50% error
    [ak_O3_airs,dofs_T_airs(ii)] = generic_ak_dof(finstr(g),Se,Sa,squeeze(jac(ii,g,iInd)),'T');
end

for ii = 1 : 5
  figure(ii); clf
end

for ii = 1 : 45
  dofs2_WV_airs(ii) = dofs_WV_airs(ii).dofs;
  dofs2_O3_airs(ii) = dofs_O3_airs(ii).dofs;
  dofs2_T_airs(ii)  = dofs_T_airs(ii).dofs;
end

figure(1); plot(p.scanang,dofs2_WV_airs,'b.-',-22,dofs_WV_airsAVG.dofs,'ro'); title('WV dof'); grid
  xlabel('scanang'); hl = legend('jac(scanang)','avg jac','location','best'); 
figure(2); plot(p.scanang,dofs2_O3_airs,'b.-',-22,dofs_O3_airsAVG.dofs,'ro'); title('O3 dof'); grid
  xlabel('scanang'); hl = legend('jac(scanang)','avg jac','location','best'); 
figure(3); plot(p.scanang,dofs2_T_airs,'b.-',-22,dofs_T_airsAVG.dofs,'ro'); title('T dof'); grid
  xlabel('scanang'); hl = legend('jac(scanang)','avg jac','location','best'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iInd = (1:97)+3*97;  %% wgt fcn
wx45deg = squeeze(jac(1,g,iInd));    wx00deg = squeeze(jac(45,g,iInd));      wx22deg = squeeze(jac(22,g,iInd));

iChan = find(fKc >= 720 & fKc <= 800); iChan = iChan(iChan < 500); iChan = iChan(1:10:length(iChan)); whos iChan
figure(4); plot(wx00deg(iChan,:),1:97,'b',wx45deg(iChan,:),1:97,'r'); title('15 um wgt (b) 0 deg (r) 45 deg')
  ylim([0 50]);

iChan = find(fKc >= 1320 & fKc <= 1620); iChan = iChan(iChan < 1500); iChan = intersect(iChan,g); iChan = iChan(1:10:length(iChan)); whos iChan
figure(5); plot(wx00deg(iChan,:),1:97,'b',wx45deg(iChan,:),1:97,'r'); title('6.7um um wgt (b) 0 deg (r) 45 deg')
  ylim([0 50]);

