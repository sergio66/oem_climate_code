function newJ = add_deltaJ_to_jacs(J,iiBin,iAllorOnlyColG_stemp)

iAllorOnlyColG_stemp = 0
iAllorOnlyColG_stemp = 0
iAllorOnlyColG_stemp = 0

if nargin == 2
  iAllorOnlyColG_stemp = +1;   %% change all (except CFC = index 4,5)
                               %% -1 : don't change WV(z),T(z),O3(z) ie only change first few
                               %%  0 : change nothing : ratio = 1
end

%{
iiBin = 15;
jac0 = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_089_M_TS_jac_all_5_97_97_97_2645.mat');
J = squeeze(jac0.M_TS_jac_all(iiBin,:,:)); renorm = jac0.qrenorm;
newJ = add_deltaJ_to_jacs(J,iiBin);
%}

load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/NewJacs_forAMT2020/ppmv_G2_4_6.mat
i2834to2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat');
i2834to2645 = i2834to2645.ichan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir0 = ['/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/NewJacs_forAMT2020/'];

if iiBin == 40
  %% oops, only did Dec of all years ... so never did norhtern most (latbin 40)
  iiBin = 39;
end

w0 = load([dir0 'PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);    w1 = load([dir0 'PAFTER/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);
a0 = load([dir0 'PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']); a1 = load([dir0 'PAFTER/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']);
c0 = load([dir0 'PBEFORE/coljacresults' num2str(iiBin,'%02d') '.mat']);     c1 = load([dir0 'PAFTER/coljacresults' num2str(iiBin,'%02d') '.mat']);

%jacERA  = [a0.jacCO2 a0.jacN2O a0.jacCH4 zeros(size(a0.jacCH4)) zeros(size(a0.jacCH4)) a0.jacST w0.jacWV a0.jacT a0.jacO3];  jacERA  = jacERA(i2834to2645,:);
%jacUMBC = [a1.jacCO2 a1.jacN2O a1.jacCH4 zeros(size(a1.jacCH4)) zeros(size(a1.jacCH4)) a1.jacST w1.jacWV a1.jacT a1.jacO3];  jacUMBC = jacUMBC(i2834to2645,:);
jacERA  = [c0.jacCO2 c0.jacN2O c0.jacCH4 c0.jacCFC11 c0.jacCFC12 c0.jacST fliplr(w0.jacWV) fliplr(a0.jacT) fliplr(a0.jacO3)];  jacERA  = jacERA(i2834to2645,:);
jacUMBC = [c1.jacCO2 c1.jacN2O c1.jacCH4 c1.jacCFC11 c1.jacCFC12 c1.jacST fliplr(w1.jacWV) fliplr(a1.jacT) fliplr(a1.jacO3)];  jacUMBC = jacUMBC(i2834to2645,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ratio = (jacUMBC ./ jacERA);
ratio(ratio < 0.75) = 0.75;
ratio(ratio > 1.25) = 1.25;
bad = find(isnan(ratio) | isinf(ratio));
ratio(bad) = 1.0;

if iAllorOnlyColG_stemp == 0
  ratio = ones(size(ratio));
end

oldJ = J;
newJ = J;

%% do trace gases
ii = 1; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);
ii = 2; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);
ii = 3; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);
ii = 4; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);
ii = 5; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);

if iAllorOnlyColG_stemp > 0
  %% do stemp, WV(z),T(z),O3(z)
  ii = 6; newJ(:,ii) = newJ(:,ii) .* ratio(:,ii);
  ind = 7:297; newJ(:,ind) = newJ(:,ind) .* ratio(:,ind);
end
