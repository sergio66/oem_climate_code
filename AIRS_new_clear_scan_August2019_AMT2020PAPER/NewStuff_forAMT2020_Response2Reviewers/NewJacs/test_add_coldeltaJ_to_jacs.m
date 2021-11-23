%function newJ = add_deltaJ_to_jacs(J,renorm,itimestep,iiBin)

i2834to2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat');
i2834to2645 = i2834to2645.ichan;

iiBin = 20;
itimestep = 89;

junk = num2str(itimestep,'%03d');
jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];

loader = ['jac0 = load(''' jacobian.filename ''');'];
eval(loader);

%jac0.qrenorm(1:6)
%jac0.qrenorm([10 110 210])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);    w1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_1.mat']);
a0 = load(['PBEFORE/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']); a1 = load(['PAFTER/jacresults' num2str(iiBin,'%02d') '_G_2356.mat']);
c0 = load(['PBEFORE/coljacresults' num2str(iiBin,'%02d') '.mat']);     c1 = load(['PAFTER/coljacresults' num2str(iiBin,'%02d') '.mat']);

qrenorm = ones(length(a0.f),1) * jac0.qrenorm;

jacERA  = [c0.jacCO2 c0.jacN2O c0.jacCH4 c0.jacCFC11 c0.jacCFC12 c0.jacST fliplr(w0.jacWV) fliplr(a0.jacT) fliplr(a0.jacO3)];  jacERA  = jacERA(i2834to2645,:);
jacUMBC = [c1.jacCO2 c1.jacN2O c1.jacCH4 c1.jacCFC11 c1.jacCFC12 c1.jacST fliplr(w1.jacWV) fliplr(a1.jacT) fliplr(a1.jacO3)];  jacUMBC = jacUMBC(i2834to2645,:);

[h,ha,pbefore,pa] = rtpread('pbeforeavg1_39.op.rtp');
[ppmvLay2] = layers2ppmv(h,pbefore,1:length(pbefore.stemp),2);
[ppmvLay4] = layers2ppmv(h,pbefore,1:length(pbefore.stemp),4);
[ppmvLay6] = layers2ppmv(h,pbefore,1:length(pbefore.stemp),6);

ppmvLay2 = ppmvLay2(90,:);
ppmvLay4 = ppmvLay4(90,:);
ppmvLay6 = ppmvLay6(90,:);

save ppmv_G2_4_6.mat ppmvLay*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(jac0.f - a0.f(i2834to2645))

figure(1); subplot(221); ii = 1; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*2.2/ppmvLay2(iiBin)*100,'r')
figure(1); subplot(222); ii = 2; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*1.0e-3/ppmvLay4(iiBin)*100,'r')
figure(1); subplot(223); ii = 3; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*5.0e-3/ppmvLay6(iiBin)*100,'r')
figure(1); subplot(224); ii = 6; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*0.1*100,'r')

figure(2); subplot(221); ii = 1; plot(jac0.f,jacUMBC(:,ii) ./ jacERA(:,ii));
figure(2); subplot(222); ii = 2; plot(jac0.f,jacUMBC(:,ii) ./ jacERA(:,ii));
figure(2); subplot(223); ii = 3; plot(jac0.f,jacUMBC(:,ii) ./ jacERA(:,ii));
figure(2); subplot(224); ii = 6; plot(jac0.f,jacUMBC(:,ii) - jacERA(:,ii));

figure(3); ind = (1:97) + 0*97; 
  subplot(211); plot(jac0.f,sum(squeeze(jac0.M_TS_jac_all(iiBin,:,ind))'),'b.-',jac0.f,sum(jacERA(:,ind)')*0.01,'r')
  subplot(212); plot(jac0.f,sum(jacUMBC(:,ind)') ./ sum(jacERA(:,ind)'),'r'); ax = axis; ax(3) = 0.9; ax(4) = 1.1; axis(ax)
figure(4); ind = (1:97) + 1*97; 
  subplot(211); plot(jac0.f,sum(squeeze(jac0.M_TS_jac_all(iiBin,:,ind))'),'b.-',jac0.f,sum(jacERA(:,ind)')*0.01,'r')
  subplot(212); plot(jac0.f,sum(jacUMBC(:,ind)') ./ sum(jacERA(:,ind)'),'r'); ax = axis; ax(3) = 0.9; ax(4) = 1.1; axis(ax)
figure(5); ind = (1:97) + 2*97; 
  subplot(211); plot(jac0.f,sum(squeeze(jac0.M_TS_jac_all(iiBin,:,ind))'),'b.-',jac0.f,sum(jacERA(:,ind)')*0.01,'r')
  subplot(212); plot(jac0.f,sum(jacUMBC(:,ind)') ./ sum(jacERA(:,ind)'),'r'); ax = axis; ax(3) = 0.9; ax(4) = 1.1; axis(ax)

ratio = (jacUMBC ./ jacERA);
ratio(ratio < 0.75) = 0.75;
ratio(ratio > 1.25) = 1.25;

oldjac = squeeze(jac0.M_TS_jac_all(iiBin,:,:));
newjac = squeeze(jac0.M_TS_jac_all(iiBin,:,:));
ii = 1; newjac(:,ii) = newjac(:,ii) .* ratio(:,ii);
ii = 2; newjac(:,ii) = newjac(:,ii) .* ratio(:,ii);
ii = 3; newjac(:,ii) = newjac(:,ii) .* ratio(:,ii);
ii = 6; newjac(:,ii) = newjac(:,ii) .* ratio(:,ii);

ind = 7:297; newjac(:,ind) = newjac(:,ind) .* ratio(:,ind);

figure(6); clf; 
pcolor(jac0.f,1:297,newjac'./oldjac'); shading flat; colorbar; colormap(usa2)
axis([640 1640 1 297]);

figure(7); subplot(221); ii = 4; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*5.0e-3/ppmvLay6(iiBin)*100,'r')
figure(7); subplot(222); ii = 5; plot(jac0.f,squeeze(jac0.M_TS_jac_all(iiBin,:,ii)),'b.-',jac0.f,jacERA(:,ii)*5.0e-3/ppmvLay6(iiBin)*100,'r')
figure(7); subplot(223); ii = 4; plot(jac0.f,jacUMBC(:,ii) ./ jacERA(:,ii));
figure(7); subplot(224); ii = 5; plot(jac0.f,jacUMBC(:,ii) ./ jacERA(:,ii));
