%% /home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res quickuse.nml rad.dat jac.dat
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

caPFname = '/home/sergio/KCARTA/IP_PROFILES/junk49.op.rtp';
[h,ha,p,pa] = rtpread(caPFname);
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,p,1,2);  %% so this is for 370 ppmv
pN = p.plevs(1:100,:)-p.plevs(2:101,:);
pD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
plays = pN./pD;
plays = flipud(plays);

[rad,w] = readkcstd('rad.dat');
[jac,w] = readkcjac('jac.dat');

addpath /home/sergio/MATLABCODE
[fc,jc] = quickconvolve(w,jac(:,1:97),0.5,0.5);

addpath /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/
airs_convolve_file_numchans
w = 1 : length(rad);
w = (w-1)*0.0025 + 605;
[fc,rc] = convolve_airs(w,rad,clist,sfile);
[fc,jc] = convolve_airs(w,jac,clist,sfile);

addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR
lstrow = load('sarta_chans_for_l1c.mat');
fc = fc(lstrow.ichan);
rc = rc(lstrow.ichan);
jc = jc(lstrow.ichan,1:97) * 2.2/370;
whos jc

figure(1); pcolor(fc,1:97,jc'); shading flat; colorbar
axis([640 940 0 97])

ix = find(fc >= 700);
figure(2); semilogy(mean(jc),plays(1:97),'b',mean(jc(ix,1:97)),plays(1:97),'r'); grid
xlabel('dBT/dq (K/2.2ppmv)')
set(gca','ydir','reverse'); axis([-0.6 +0.6 1e1 1000]);

figure(3); plot(fc,sum(jc')); title('col jac')
axis([640 940 -0.1 +0.1])

%{
topts.numchan = 2645;
topts.chan_LW_SW =  0;  %% just LW/MW DEFAULT, keep strat T chans                       470 chans
%% topts.chan_LW_SW = -1;  %% all chans,LW/MW and SW, keep strat T chans <<<<<<<           769 chans
%% topts.chan_LW_SW = +1;  %% LW/MW but avoid deep 15 um                                   436 chans
%% topts.chan_LW_SW = +2;  %% LW/MW/SW but avoid deep 15 um not yet coded                  735 chans
%% topts.chan_LW_SW = -2;  %% SW only                                                      299 chans

-rw-rw-r-- 1 sergio pi_strow 291177 Dec 15 16:52 SAVE_chanselections/chan_LW_SW_m2.mat
-rw-rw-r-- 1 sergio pi_strow 513483 Dec 15 16:51 SAVE_chanselections/chan_LW_SW_p2.mat
-rw-rw-r-- 1 sergio pi_strow 359143 Dec 15 16:50 SAVE_chanselections/chan_LW_SW_p1.mat
-rw-rw-r-- 1 sergio pi_strow 538791 Dec 15 16:50 SAVE_chanselections/chan_LW_SW_m1.mat
-rw-rw-r-- 1 sergio pi_strow 376509 Dec 15 16:48 SAVE_chanselections/chan_LW_SW_0.mat
%}

allLWchans = load('../SAVE_chanselections/chan_LW_SW_0.mat');
  allLWchans = allLWchans.driver.jacobian.chanset;
nostratLWchans = load('../SAVE_chanselections/chan_LW_SW_p1.mat');
  nostratLWchans = nostratLWchans.driver.jacobian.chanset;

iall = allLWchans;
ix = nostratLWchans;

figure(2); semilogy(mean(jc(iall,1:97)),plays(1:97),'b+-',mean(jc(ix,1:97)),plays(1:97),'r+-'); grid
xlabel('dBT/dq (K/2.2ppmv)'); ylabel('p(mb)')
set(gca','ydir','reverse'); axis([-0.005 +0.0005 1e1 1000]);
hl = legend('Channels used','Remove channels < 700 cm^-1','location','best'); set(hl,'fontsize',10);

figure(2); semilogy(sum(jc(iall,1:97)),plays(1:97),'b+-',sum(jc(ix,1:97)),plays(1:97),'r+-'); grid
xlabel('dBT/dq (K/2.2ppmv)'); ylabel('p(mb)')
set(gca','ydir','reverse'); axis([-0.1 +0.1 1e1 1000]);
hl = legend('Channels used','Remove channels < 700 cm^-1','location','best'); set(hl,'fontsize',10);

figure(2); semilogy(370/2.2*mean(jc(iall,1:97)),plays(1:97),'b+-',370/2.2*mean(jc(ix,1:97)),plays(1:97),'r+-'); grid
xlabel('dBT/dq (K)'); ylabel('p(mb)')
set(gca','ydir','reverse'); axis([-0.005 +0.0005 1e1 1000]);
hl = legend('Channels used','Remove channels < 700 cm^-1','location','best'); set(hl,'fontsize',10);

%{
aslprint('co2jac.png',1);
%}
