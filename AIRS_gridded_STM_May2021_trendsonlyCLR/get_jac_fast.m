function [m_ts_jac,nlays,qrenorm,freq2645,colo3,profilejunk] = get_jac_fast(fname,iWhichLatLonBin,iLonBin,iLatBin,iVersJac,iOldORNew,topts);

iDoSlow = -1;
if iDoSlow > 0
  tic;
  this_is_slow_jac_load  %% all it needs is fname which is set in strow_override_defaults_latbins_AIRS_fewlays.m (ie does not really need iVersJac)
  'done1'
  timervalx1 = toc;
end

if iDoSlow > 0
  tic;
end

if iVersJac == 2019
  see_clust_put_together_jacs_clr
elseif iVersJac == 2022
  %% this one gives colo3
  see_clust_put_together_jacs_cldERA5_2022
elseif iVersJac == 2021
  see_clust_put_together_jacs_clrERA5_2021
elseif iVersJac == 2012 | iVersJac == 2015
  see_clust_put_together_jacs_clrERA5_2012
elseif iVersJac == 2014
  see_clust_put_together_jacs_clrERA5_2014
else
  iVersJac
  error('iVersJac = 2012,2014,2019,2021,2022 only')
end

if ~exist('colo3')
 colo3 = zeros*freq2645;
end

if iDoSlow > 0
  'done2'
  timervalx2 = toc;

  timervalx1
  timervalx2

  [mm,nnS] = size(m_ts_jac_slow);
  [mm,nnF] = size(m_ts_jac_fast);

  if nnS == nnF
    plot(sum(m_ts_jac_slow-m_ts_jac_fast),'o-'); pause(0.1); grid
  end
else
  [mm,nnF] = size(m_ts_jac_fast);
end
qrenorm = ones(1,nnF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see ../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/compute_sarta_cld_jacobians_onedelta.m
% dq_CO2_nominal = (tx2p(:,j2) - t0)*(2.2/370) /(log(1+normer.dQ));  %% for 2.2/370 frac change
% dq_N2O_nominal = (tx2p(:,j4) - t0)*(1/300) /(log(1+normer.dQ));    %% for 1/300   frac change
% dq_CH4_nominal = (tx2p(:,j6) - t0)*(5/1860) /(log(1+normer.dQ));   %% for 5/1860  frac change 
% dq_CNG1_nominal  = (tx2p(:,jCNG1) - t0)*normer.normCNG/(log(1+normer.dQ));
% dq_CNG2_nominal  = (tx2p(:,jCNG2) - t0)*normer.normCNG/(log(1+normer.dQ));
%
% see ../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/put_together_sarta_jacs_anom40.m
% str1 = {'CO2'  'N2O'  'CH4'  'CFC11'  'ST'  'WV(97)'  'T(97)' 'O3(97)'};
% str2 = '[ [2.2 1.0 5 1 0.1] [0.01*ones(1,97)] [0.01*ones(1,97)] [0.01*ones(1,97)]]';
% qrenorm = [[2.2 1.0 5 1 0.1] [0.01*ones(1,97)] [0.01*ones(1,97)] [0.01*ones(1,97)]];

% %% NEW, see normer.cld
% str1 = {'CO2' 'N2O'  'CH4' 'CFC11'  'CFC12' 'ST'  'CNG1' 'CNG2' 'CSZ1' 'CSZ2' 'CPR1' 'CPR2' 'WV(97)'  'T(97)'  ' O3(97)'};
% str2 = '[ [2.2 1.0 5 1 1 0.1][0.001 0.001 0.001 0.001 0.001 0.001] [0.01*ones(1,97)] [0.01*ones(1,97)]] [0.01*ones(1,97)]]';
% qrenorm = [[2.2 1.0 5 1 1 0.1] [0.01*ones(1,6)] [0.01*ones(1,97)] [0.01*ones(1,97)] [0.01*ones(1,97)]];

%qrenorm(1) = 2.2/subjac.ppmv2(iLonBin);         %% 2.2 ppmv/400 ppmv
%qrenorm(2) = 1.0/(subjac.ppmv4(iLonBin)*1000);  %% 1 ppb/300 ppb
%qrenorm(3) = 5.0/(subjac.ppmv6(iLonBin)*1000);  %% 5 ppm/1840 ppb
qrenorm(1) = 2.2/kcarta.subjac.ppmv2;         %% 2.2 ppmv/400 ppmv
qrenorm(2) = 1.0/(kcarta.subjac.ppmv4*1000);  %% 1 ppb/300 ppb
qrenorm(3) = 5.0/(kcarta.subjac.ppmv6*1000);  %% 5 ppm/1840 ppb
qrenorm(4) = 1/300;  %% see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/driver_put_together_kcarta_jacs.m
                     %% https://agage.mit.edu/data/agage-data    CFC11 = 242 ppt
qrenorm(5) = 1/600;  %% see /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/driver_put_together_kcarta_jacs.m
                     %% https://agage.mit.edu/data/agage-data    CFC11 = 242 ppt
qrenorm(6) = 0.1;

iVersQRenorm = 1; %% DEFAULT
iVersQRenorm = 3;
iVersQRenorm = topts.iVersQRenorm;

%% see iLatX in build_cov_matrices.m
iLatX = topts.iLatX;
if iLatBin <= iLatX | iLatBin >= 64-iLatX
  iVersQRenorm = 4; %% TESTING POALR REGIONS JAN 2023
end

fprintf(1,'  ....  get_jac_fast.m : iVersQRenorm = %2i (default = 1) ....  \n',iVersQRenorm)

if iVersQRenorm == 1
  %% gives good T(z) rates, small WV rates, till AUg 2021
  qrenorm((1:nlays)+6+nlays*0) = 0.01; %% WV
  qrenorm((1:nlays)+6+nlays*1) = 0.01; %% T
  qrenorm((1:nlays)+6+nlays*2) = 0.01; %% O3
elseif iVersQRenorm == 2
  %% try this Aug 26, 2021
  qrenorm((1:nlays)+6+nlays*0) = 0.1; %% WV
  qrenorm((1:nlays)+6+nlays*1) = 0.01; %% T
  qrenorm((1:nlays)+6+nlays*2) = 0.1; %% O3
elseif iVersQRenorm == 3
  qrenorm((1:nlays)+6+nlays*0) = 0.01; %% WV
  qrenorm((1:nlays)+6+nlays*1) = 0.10; %% T
  qrenorm((1:nlays)+6+nlays*2) = 0.01; %% O3
elseif iVersQRenorm == 4
  qrenorm((1:nlays)+6+nlays*0) = 0.1; %% WV
  qrenorm((1:nlays)+6+nlays*1) = 1.00; %% T
  qrenorm((1:nlays)+6+nlays*2) = 0.1; %% O3

  qrenorm((1:nlays)+6+nlays*0) = 0.01; %% WV
  qrenorm((1:nlays)+6+nlays*1) = 0.10; %% T
  qrenorm((1:nlays)+6+nlays*2) = 0.01; %% O3

  qrenorm((1:nlays)+6+nlays*1) = 0.05; %% T
  qrenorm((1:nlays)+6+nlays*1) = 0.50; %% T
  qrenorm((1:nlays)+6+nlays*1) = 1.00; %% T
  qrenorm(6) = 1.0;

%  qrenorm((1:nlays)+6+nlays*1) = 0.10; %% T
%  qrenorm(6) = 0.1;
end

qrenormUSE = qrenorm;
m_ts_jac = (ones(2645,1) * qrenorm) .* m_ts_jac_fast;  %%% SO THIS IS NORMALIZED JAC!!!!!!!!!!!!!!!!! JTRUE * RENORM

%% but really
qrenorm(1) = 2.2;         %% 2.2 ppmv/400 ppm
qrenorm(2) = 1.0;         %% 1 ppb/300 ppb
qrenorm(3) = 5.0;         %% 5 ppm/1840 ppb
qrenorm(4) = 1.0;         %% 1 ppt/300 ppt
qrenorm(5) = 1.0;         %% 1 ppt/600 ppt
qrenormACTUAL = qrenorm;

iSaveJac = -1;
if iSaveJac > 0
  commentJacX = 'see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/get_jac_fast.m';
  if iVersJac == 2019
    saverjac = ['save /asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_latbin_' num2str(iLatBin,'%02d') '_lonbin_' num2str(iLonBin,'%02d') '.mat commentJacX qrenormUSE qrenormACTUAL nlays m_ts_jac'];
  elseif iVersJac == 2021
    saverjac = ['save /asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_ERA5_2021_latbin_' num2str(iLatBin,'%02d') '_lonbin_' num2str(iLonBin,'%02d') '.mat commentJacX qrenormUSE qrenormACTUAL nlays m_ts_jac'];
  elseif iVersJac == 2014
    saverjac = ['save /asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_ERA5_2014_latbin_' num2str(iLatBin,'%02d') '_lonbin_' num2str(iLonBin,'%02d') '.mat commentJacX qrenormUSE qrenormACTUAL nlays m_ts_jac'];
  end
  eval(saverjac);
end

oo = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g2_jac.mat');
plot(h.vchan,m_ts_jac(:,1),oo.fout,sum(oo.jout')*2.2/370)

