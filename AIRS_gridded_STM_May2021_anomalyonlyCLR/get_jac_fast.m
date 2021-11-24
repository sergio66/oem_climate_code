function [m_ts_jac,nlays,qrenorm,freq2645] = get_jac(fname,iWhichLatLonBin,iLonBin,iLatBin);

iDoSlow = -1;
if iDoSlow > 0
 tic;
  this_is_slow_jac_load
  'done1'
  timervalx1 = toc;
end

if iDoSlow > 0
  tic;
end
see_clust_put_together_jacs_clr
if iDoSlow > 0
  'done2'
  timervalx2 = toc;

  timervalx1
  timervalx2

  [mm,nnS] = size(m_ts_jac_slow)
  [mm,nnF] = size(m_ts_jac_fast)

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
qrenorm((1:nlays)+6+nlays*0) = 0.01; %% WV
qrenorm((1:nlays)+6+nlays*1) = 0.01; %% T
%qrenorm((1:nlays)+6+nlays*1) = 0.10; %% T
qrenorm((1:nlays)+6+nlays*2) = 0.01; %% O3

m_ts_jac = (ones(2645,1) * qrenorm) .* m_ts_jac_fast;

%% but really
qrenorm(1) = 2.2;         %% 2.2 ppmv/400 ppm
qrenorm(2) = 1.0;         %% 1 ppb/300 ppb
qrenorm(3) = 5.0;         %% 5 ppm/1840 ppb
qrenorm(4) = 1.0;         %% 1 ppt/300 ppt
qrenorm(5) = 1.0;         %% 1 ppt/600 ppt

oo = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g2_jac.mat');
plot(h.vchan,m_ts_jac(:,1),oo.fout,sum(oo.jout')*2.2/370)
