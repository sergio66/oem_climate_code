function [m_ts_jac,nlays,qrenorm] = get_jac(fname,iWhichLonBin);

loader = ['load ' fname];
eval(loader);

%% no need to add in extra column for CFC12 >>>>>>
m_ts_jac = [];

boo = squeeze(subjac.jacCO2z(:,:,iWhichLonBin)); boo = sum(boo,1);  m_ts_jac = [m_ts_jac boo'];
boo = subjac.coljacN2O(:,iWhichLonBin);                             m_ts_jac = [m_ts_jac boo];
boo = subjac.coljacCH4(:,iWhichLonBin);                             m_ts_jac = [m_ts_jac boo];
boo = subjac.jacCld1(:,iWhichLonBin);                               m_ts_jac = [m_ts_jac boo];
boo = subjac.jacCld2(:,iWhichLonBin);                               m_ts_jac = [m_ts_jac boo];
boo = subjac.jacST(:,iWhichLonBin);                                 m_ts_jac = [m_ts_jac boo];

nlays = subjac.nlevs(iWhichLonBin)-1;
boo = squeeze(subjac.jacWV(1:nlays,:,iWhichLonBin));                 m_ts_jac = [m_ts_jac boo'];
boo = squeeze(subjac.jacT(1:nlays,:,iWhichLonBin));  boo(1:3,:) = 0; m_ts_jac = [m_ts_jac boo'];
boo = squeeze(subjac.jacO3(1:nlays,:,iWhichLonBin));                 m_ts_jac = [m_ts_jac boo'];

[mm,nn] = size(m_ts_jac);
qrenorm = ones(1,nn);

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

qrenorm(1) = 2.2/subjac.ppmv2(iWhichLonBin);         %% 2.2 ppmv/400 ppmv
qrenorm(2) = 1.0/(subjac.ppmv4(iWhichLonBin)*1000);  %% 1 ppb/300 ppb
qrenorm(3) = 5.0/(subjac.ppmv6(iWhichLonBin)*1000);  %% 5 ppm/1840 ppb
qrenorm(4) = 0.001;
qrenorm(5) = 0.001;
qrenorm(6) = 0.1;
qrenorm((1:nlays)+6+nlays*0) = 0.01;
qrenorm((1:nlays)+6+nlays*1) = 0.01;
qrenorm((1:nlays)+6+nlays*2) = 0.01;

m_ts_jac = (ones(2645,1) * qrenorm) .* m_ts_jac;

%% but really
qrenorm(1) = 2.2;         %% 2.2 ppmv/400 ppmv
qrenorm(2) = 1.0;         %% 1 ppb/300 ppb
qrenorm(3) = 5.0;         %% 5 ppm/1840 ppb

