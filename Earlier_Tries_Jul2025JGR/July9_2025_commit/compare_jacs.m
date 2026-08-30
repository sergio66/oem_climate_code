%% this is for rates, constant jac
jacobian.filename1 = '../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat';

iLatBin = 35;  %% at least in spectral shape they look like CO2, though too big in magnitude (esp for v2)
iLatBin = 16;  %% too big in magnitude and loo spectrally incorrect

i16daytimestep = input('enter timestep (1 : 365) : ');
junk = num2str(i16daytimestep,'%03d');

%% sarta time vary jacs
jacobian.filename2 = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16/RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2645.mat'];

%% kcarta time vary jacs, dQ = 0.1
jacobian.filename3 = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_dQpert0.1_dTpert1.0/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];

%% kcarta time vary jacs, dQ = 0.001, bad tracegas profiles
jacobian.filename4 = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_tillJuly01_2019/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now have to worry about jacobian directories; from here on kcarta finite gas jacobians used dQ = 0.001
JacDir = 'Anomaly365_16_12p8_July12_2019_Great_But_with_Seasonal';   %% this used the nice tracegas profiles (CO2/CH4/N2O) ... worked fine ... but the T WV O stemp had seasonal in them

%% kcarta time vary jacs, dQ = 0.001, better tracegas profiles
jacobian.filename5 = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/' JacDir '/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];

%% see replace_time_co2jac.m
newjacname1 = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/' JacDir '/RESULTS_FiniteDiff_Try1/'];
newjacname1 = [newjacname1 '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645.mat'];   %% June 30/July 1, bad CO2/CH4/N2O prof

%% see replace_time_co2jac.m
newjacname2 = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/' JacDir '/RESULTS_FiniteDiff_Try2/'];
newjacname2 = [newjacname2 '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645_V2.mat']; %% July 2,  better CO2/CH4/N2O prof, screwed at 2002/09

%% see replace_time_co2jac.m
newjacname3 = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/' JacDir '/RESULTS_FiniteDiff_Try3/'];
newjacname3 = [newjacname3 '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645_V3.mat']; %% July 12,  better CO2/CH4/N2O prof, correct at 2002/09

if ~exist('a1')
  a1 = load(jacobian.filename1);
end
if ~exist('a2')
  a2 = load(jacobian.filename2);
end
if ~exist('a3')
  a3 = load(jacobian.filename3);
end
if ~exist('a4')
  a4 = load(jacobian.filename4);
end
if ~exist('a5')
  a5 = load(jacobian.filename5);
end
if ~exist('a6')
  a6 = load(newjacname1);
  a6 = load(newjacname2);
  a6 = load(newjacname3);
end

figure(1); 
plot(a1.f,squeeze(a1.M_TS_jac_all(iLatBin,:,1)),'c',a1.f,squeeze(a2.M_TS_jac_all(iLatBin,:,1)),'m',a1.f,squeeze(a3.M_TS_jac_all(iLatBin,:,1)),'g',...
     a1.f,squeeze(a4.M_TS_jac_all(iLatBin,:,1)),'b.-',a1.f,squeeze(a5.M_TS_jac_all(iLatBin,:,1)),'r',a1.f,squeeze(a6.tracegas(iLatBin,:,1)),'k','linewidth',2);
hl = legend('kcarta constant','sarta time vary','kcarta timevary, dq=0.1 bad tracegas','kcarta timevary, dq=0.001 BAD tracegas','kcarta timevary, dq=0.001 GOOD tracegas','strow','location','best');
set(hl,'fontsize',8)
axis([640 840 -0.1 +0.05]); grid
