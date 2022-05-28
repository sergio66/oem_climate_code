addpath /home/sergio/MATLABCODE/CRIS_Hi2Lo/
i1305 = load('/asl/matlib/cris/ch_std_from1317.mat');

for i16daytimestep = 1 : 157
  newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS_FiniteDiff_Try3/'];       %% no seasonal, redone
  newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_5_2235_V4.mat']; %% Aug 21,         better CO2/CH4/N2O/CFC11/CFC12 prof, fixed 2002/09 yay
  junk = load(newjacname);
  co2cris(:,i16daytimestep) = squeeze(junk.tracegas(20,:,1));
end
[foutC,co2jacoutC] = translate_hi2lo(junk.f,co2cris);
foutC.ichan = foutC.ichan(i1305.ch_std_i);
foutC.vchan = foutC.vchan(i1305.ch_std_i);
co2jacoutC = co2jacoutC(i1305.ch_std_i,:); grid
plot(C.okdates,co2jacoutC(iC791-1:iC791,:)); grid

for i16daytimestep = 1 : 157
  if mod(i16daytimestep,25) == 0
    fprintf(1,'.');
  end
  junk = num2str(i16daytimestep,'%03d');
  newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2235.mat']; 
  junk = load(newjacname);
  stempcris(:,i16daytimestep) = squeeze(junk.M_TS_jac_all(20,:,6));
  co2cris(:,i16daytimestep) = squeeze(junk.M_TS_jac_all(20,:,1));
end
fprintf(1,'\n');

[foutC,stempjacoutC] = translate_hi2lo(junk.f,stempcris);
[foutC,co2jacoutC]   = translate_hi2lo(junk.f,co2cris);
foutC.ichan = foutC.ichan(i1305.ch_std_i);
foutC.vchan = foutC.vchan(i1305.ch_std_i);

stempjacoutC = stempjacoutC(i1305.ch_std_i,:); grid
plot(C.okdates,stempjacoutC(iC1231,:)/0.1); grid    %% 0.1 is the stemp normalization

co2jacoutC = co2jacoutC(i1305.ch_std_i,:); grid
plot(C.okdates,co2jacoutC(iC791,:)/0.1-co2jacoutC(iC792,:)/0.1); grid    %% 0.1 is the co2 normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i16daytimestep = 1 : 365
  if mod(i16daytimestep,25) == 0
    fprintf(1,'.');
  end
  junk = num2str(i16daytimestep,'%03d');
  newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];
  junk = load(newjacname);
  stemp0airs(:,i16daytimestep) = squeeze(junk.M_TS_jac_all(20,:,6));
  co2airs(:,i16daytimestep) = squeeze(junk.M_TS_jac_all(20,:,1));
end
fprintf(1,'\n');

figure(14); plot(C.okdates,stempjacoutC(iC1231,:)/0.1,'b',A.okdates,stempairs(iA1231,:)/0.1,'r'); grid    %% 0.1 is the stemp normalization
  title('stemp jacobians using 1231 cm-1');
  hl = legend('cris','airs','location','best');

figure(15); plot(iC.okdates,co2jacoutC(iC791,:)/0.1-co2jacoutC(iC792,:)/0.1,'b',...
                 iA.okdates,co2jacoutA(iA791,:)/0.1-co2jacoutA(iA792,:)/0.1,'r'); grid    %% 0.1 is the co2 normalization
  title('co2 jacobians using 1231 cm-1');
  hl = legend('cris','airs','location','best');

keyboard_nowindow
