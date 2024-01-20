addpath /asl/matlib/science/
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

[h5,ha,p5,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp');
[p5.salti,p5.landfrac] = usgs_deg10_dem(p5.rlon,p5.rlat);

iProf = 10; %% midlats
iProf = 02; %% polar
iProf = 20; %% default, tropics

[hnew,pnew] = replicate_rtp_headprof(h5,p5,iProf,5);  %% raw, increase WV, CO2, T, all three
pnew.upwell = 2 * ones(size(pnew.stemp));
pnew.scanang = 0 * pnew.scanang;
pnew.satzen  = 0 * pnew.satzen;
pnew.solzen  = 150 * ones(size(pnew.satzen));

pnew.gas_1(:,2) = pnew.gas_1(:,2) * 1.001;
pnew.gas_2(:,3) = pnew.gas_2(:,3) * (1 + 2.2/380);
pnew.ptemp(:,4) = pnew.ptemp(:,4) + 0.02;

pnew.gas_1(:,5) = pnew.gas_1(:,5) * 1.001;
pnew.gas_2(:,5) = pnew.gas_2(:,5) * (1 + 2.2/380);
pnew.ptemp(:,5) = pnew.ptemp(:,5) + 0.02;

rtpwrite('upwell.op.rtp',hnew,ha,pnew,pa);

pdn = pnew;
pdn.upwell = ones(size(pdn.upwell));
rtpwrite('dnwell.op.rtp',hnew,ha,pdn,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% keep updating iRTP and then
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up_jacZ.nml   UPLOOKJAC/uplook1.dat UPLOOKJAC/uplook1.jac']; eval(kcer)
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up_coljac.nml UPLOOKJAC/uplook1.dat UPLOOKJAC/uplook1.jac']; eval(kcer)

kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up.nml UPLOOKJAC/uplook1.dat']; eval(kcer)
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up.nml UPLOOKJAC/uplook2.dat']; eval(kcer)
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up.nml UPLOOKJAC/uplook3.dat']; eval(kcer)
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up.nml UPLOOKJAC/uplook4.dat']; eval(kcer)
kcer = ['!time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20 quickuse_up.nml UPLOOKJAC/uplook5.dat']; eval(kcer)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spacing = 2; width = 2;

addpath /home/sergio/KCARTA/MATLAB
[dradall(:,1),wall] = readkcstd('uplook1.dat');
if ~exist('djacZ')
  [djacZ,wall] = readkcjac('uplook1.jac');
  [fc,jacc] = quickconvolve(wall,djacZ,spacing,width);
end

%% order is Gas1(z), Gas2(z), Gas101(z), Gas102(z), Gas103(z), T(z), WgtFcn(z)
[sum(sum(jacc(:,(1:97)+97*0),2)*0.001) sum(sum(jacc(:,(1:97)+97*1),2)*2.2/380) sum(sum(jacc(:,(1:97)+97*2),2)*0.02)]
jaccSum(:,1) = (sum(jacc(:,(1:97)+0*97)+jacc(:,(1:97)+2*97) + jacc(:,(1:97)+3*97) + jacc(:,(1:97)+4*97),2))*0.001;
jaccSum(:,2) = (sum(jacc(:,(1:97)+1*97),2))*2.2/380;
jaccSum(:,3) = (sum(jacc(:,(1:97)+5*97),2))*0.02;
figure(1); clf; plot(fc,jaccSum,'linewidth',2); 
  xlabel('Wavenumber cm-1'); ylabel('\delta radiance(mW/m2/sr-1/cm-1)');   
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','location','best'); title('kCARTA 97 layer jacs, summed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xdcoljac,wall] = readkcBasic('uplook1.jac_COL');
Tdcoljac = rad2bt(wall,xdcoljac) - rad2bt(wall,dradall(:,1)*ones(1,7));
xdcoljac = xdcoljac - dradall(:,1)*ones(1,7);

iLinOrLog = +1;
iLinOrLog = -1;

  % Revisions.txt
  % 6/25/19        Changed the column jac mult    from 0.1 to 0.001
  %                Changed the column temp offset from 1.0 to 0.01
  %                column jac perturbations : see kDefaultToffset and kDefaultColMult in INLUDE/*.param
dT = 0.01;
dQ = 0.001;
dcoljac = zeros(890000,3);
dcoljac(:,1) = xdcoljac(:,1) + xdcoljac(:,3) + xdcoljac(:,4) + xdcoljac(:,5);
dcoljac(:,2) = xdcoljac(:,2);
dcoljac(:,3) = xdcoljac(:,6);
if iLinOrLog == -1
  dcoljac(:,1:2) = dcoljac(:,1:2)/log(1+dQ);
else
  dcoljac(:,1:2) = dcoljac(:,1:2)/(dQ);
end
dcoljac(:,3)   = dcoljac(:,3)/dT;

tcoljac = zeros(890000,3);
tcoljac(:,1) = Tdcoljac(:,1) + Tdcoljac(:,3) + Tdcoljac(:,4) + Tdcoljac(:,5);
tcoljac(:,2) = Tdcoljac(:,2);
tcoljac(:,3) = Tdcoljac(:,6);
if iLinOrLog == -1
  tcoljac(:,1:2) = tcoljac(:,1:2)/log(1+dQ);
else
  tcoljac(:,1:2) = tcoljac(:,1:2)/(dQ);
end
tcoljac(:,3)   = tcoljac(:,3)/dT;

figure(2); clf; plot(wall,dcoljac);
figure(2); clf; plot(wall,tcoljac);

dcoljacX(:,1) = dradall(:,1);
if iLinOrLog == -1
  dcoljacX(:,2) = tcoljac(:,1)*log(1.001); 
  dcoljacX(:,3) = tcoljac(:,2)*log(1+2.2/380);
else
  dcoljacX(:,2) = tcoljac(:,1)*(0.001); 
  dcoljacX(:,3) = tcoljac(:,2)*(0+2.2/380);
end
dcoljacX(:,4) = tcoljac(:,3)*0.02;    
[fc,qcTcol] = quickconvolve(wall,dcoljacX,spacing,width);
%figure(2); clf; plot(fc,rad2bt(fc,qc(:,2:4))-rad2bt(fc,qc(:,1))*ones(1,3),'linewidth',2); xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)'); 
figure(2); clf; plot(fc,qcTcol(:,2:4),'linewidth',2); xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)'); title('kCARTA col jacs')
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','location','best');

dcoljacX(:,1) = dradall(:,1);
if iLinOrLog == -1
  dcoljacX(:,2) = dcoljac(:,1)*log(1.001); 
  dcoljacX(:,3) = dcoljac(:,2)*log(1+2.2/380);
else
  dcoljacX(:,2) = dcoljac(:,1)*(0.001); 
  dcoljacX(:,3) = dcoljac(:,2)*(0+2.2/380);
end
dcoljacX(:,4) = dcoljac(:,3)*0.02;    
[fc,qcRcol] = quickconvolve(wall,dcoljacX,spacing,width);
figure(3); clf; plot(fc,qcRcol(:,2:4),'linewidth',2); xlabel('Wavenumber cm-1'); ylabel('\delta radiance(mW/m2/sr-1/cm-1)');  title('kCARTA col jacs')
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','location','best');

dcoljacflux = sum(dcoljacX(:,2:4),1)*0.0025*pi/1000;
dcoljacflux(4) = sum(dcoljacflux);
junk = [dcoljacflux; dcoljacflux/dcoljacflux(4)];
disp('    WV   CO2   T   Total')
fprintf(1,'%8.2f %8.2f %8.2f %8.2f \n',junk');
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%

[dradall(:,2),wall] = readkcstd('uplook2.dat');
[dradall(:,3),wall] = readkcstd('uplook3.dat');
[dradall(:,4),wall] = readkcstd('uplook4.dat');
[dradall(:,5),wall] = readkcstd('uplook5.dat');

figure(4); clf; plot(wall,rad2bt(wall,dradall(:,2:5))-rad2bt(wall,dradall(:,1))*ones(1,4)); xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)');  title('finite changes kCARTA')
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','all three','location','best');

raddiff = sum(dradall(:,2:5) - dradall(:,1)*ones(1,4),1)*2*pi/2*0.0025/1000;
junk = [raddiff; raddiff/raddiff(4)];
disp('    WV   CO2   T   Total')
fprintf(1,'%8.2f %8.2f %8.2f %8.2f \n',junk');
disp(' ')

figure(5); clf; plot(wall,dradall(:,2:5) - dradall(:,1)*ones(1,4)); xlabel('Wavenumber cm-1'); ylabel('\delta radiance(mW/m2/sr-1/cm-1)'); title('finite changes kCARTA')
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','all three','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%
[fc,qc] = quickconvolve(wall,dradall,spacing,width);
figure(4); clf; plot(fc,rad2bt(fc,qc(:,2:5))-rad2bt(fc,qc(:,1))*ones(1,4),'linewidth',2); xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)'); title('finite changes kCARTA')
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','all three','location','best');

figure(5); clf; plot(fc,qc(:,2:5) - qc(:,1)*ones(1,4),'linewidth',2); xlabel('Wavenumber cm-1'); ylabel('\delta radiance(mW/m2/sr-1/cm-1)'); title('finite changes kCARTA')
  xlim([605 1640])
  hl = legend('WV x1.001','CO2 x 1.0059','T+0.02','all three','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); plot(fc,qcRcol(:,2:4),fc,qc(:,2:4) - qc(:,1)*ones(1,3))
figure(6); plot(fc,qcRcol(:,3),'b.-',fc,qc(:,3)-qc(:,1),'gx-',fc,jaccSum(:,2),'r','linewidth',2); xlim([600 1600]); legend('coljac','sum(layer) jac','finite diff','location','best')
figure(6); plot(fc,qcRcol(:,4),'b.-',fc,qc(:,4)-qc(:,1),'gx-',fc,jaccSum(:,3),'r','linewidth',2); xlim([600 1600]); legend('coljac','sum(layer) jac','finite diff','location','best')
figure(6); plot(fc,qcRcol(:,2),'b.-',fc,qc(:,2)-qc(:,1),'gx-',fc,jaccSum(:,1),'r','linewidth',2); xlim([600 1600]); legend('coljac','sum(layer) jac','finite diff','location','best')

