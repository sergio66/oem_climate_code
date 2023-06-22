dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/UMBC/';

iDo = +1;
for ii = 1 : 64
  fin  = [dir0 'simulate64binsERA5_' num2str(ii) '_2002_09_2022_08.ip.rtp'];
  fout = [dir0 'simulate64binsUMBC_' num2str(ii) '_2002_09_2022_08.ip.rtp'];
  mver = ['!/bin/mv ' fin '  ' fout];
  eval(mver);

  fin  = [dir0 'simulate64binsERA5_' num2str(ii) '_2002_09_2022_08.op.rtp'];
  fout = [dir0 'simulate64binsUMBC_' num2str(ii) '_2002_09_2022_08.op.rtp'];
  mver = ['!/bin/mv ' fin '  ' fout];
  eval(mver);

  fin  = [dir0 'simulate64binsERA5_' num2str(ii) '_2002_09_2022_08.rp.rtp'];
  fout = [dir0 'simulate64binsUMBC_' num2str(ii) '_2002_09_2022_08.rp.rtp'];
  mver = ['!/bin/mv ' fin '  ' fout];
  eval(mver);

end

