addpath /home/sergio/MATLABCODE

load climateQA_LW_noCFC.mat
load ../f2645.mat

fL1c = f2645(chanset);

fairs = instr_chans;

for ii = 1 : length(fL1c)
  boo = abs(fL1c(ii)-fairs);
  chansetL1b(ii) = find(boo == min(boo),1);
end

array = [fL1c    chanset  chansetL1b'];

fid = fopen('climateQA_LW_noCFC.txt','w');
fprintf(fid,'L1c center freq       L1c chanID      L1b chanID \n');
fprintf(fid,'-------------------------------------------------\n');
fprintf(fid,'  %08.3f             %4i             %4i \n',array');
fclose(fid);
