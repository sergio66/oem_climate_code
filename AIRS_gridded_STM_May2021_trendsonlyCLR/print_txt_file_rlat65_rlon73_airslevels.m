%% For Stephen Leroy
clear all
do_XX_YY_from_X_Y
fid = fopen('rlat65.txt','w'); fprintf(fid,'%8.6f \n',rlat65); fclose(fid)
fid = fopen('rlon73.txt','w'); fprintf(fid,'%8.6f \n',rlon73); fclose(fid)
boo = load('/home/sergio/MATLABCODE/airslevels.dat');
fid = fopen('airslevels.txt','w'); fprintf(fid,'%8.6f \n',boo); fclose(fid)
