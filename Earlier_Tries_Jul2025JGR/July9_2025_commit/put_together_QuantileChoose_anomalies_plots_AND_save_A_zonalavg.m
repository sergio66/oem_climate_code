a1 = load(fname,'yy_desc','mm_desc','dd_desc','hh_desc','rtime_desc');
yyD = a1.yy_desc;
mmD = a1.mm_desc;
ddD = a1.dd_desc;
hhD = a1.hh_desc;
rtimeD = a1.rtime_desc;

a1 = load(fname,'yy_asc','mm_asc','dd_asc','hh_asc','rtime_asc');
yyA = a1.yy_asc;
mmA = a1.mm_asc;
ddA = a1.dd_asc;
hhA = a1.hh_asc;
rtimeA = a1.rtime_asc;

coslat = cos(rlat*pi/180)*ones(1,iNumAnomTimeSteps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iSave = input('save (-1/+1) : ');
iSave = +1;
if iSave > 0
  comment = 'look at driver_put_together_QuantileChoose_anomalies_zonalavg.m';
  f_ind = h.vchan(ind);
  if iHuge < 0
    fout = ['anomaly_zonalavg_chID_' num2str(iChID,'%04d') '_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat'];  %%orig
    fout = ['anomaly_iQAX_' num2str(iQAX) '_zonalavg_chID_' num2str(iChID,'%04d') '_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '.mat'];
    saver = ['save ' fout ' bt1231_D btChID_D btanomD yyD mmD ddD hhD rtimeD RRTM_bands fluxanomD    btanomA yyA mmA ddA hhA rtimeA bt1231_A btChIDn_A  fluxanomA  comment ind f_ind'];
  else
    fout = ['/asl/s1/sergio/JUNK/anomaly_zonalavg_ALL_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '_A.mat'];  %% orig
    fout = ['/asl/s1/sergio/JUNK/anomaly_iQAX_' num2str(iQAX) '_zonalavg_ALL_Q' num2str(iQuant,'%02d') '_numyears_' num2str(iNumYears,'%6.2f') '_iNumAnomTimeSteps_' num2str(iNumAnomTimeSteps) '_A.mat'];
    saver = ['save -v7.3 ' fout ' btanomA yyA mmA ddA hhA rtimeA comment bt1231_A RRTM_bands fluxanomA ind f_ind'];
  end
  eval(saver)
end

