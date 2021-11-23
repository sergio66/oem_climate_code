function m_ts_jac_coljac = replace_time_co2_3layjac(m_ts_jac_coljac,iibin,i16daytimestep);

if i16daytimestep < 0
  %% disp('use constant jac')
else
  fprintf(1,'replacing co2 jac with that at timestep %4i of 365 \n',i16daytimestep)
  lstrow = load('sarta_chans_for_l1c.mat');
  newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CO2_370_385_400_415_12p8/co2_3layjac_2834_latbin' num2str(iibin) '.mat'];
  new = load(newjacname);
  anew_co2_jac_alltime = new.qclowerx(lstrow.ichan,:);
  bnew_co2_jac_alltime = new.qcmidx(lstrow.ichan,:);
  cnew_co2_jac_alltime = new.qcupperx(lstrow.ichan,:);
  new_iaCO2 = new.iaCO2;

  %% we know the nominal is 370 + (yy-2002)*2.2
  yyjunk = 2002:2030;
  co2junk = 370 + (yyjunk-2002)*2.2;
  new_timeCO2 = interp1(co2junk,yyjunk,new_iaCO2);       %% these are the jac times
%  plot(yyjunk,co2junk,'.-',new_timeCO2,new_iaCO2,'o-')
%  whos new_co2_jac_alltime new_timeCO2

%{
>> load anomaly_0dayavg_results_spectra_cal_constantCO2jac.mat
>> whos
  Name              Size               Bytes  Class     Attributes

  chanset         532x1                 4256  double
  iaTropics         1x20                 160  double
  okdates           1x365               2920  double
  okrtime           1x365               2920  double
  raaCal         2645x365            7723400  double
  raaObs         2645x365            7723400  double

>> save ok365times okdates okrtime
%}

  oktimes = load('ok365times.mat');
  oktimes = oktimes.okdates(i16daytimestep);   %% need to get the CO2 jac at this time!!!!

  for ii = 1 : 2645
    anew_co2_jac(ii) = interp1(new_timeCO2,anew_co2_jac_alltime(ii,:),oktimes);
    bnew_co2_jac(ii) = interp1(new_timeCO2,bnew_co2_jac_alltime(ii,:),oktimes);
    cnew_co2_jac(ii) = interp1(new_timeCO2,cnew_co2_jac_alltime(ii,:),oktimes);
  end

%  plot(1:2645,m_ts_jac_coljac(:,1),'k.-',1:2645,anew_co2_jac,'r',1:2645,anew_co2_jac_alltime(:,1),'b',1:2645,anew_co2_jac_alltime(:,10),'c')
%    hl = legend('orig col','new lower atm','lower atm370','lower atm405','location','best');

  new_m_ts_jac_coljac(:,1) = anew_co2_jac;
  new_m_ts_jac_coljac(:,2) = bnew_co2_jac;
  new_m_ts_jac_coljac(:,3) = cnew_co2_jac;
  new_m_ts_jac_coljac(:,4:7) = m_ts_jac_coljac(:,2:5);

  m_ts_jac_coljac = new_m_ts_jac_coljac;
end
