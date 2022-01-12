addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat
[h,ha,p,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/simulate64binsERA5_32.rp.rtp');

for JOB = 1:72;
  fprintf(1,'JOB = %2i \n',JOB)

  ind = (1:72:16416);
  ind = (ind) + (JOB-1);
  
  [yy,mm,dd,hh] = tai2utcSergio(p.rtime(ind));
  dayOFtime = change2days(yy,mm,dd,2002);
  
  loader = ['load /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/ERA5_LATBIN32_72lonbinsx228timesteps/timeseries_latbin32_lonbin' num2str(JOB,'%02d') '.mat'];
  eval(loader);
  fKc = fKc(ichan);
  raaData = raaData(ichan,:);
  raaData = rad2bt(fKc,raaData);
  
  warning off
  for iii = 1 : 2645
    if mod(iii,1000) == 0
      fprintf(1,'+');
    elseif mod(iii,100) == 0
      fprintf(1,'.');
    end
    data = double(raaData(iii,:));
    zoo = find(isfinite(data));
    if length(zoo) > 20
      junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      thesave.xtrendSpectral(iii) = junk(2);
      xconstr72(iii) = junk(1);
    else
      thesave.xtrendSpectral(iii) = NaN;
      xconstr72(iii) = NaN;
    end
  end
  
  kctrend = thesave.xtrendSpectral;
  saver = ['save KCARTA_latbin32_spectral_trends/kcarta_spectraltrends_latbin32_lonbin' num2str(JOB,'%02d') '.mat fKc kctrend'];
  eval(saver);
  plot(fKc,thesave.xtrendSpectral); title(num2str(JOB)); pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now compare to sarta

main_compare_sarta_kcarta_latbin32_trends
