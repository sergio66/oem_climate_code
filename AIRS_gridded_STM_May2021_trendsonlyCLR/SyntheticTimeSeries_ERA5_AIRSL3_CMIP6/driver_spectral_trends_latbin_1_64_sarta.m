addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

system_slurm_stats
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 1

frp = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/simulate64binsERA5_' num2str(JOB) '.rp.rtp'];
[h,ha,pppp,pa] = rtpread(frp);

for xx = 1 : 72
  ind = (1:72:16416);
  ind = ind + (xx-1);
  boo = find(ind <= 19*12*72);
  ind = ind(boo); 
  plot(pppp.rlon(ind),'o-'); title(num2str(xx)); pause(0.1);

  [yy,mm,dd,hh] = tai2utcSergio(pppp.rtime(ind));
  dayOFtime = change2days(yy,mm,dd,2002);

  raaData = rad2bt(h.vchan,pppp.rcalc(:,ind));
  for iii = 1 : 2645
    if mod(iii,1000) == 0
      fprintf(1,'+');
    elseif mod(iii,100) == 0
      fprintf(1,'.');
    end
    data = double(raaData(iii,:));
    zoo = find(isfinite(data));
    if length(zoo) > 20
      [junk err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      thesave.xtrend(iii,xx)    = junk(2);
      thesave.xtrendErr(iii,xx) = err.se(2);
      thesave.const(iii,xx)     = junk(1);
    else
      thesave.xtrend(iii,xx)    = NaN;
      thesave.xtrendErr(iii,xx) = NaN;
      thesave.const(iii,xx)     = NaN;
    end
  end
end
dirout = ['ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
if ~exist(dirout)
  mker = ['!mkdir ' dirout];
  eval(mker);
end
comment = 'see MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_spectral_trends_latbin_1_64_sarta.m';
saver = ['save ' dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d') '.mat thesave comment'];
eval(saver);
