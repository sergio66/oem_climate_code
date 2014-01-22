% do_retrieval_hand.m
%
% Script to do OEM and LLS retrieval of radiance rates
% Retrieval parameters set in get_struct
%
% V0.1 Sergio, Larrabee: Sept 11, 2012
%
% can easily do all latbins via
%     for JOB = 1 : 36
%       run_retrieval_cluster_AIRS
%     end
% 
% /bin/rm slurm*; /bin/rm ../Output/*.mat; clustcmd -q short -n 36 run_retrieval_latbins_AIRS.m 1:N

ix = JOB;

add_path_to_oem_pkg

%% Define default driver structure
driver = set_struct;

%% Open debug file if desired
if driver.debug
  writelog('open');
end;

%% Change some defaults
%driver = override_defaults_airs_test(driver,ix);      %% this is TEST and should always work!
driver = override_defaults_latbins_AIRS(driver,ix);  %% this is YOUR settings

% Get rate data and Jacobians
driver            = get_rates(driver);
[driver,m_ts_jac] = get_jacs(driver);          %% strow's new renorm
%[driver,m_ts_jac] = get_jacs_NOrenorm(driver); %% no renorm


%  Adjust the rates?
if driver.rateset.adjust
  % m_ts_jac0 = get_jacs_6_97_97(driver);
  m_ts_jac0 = get_jacs0(driver);
  driver = adjust_rates(driver,m_ts_jac0);
end
figure(1)
% Do the retrieval
driver = retrieval(driver,m_ts_jac);

%% Save retrieval output
save(driver.filename,'-struct','driver');

%% Close debug file
if driver.debug
  writelog('close')
end

profile=['/home/tangborn/Matlab/oem_pkg_run/Tangborn/prof',num2str(ix),'.dat']
fid=fopen(profile,'w');
fprintf(fid,'%12.4e\n',driver.oem.finalrates); 
profsigs=['/home/tangborn/Matlab/oem_pkg_run/Tangborn/sigs',num2str(ix),'.dat']
fid2=fopen(profsigs,'w');
fprintf(fid2,'%12.4e\n',driver.oem.finalsigs); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa = load(driver.rateset.datafile);

figure(1);
  g  = dogoodchan; ff = instr_chans;
  g1 = driver.jacobian.chanset_used;
  plot(ff(g1),driver.rateset.rates(g1),ff(g1),squeeze(aa.b_obs(ix,g1,2)),'k',ff(g1),driver.oem.fit(g1),'r.-'); grid; 
    axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','DEFINITELY OBS','fits'); set(hl,'fontsize',10)
figure(2);
  g  = dogoodchan; ff = instr_chans;
  g1 = driver.jacobian.chanset_used;
  plot(ff(g1),driver.rateset.rates(g1),ff(g1),driver.rateset.rates(g1)-driver.oem.fit(g1)','r'); grid; axis([500 3000 -0.15 +0.15])
  title('AIRS'); hl=legend('data','fits'); set(hl,'fontsize',10)
if length(driver.oem.finalrates) == 200
%  figure(3)
%  plot(driver.oem.finalrates(7:103),1:97,driver.oem.finalrates(104:200),1:97,'r')
%  set(gca,'ydir','reverse');
%  title('AIRS (b) : WV frac/yr (r) T K/yr'); grid
end

plot_retrieval_latbins

