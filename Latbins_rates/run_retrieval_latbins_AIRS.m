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

% Do the retrieval
driver = retrieval(driver,m_ts_jac);

%% Save retrieval output
save(driver.filename,'-struct','driver');

%% Close debug file
if driver.debug
  writelog('close')
end

figure(1); 
  g=dogoodchan; ff = instr_chans;
  plot(ff(g),driver.rateset.rates(g),ff(g),driver.oem.fit(g),'r'); grid; axis([500 3000 -0.15 +0.15])
  title('AIRS');
if length(driver.oem.finalrates) == 200
  figure(2)
  plot(driver.oem.finalrates(7:103),1:97,driver.oem.finalrates(104:200),1:97,'r')
  set(gca,'ydir','reverse');
  title('AIRS (b) : WV frac/yr (r) T K/yr'); grid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
end
if strcmp(user_name(1:end-1),'sergio')
  ecmfile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/';
  ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012.mat'];
  load(ecmfile);

  if length(driver.oem.finalrates) == 200
    water = driver.oem.finalrates(7:103);
    watersigs = driver.oem.finalsigs(7:103);
    temp = driver.oem.finalrates(104:200);
    tempsigs = driver.oem.finalsigs(104:200);
 elseif length(driver.oem.finalrates) == 103 & length(driver.jacobian.Q1jacindex) < 1
    temp = driver.oem.finalrates(7:103);
    tempsigs = driver.oem.finalsigs(7:103);
    water = temp * 0;
    watersigs = temp*0;
 elseif length(driver.oem.finalrates) > 103 & length(driver.jacobian.Q1jacindex) < 97
    bonk = length(driver.jacobian.Q1jacindex);
    temp = driver.oem.finalrates((7:103)+bonk);
    tempsigs = driver.oem.finalsigs((7:103)+bonk);
    water = temp * 0;      water(driver.jacobian.Q1jacindex) = driver.oem.finalrates(7:7+length(driver.jacobian.Q1jacindex)-1);
    watersigs = temp * 0;  watersigs(driver.jacobian.Q1jacindex) = driver.oem.finalsigs(7:7+length(driver.jacobian.Q1jacindex)-1);
  else
    return
  end
  figure(6); clf
  subplot(121)
  shadedErrorBarY(water,1:97,watersigs,'bo-',1);
  hold on
  shadedErrorBarY(waterrate(ix,:),1:97,waterratestd(ix,:),'rx-',1);
  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 100]); 

  subplot(122)
  shadedErrorBarY(temp,1:97,tempsigs,'bo-',1);
  hold on
  shadedErrorBarY(ptemprate(ix,:),1:97,ptempratestd(ix,:),'rx-',1);
  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 100]);
  
end
