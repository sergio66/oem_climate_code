%% this is plotting

ecmfile = '/asl/s1/rates/clear/Aug2013/';
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

plevs = load('airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

%% shadedErrorBarYLog(temp,plays,tempsigs*1e4,'bo-',1);

%figure(2); clf
%  subplot(121)
%  shadedErrorBarYLog(water,plays,watersigs,'bo-',1);
%  hold on
%  shadedErrorBarYLog(waterrate(ix,:),plays,waterratestd(ix,:),'rx-',1);
%  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
%  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

%  subplot(122)
%  shadedErrorBarYLog(temp,plays,tempsigs,'bo-',1);
%  hold on
%  shadedErrorBarYLog(ptemprate(ix,:),plays,ptempratestd(ix,:),'rx-',1);
%  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10); 
%  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);

figure(5); clf
'test 7' 
  subplot(121)
  shadedErrorBarYLog10(water,plays,watersigs,'bo-');
  hold on
  shadedErrorBarYLog10(waterrate(ix,:),plays,waterratestd(ix,:),'rx-');
  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 1 3]); 
  rms_wtr_strat=rms(water(1:49)'-waterrate(ix,1:49))./rms(waterrate(ix,1:49)) 
  rms_wtr_trop=rms(water(50:97)'-waterrate(ix,50:97))./rms(waterrate(ix,50:97)) 
  subplot(122)
  shadedErrorBarYLog10(temp,plays,tempsigs,'bo-');
  hold on
  shadedErrorBarYLog10(ptemprate(ix,:),plays,ptempratestd(ix,:),'rx-');
  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.2 +0.15 1 3]);
  rms_tmp_strat=rms(temp(1:49)'-ptemprate(ix,1:49))./rms(ptemprate(ix,1:49)) 
  rms_tmp_trop=rms(temp(50:97)'-ptemprate(ix,50:97))./rms(ptemprate(ix,50:97))
'test 9' 
%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
end

if strcmp(user_name(1:end-1),'sergio')
  hz = p2h(plays)/1000; %% change to km
  firstderW = gradient(water,hz); secondderW = gradient(firstderW,hz);
  firstderT = gradient(temp,hz);  secondderT = gradient(firstderT,hz);
  p50  = find(plays >= 50);
  p100 = find(plays >= 100);

  W1x             = zerocross(firstderW);   T1x             = zerocross(firstderT); 
  [MaxW1x,MinW1x] = localmaxmin(firstderW); [MaxT1x,MinT1x] = localmaxmin(firstderT);

  second_deriv = sqrt([nanmean(secondderW(p50).*secondderW(p50)) nanmean(secondderT(p50).*secondderT(p50))]);
  fprintf(1,'  ---->>>> 50mb-1000mb average second deriv (squared) of W, T = %8.6f %8.6f \n',second_deriv);

  % figure(3); clf; plot(secondderW,1:97,secondderT,1:97,'r')
  % figure(3); clf; semilogy(secondderW,plays,secondderT,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  % figure(3); clf; semilogy(firstderW,plays,secondderW,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  % figure(3); clf; semilogy(firstderT,plays,secondderT,plays,'r')
  % set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 1 1000]);

  figure(6); clf; 
  subplot(121); semilogy(temp,plays)
    set(gca,'ydir','reverse'); grid; axis([-0.1 +0.1 20 1000]);
    hl = legend('T(z)','Location','northeast'); set(hl,'Fontsize',8);
  subplot(122); semilogy(firstderT,plays,...
                         secondderT,plays,'r',...
                         zeros(size(T1x)),plays(ceil(T1x)),'bo',...
                         firstderT(MaxT1x),plays(MaxT1x),'bd',...
                         firstderT(MinT1x),plays(MinT1x),'bs')

    set(gca,'ydir','reverse'); grid; axis([-0.1 +0.1 20 1000]);
    hl = legend('dT/dz','d2T/dz2','Location','northeast'); set(hl,'Fontsize',8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(4); clf; imagesc(log10(abs(log10(abs(driver.oem.inv_se))))); caxis([-1 2]); colorbar

