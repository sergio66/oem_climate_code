addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
iaTropics = find(abs(latbinsx) <= 30);
iaTropics = find(abs(latbinsx) <= 7.5);

iCnt = 0;
for yy = 2003 : 2017
  for doy = 1 : 366
    fname = ['/asl/rtp/rtp_airicrad_v6/clear/' num2str(yy) '/era_airicrad_day' num2str(doy,'%03d') '_clear.rtp'];
    if exist(fname)      
      iCnt = iCnt + 1;
      fprintf(1,'iCnt = %5i  %s \n',iCnt,fname)

      sergio.year(iCnt) = yy;
      sergio.day_of_year(iCnt) = doy;

      [h,ha,p,pa] = rtpread(fname);

      for jj = 1 : length(iaTropics)-1
        iaX = find(p.solzen > 90 & p.landfrac == 0 & p.rlat > latbins(iaTropics(jj)) & p.rlat <= latbins(iaTropics(jj+1))); 
        if length(iaX) > 3
          sergio.count(iCnt,jj) = length(iaX);

          tobs = rad2bt(h.vchan,p.robs1(:,iaX));
          sergio.mean_tobs(iCnt,jj,:) = mean(tobs');
          sergio.std_tobs(iCnt,jj,:) = std(tobs');
      
          tcal = rad2bt(h.vchan,p.rclr(:,iaX));
          sergio.mean_tcal(iCnt,jj,:) = mean(tcal');
          sergio.std_tcal(iCnt,jj,:) = std(tcal');
        else
          sergio.count(iCnt,jj) = 0;
          sergio.mean_tobs(iCnt,jj,:) = nan * ones(1,2645);
          sergio.std_tobs(iCnt,jj,:) = nan * ones(1,2645);
          sergio.mean_tcal(iCnt,jj,:) = nan * ones(1,2645);
          sergio.std_tcal(iCnt,jj,:) = nan * ones(1,2645);
        end  %% if iaX > 0
      end    %% for jj = 1 : length(iaTropics)
    end      %% if fname exist
  end        %% loop doy
end          %% loop yy
      
comment = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_April2019/run_quick_stats.m'
save /asl/s1/sergio/USEFUL_LARGE_MATFILES/Stats_Mean_StdDev/mean_stddev_clear_14year.mat   sergio comment
error('oo')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

load /asl/s1/sergio/USEFUL_LARGE_MATFILES/Stats_Mean_StdDev/mean_stddev_clear_14year.mat
yy=2003; doy=1;
fname = ['/asl/rtp/rtp_airicrad_v6/clear/' num2str(yy) '/era_airicrad_day' num2str(doy,'%03d') '_clear.rtp'];
[h,ha,p,pa] = rtpread(fname);

plot(squeeze(sergio.mean_tobs(:,1,1520)))
plot(squeeze(sergio.std_tobs(:,1,1520)))

oo = find(sergio.day_of_year == 1); errorbar(1:length(oo),squeeze(sergio.mean_tobs(oo,1,1520)),squeeze(sergio.std_tobs(oo,1,1520))); grid
errorbar(1:length(oo),squeeze(sergio.mean_tobs(oo,1,1520))-mean(squeeze(sergio.mean_tobs(oo,1,1520))),...
         squeeze(sergio.std_tobs(oo,1,1520)),'o-'); grid
signal = squeeze(sergio.mean_tobs(:,1,1520));


for ii = 1 : 365
  oo = find(sergio.day_of_year == ii);
  newsignal(oo) = signal(oo)-mean(signal(oo));
end
for ii = 366
  oo = find(sergio.day_of_year == 365 | sergio.day_of_year == 366);
  newsignal(oo) = signal(oo)-mean(signal(oo));
end
plot(smooth(newsignal,180)); grid

for ii = 1 : 365
  oo = find(sergio.day_of_year == ii);
  newsignal(oo) = signal(oo)-mean(signal(oo));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE
nedt = instr_chans('airs',2);
plot(sergio.count(:,1)); ylabel('count for tropical latbin');
plot(nedt)
plot(nedt); axis([1 2378 0 0.5])
plot(nedt(1291)./sqrt(sergio.count(:,1)))
plot(nedt(1291)./sqrt(sergio.count(:,1))); ylabel('NeDT for 1231 cm-1/sqrt(N)');
hist(nedt(1291)./sqrt(sergio.count(:,1))); ylabel('NeDT for 1231 cm-1/sqrt(N)');
hist(nedt(1291)./sqrt(sergio.count(:,1)),100); ylabel('NeDT for 1231 cm-1/sqrt(N)');

%function  y = nedt_T0_T1(nu,nedt0,BT0,BTX);
y = nedt_T0_T1(1231*ones(1,5427),nedt(1291),250*ones(1,5427),(squeeze(sergio.mean_tobs(:,1,1520)))');
plot(1:length(y),ones(size(y))*nedt(1291),1:length(y),y); 
  ylabel('NeDT (K)'); legend('250 K','obs temp');

plot(1:length(y),nedt(1291)./sqrt(sergio.count(:,1)),'b',1:length(y),y./sqrt(sergio.count(:,1)),'r');
semilogy(1:length(y),nedt(1291)./sqrt(sergio.count(:,1)),'b',1:length(y),y./sqrt(sergio.count(:,1)),'r');
  ylabel('NeDT/sqrt(N) (K)'); legend('250 K','obs temp');
hist(y./sqrt(sergio.count(:,1)),100); ylabel('NeDT for anomaly BT adjust 1231 cm-1/sqrt(N)');

%% now let us model the noise at 250 K
NedT_2645 = interp1(instr_chans,nedt,h.vchan);
for cc = 1 : 2645
  if mod(cc,500) == 1
    fprintf(1,'.');
  end
  yall(cc,:) = nedt_T0_T1(h.vchan(cc)*ones(1,5427),NedT_2645(cc),250*ones(1,5427),(squeeze(sergio.mean_tobs(:,1,cc)))');
end
fprintf(1,'\n');

%% we do 16 day averages so another factior of sqrt(16)
reduced_noise = nanmean(yall') ./ sqrt(max(sergio.count(:,1))) / sqrt(16);
wah = sergio.count(:,1); wah = wah(wah > 0); min(wah); reduced_noise = nanmean(yall') ./ sqrt(min(wah))  / sqrt(16);
reduced_noise = nanmean(yall') ./ sqrt(mean(sergio.count(:,1))) / sqrt(16);

plot(h.vchan,reduced_noise./NedT_2645','k'); axis([min(h.vchan) max(h.vchan) 0 0.1])
plot(h.vchan,NedT_2645'./reduced_noise,'k'); axis([min(h.vchan) max(h.vchan) 0 100]); title('Noise reduction (on average)')
  axis([min(h.vchan) 1620 0 60]); grid

semilogy(h.vchan,NedT_2645,'b',h.vchan,reduced_noise,'r',h.vchan,0.0025*ones(1,length(h.vchan)),'k'); 
axis([min(h.vchan) max(h.vchan) 0 1]); grid on
axis([min(h.vchan) 1620         0 1]); grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% from larrabee : daily
load noise_16day_avg_mission.mat
plot(mean(btn_av')); axis([2 40 0 0.01])
plot(mean(btn_av)); axis([1 2645 0 0.01])

load btn_avg.ma
whos
plot(1:362,squeeze(btn_avg(20,:,[62 449 1520 1700])),1:362,ones(1,362)*0.0025,'k'); legend(num2str([62 449 1520 1700]'))


