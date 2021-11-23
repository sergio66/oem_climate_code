%% based on driver_plot_T_WV_anomaly_retrievals.m

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP

latbins = equal_area_spherical_bands(21);
latbinsx = 0.5*(latbins(1:end-1)+latbins(2:end));
tropics = find(abs(latbinsx) <= 30);
tropics = 11 : 30;
extralats_tropics = find(abs(latbinsx) <= 50);   %% this is 6 - 37

iOBSorCAL = 0;
iAvgNumDays = 0;
dataset = 1; %% AIRS

iLatbin = 20;
if iOBSorCAL == 0 & dataset == 1
  fin = ['../../ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 1
  fin = ['../../ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
elseif iOBSorCAL == 0 & dataset == 2
  fin = ['../../IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 2
  fin = ['../../IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
end

load(fin)
save_days = avg_doy_since2002;
save_rtime = avg16_rtime;
save_dat_1231 = zeros(length(save_days),40);
save_days_map = zeros(length(save_days),40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('loading up the data')
if ~exist('jacmean')
  if ~exist('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/extralats_tropics_jacs.mat')
    disp('doing extra lats')
    for ii = 1 : 365
      fprintf(1,'ii = %3i \n',ii)
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' num2str(ii,'%03d') '_M_TS_jac_all_5_97_97_97_2645.mat'];
      load(fname)
      %jacall(ii,:,:,:) = M_TS_jac_all(extralats_tropics,:,:);
      jacmean(ii,:,:)  = squeeze(nanmean(M_TS_jac_all(extralats_tropics,:,:),1));
    end
    save -v7.3 /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/extralats_tropics_jacs.mat f qrenorm jac*
  else
    load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/extralats_tropics_jacs.mat
  end
end

disp('loading up the jacs')
if ~exist('jacmeanTrop')
  if ~exist('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_jacs.mat')
    disp('doing trops only')
    for ii = 1 : 365
      fprintf(1,'ii = %3i \n',ii)
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' num2str(ii,'%03d') '_M_TS_jac_all_5_97_97_97_2645.mat'];
      load(fname)
      %jacall(ii,:,:,:) = M_TS_jac_all(tropics,:,:);
      jacmeanTrop(ii,:,:)  = squeeze(nanmean(M_TS_jac_all(tropics,:,:),1));
    end
    save -v7.3 /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_jacs.mat f qrenorm jacmeanTrop
  else
    load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/only_tropics_jacs.mat
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('loading up the retrievals')   %% based on compare_anomaly_runs
%file1 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/nopert_jac_anomaly_0dayavg_results';

%% based on running driver_plot_T_WV_anomaly_retrievals for the "nopert" runs then
%% save NewStuff_forAMT2020_Response2Reviewers/PartitionAnomSignal_by_Constituent/complete_jac_anomaly_0dayavg_results co2 ch4 n2o stemp cfc* tz wv o3
file1 = 'NewStuff_forAMT2020_Response2Reviewers/PartitionAnomSignal_by_Constituent/complete_jac_anomaly_0dayavg_results';

results = load(file1);

iWhich = -1;       
iWhich = +1;       
if iWhich == -1
  disp('using mean jac over extra tropical lats')
  jacX = jacmean;     %% extra trops and trops
elseif iWhich == +1
  disp('using mean jac over only tropical lats')
  jacX = jacmeanTrop; %% trop only
end

co2 = (nanmean(results.co2(extralats_tropics,:),1)' * ones(1,2645)) .* squeeze(jacX(:,:,1)) / qrenorm(1);
n2o = (nanmean(results.n2o(extralats_tropics,:),1)' * ones(1,2645)) .* squeeze(jacX(:,:,2)) / qrenorm(2);
ch4 = (nanmean(results.ch4(extralats_tropics,:),1)' * ones(1,2645)) .* squeeze(jacX(:,:,3)) / qrenorm(3);
cfc12 = (nanmean(results.cfc12(extralats_tropics,:),1)' * ones(1,2645)) .* squeeze(jacX(:,:,5)) / qrenorm(5);
stemp = (nanmean(results.stemp(extralats_tropics,:),1)' * ones(1,2645)) .* squeeze(jacX(:,:,6)) / qrenorm(6);

wvjacmean = zeros(365,2645,20);
tzjacmean  = zeros(365,2645,20);
o3jacmean = zeros(365,2645,20);
for ii = 1 : 20
  junk = results.a.jacobian.wvjaclays_used{ii};
  wvjacmean(:,:,ii) = sum(squeeze(jacX(:,:,junk+0*97)),3);
  tzjacmean(:,:,ii) = sum(squeeze(jacX(:,:,junk+1*97)),3);
  o3jacmean(:,:,ii) = sum(squeeze(jacX(:,:,junk+2*97)),3);
end
wv = zeros(365,2645);
tz  = zeros(365,2645);
o3 = zeros(365,2645);
for ii = 1 : 20
  junk = (nanmean(results.wv(extralats_tropics,:,ii),1)' * ones(1,2645)) .* squeeze(wvjacmean(:,:,ii)) / qrenorm(100);
    wv = wv + junk;
  junk = (nanmean(results.tz(extralats_tropics,:,ii),1)' * ones(1,2645)) .* squeeze(tzjacmean(:,:,ii)) / qrenorm(200);
    tz = tz + junk;
  junk = (nanmean(results.o3(extralats_tropics,:,ii),1)' * ones(1,2645)) .* squeeze(o3jacmean(:,:,ii)) / qrenorm(297);
    o3 = o3 + junk;
end

figure(1); pcolor(f,2002+save_days/365,co2); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('co2')
figure(2); pcolor(f,2002+save_days/365,n2o); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('n2o')
figure(3); pcolor(f,2002+save_days/365,ch4); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('ch4')
figure(4); pcolor(f,2002+save_days/365,stemp); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('stemp')

figure(5); pcolor(f,2002+save_days/365,wv); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('wv')
figure(6); pcolor(f,2002+save_days/365,tz); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('t')
figure(7); pcolor(f,2002+save_days/365,o3); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('o3')

total = co2 + n2o + ch4 + stemp + wv + tz + o3;
figure(8); pcolor(f,2002+save_days/365,total); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('retrieval total')
cxA = caxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spectral_results = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/anomaly_0dayavg_results_spectra.mat');  %% this is 11 -- 30
figure(9); pcolor(f,2002+save_days/365,spectral_results.raaObs'); shading flat; colormap jet; colorbar; cx = caxis; cxx = max(abs(cx)); cx = [-cxx +cxx]; caxis(cx); colormap(usa2)
  title('signal')
cxB = caxis;

figure(8); caxis(cxA);
figure(9); caxis(cxA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meansignal = nanmean(spectral_results.raaObs');
meanfit    = nanmean(spectral_results.raaCal');
meanco2    = nanmean(co2);
meann2o    = nanmean(n2o);
meanch4    = nanmean(ch4);
meancfc12  = nanmean(cfc12);
meanstemp  = nanmean(stemp);
meanwv     = nanmean(wv);
meantz     = nanmean(tz);
meano3     = nanmean(o3);

figure(10); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,meanfit,'color',[0.6 0.6 0.6],'linewidth',2); junk = [meanco2; meann2o; meanch4; meanstemp; meanwv; meantz; meano3]; plot(f,sum(junk),'r'); hold off
hl = legend('signal','fit','sum(trace+T+ST+gas','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])
plotaxis2;

figure(11); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,meanfit,'color',[0.6 0.6 0.6],'linewidth',2); hold on; plot(f,[meanco2; meann2o; meanch4; meanstemp; meanwv; meantz; meano3]); hold off
hl = legend('signal','fit','co2','n2o','ch4','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])
plotaxis2;

figure(11); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,[meanco2; meann2o; meanch4; meanstemp; meanwv; meantz; meano3]); hold off
hl = legend('signal','co2','n2o','ch4','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])
plotaxis2;

chanset = results.a.jacobian.chanset;
co2chan = find(f(chanset) >= 722 & f(chanset) <= 800 & (meanco2(chanset))' < -0.6,1); co2chan = chanset(co2chan); fprintf(1,'co2chan = %4i f(co2chan) = %8.6f \n',co2chan,f(co2chan));
n2ochan = find(f(chanset) >= 1280,1);                                                 n2ochan = chanset(n2ochan); fprintf(1,'n2ochan = %4i f(n2ochan) = %8.6f \n',n2ochan,f(n2ochan));
n2ochan = find(f >= 1280,1);                                                                                      fprintf(1,'n2ochan = %4i f(n2ochan) = %8.6f \n',n2ochan,f(n2ochan));
ch4chan = find(f(chanset) >= 1304,1);                                                 ch4chan = chanset(ch4chan); fprintf(1,'ch4chan = %4i f(ch4chan) = %8.6f \n',ch4chan,f(ch4chan));
cfc12chan = find(f(chanset) >= 923,1);                                                cfc12chan = chanset(cfc12chan); fprintf(1,'cfc12chan = %4i f(cfc12chan) = %8.6f \n',cfc12chan,f(cfc12chan));
stempchan = find(f(chanset) >= 900,1);                                                stempchan = chanset(stempchan); fprintf(1,'stempchan = %4i f(stempchan) = %8.6f \n',stempchan,f(stempchan));
wvchan = find(f(chanset) >= 1418,1);                                                  wvchan = chanset(wvchan); fprintf(1,'wvchan = %4i f(wvchan) = %8.6f \n',wvchan,f(wvchan));
tzchan = find(f(chanset) >= 1418,1);   tzchan = co2chan;                                                        fprintf(1,'tzchan = %4i f(tzchan) = %8.6f \n',tzchan,f(tzchan));
o3chan = find(f(chanset) >= 1044,1);                                                  o3chan = chanset(o3chan); fprintf(1,'o3chan = %4i f(o3chan) = %8.6f \n',o3chan,f(o3chan));

figure(11); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,[meanco2; meann2o; meanch4; meancfc12; meanstemp; meanwv; meantz; meano3]); hold off
chans = [co2chan n2ochan ch4chan cfc12chan stempchan wvchan tzchan o3chan];
chans = [co2chan n2ochan ch4chan cfc12chan stempchan wvchan        o3chan];
hold on; 
  plot(f(co2chan),meanco2(co2chan),'rx','linewidth',2);
  plot(f(n2ochan),meann2o(n2ochan),'rx','linewidth',2);
  plot(f(ch4chan),meanch4(ch4chan),'rx','linewidth',2);
  plot(f(stempchan),meanstemp(stempchan),'rx','linewidth',2);
  plot(f(wvchan),meanwv(wvchan),'rx','linewidth',2);
  plot(f(tzchan),meantz(tzchan),'rx','linewidth',2);
  plot(f(o3chan),meano3(o3chan),'rx','linewidth',2);
hold off
hl = legend('signal','co2','n2o','ch4','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])

figure(12);
plot(2002+save_days/365,spectral_results.raaObs(co2chan,:),'b',2002+save_days/365,co2(:,co2chan),'b--','linewidth',2);
hold on
plot(2002+save_days/365,spectral_results.raaObs(n2ochan,:),'r',2002+save_days/365,n2o(:,n2ochan),'r--','linewidth',2);
hold off

figure(12)
plot(2002+save_days/365,co2(:,co2chan),'b',2002+save_days/365,n2o(:,n2ochan),'g',2002+save_days/365,ch4(:,ch4chan),'r',2002+save_days/365,cfc12(:,ch4chan),'b--',...
     2002+save_days/365,stemp(:,stempchan),'c',2002+save_days/365,wv(:,wvchan),'m',2002+save_days/365,tz(:,tzchan),'k',2002+save_days/365,o3(:,o3chan),'y')
hl = legend(['co2 ' num2str(f(co2chan)) 'cm-1'],['n2o ' num2str(f(n2ochan)) 'cm-1'],['ch4 ' num2str(f(ch4chan)) 'cm-1'],['cfc12 ' num2str(f(ch4chan)) 'cm-1'],...
            ['stemp ' num2str(f(stempchan)) 'cm-1'],['wv ' num2str(f(wvchan)) 'cm-1'],['tz ' num2str(f(tzchan)) 'cm-1'],['o3 ' num2str(f(o3chan)) 'cm-1'],'location','best');
set(hl,'fontsize',8);
grid; ax = axis; ax(1) = min(2002+save_days/365); ax(2) = max(2002+save_days/365); axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /asl/matlib/plotutils

iPrint = -1;
figure(11); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,[meanco2; meann2o; meanch4; meancfc12; meanstemp; meanwv; meantz; meano3]); hold off
hl = legend('signal','co2','n2o','ch4','cfc12','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])
if iPrint > 0
  aslprint('spectral_without_x.png')
end

figure(11); 
plot(f,meansignal,'k.-','linewidth',2); hold on; plot(f,[meanco2; meann2o; meanch4; meancfc12; meanstemp; meanwv; meantz; meano3]); hold off
chans = [co2chan n2ochan ch4chan stempchan wvchan tzchan o3chan];
chans = [co2chan n2ochan ch4chan stempchan wvchan        o3chan];
hold on; 
  plot(f(co2chan),meanco2(co2chan),'rx','linewidth',2);
  plot(f(n2ochan),meann2o(n2ochan),'rx','linewidth',2);
  plot(f(ch4chan),meanch4(ch4chan),'rx','linewidth',2);
  plot(f(stempchan),meanstemp(stempchan),'rx','linewidth',2);
  plot(f(wvchan),meanwv(wvchan),'rx','linewidth',2);
  plot(f(tzchan),meantz(tzchan),'rx','linewidth',2);
  plot(f(o3chan),meano3(o3chan),'rx','linewidth',2);
hold off
hl = legend('signal','co2','n2o','ch4','cfc12','stemp','wv','tz','o3','location','best'); set(hl,'fontsize',10); grid
axis([640 1640 -0.8 +0.4])
if iPrint > 0
  aslprint('spectral_with_x.png')
end

figure(12)
plot(2002+save_days/365,co2(:,co2chan),'b',2002+save_days/365,n2o(:,n2ochan),'g',2002+save_days/365,ch4(:,ch4chan),'r',...
     2002+save_days/365,cfc12(:,cfc12chan),'b--',2002+save_days/365,stemp(:,stempchan),'c',2002+save_days/365,wv(:,wvchan),'m',2002+save_days/365,tz(:,tzchan),'k',2002+save_days/365,o3(:,o3chan),'y')
hl = legend(['co2 ' num2str(f(co2chan)) 'cm-1'],['n2o ' num2str(f(n2ochan)) 'cm-1'],['ch4 ' num2str(f(ch4chan)) 'cm-1'],['cfc12 ' num2str(f(cfc12chan)) 'cm-1'],...
            ['stemp ' num2str(f(stempchan)) 'cm-1'],['wv ' num2str(f(wvchan)) 'cm-1'],['tz ' num2str(f(tzchan)) 'cm-1'],['o3 ' num2str(f(o3chan)) 'cm-1'],'location','best');
set(hl,'fontsize',8);
grid; ax = axis; ax(1) = min(2002+save_days/365); ax(2) = max(2002+save_days/365); axis(ax);
if iPrint > 0
  %aslprint('individual_chans_marked_x_rates.png')
  aslprint('individual_chans_marked_x_rates_v2.pdf')
end
